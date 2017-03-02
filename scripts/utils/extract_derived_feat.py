#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as js
import collections as col
import statistics as st
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--run-md-file', '-tf', type=str, nargs='+', dest='runmdfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--feat-name', '-ftn', type=str, dest='featname')
    args = parser.parse_args()
    return args


def extract_class_probabilities(feat, md, training):
    """
    :param feat:
    :param md:
    :param training:
    :return:
    """
    sample_names = md['sample_info']['names']
    if training:
        class_probs = md['training_info']['true_class_prob']
        class_label = md['sample_info']['targets']
    else:
        class_probs = md['testing_info']['probabilities']['pred']
        class_label = md['testing_info']['targets']['pred']
    for n, p, l in zip(sample_names, class_probs, class_label):
        if l == 1:
            feat[n].append(p * 100.)
        else:
            feat[n].append((1 - p) * 100.)
    return feat


def collect_values(mdfile, feat):
    """
    :param mdfile:
    :param feat:
    :return:
    """
    mtype = None
    for mdf in mdfile:
        if mdf.endswith('.pck'):
            load_file = mdf.replace('.pck', '.json')
        else:
            load_file = mdf
        with open(load_file, 'r') as dumped:
            md = js.load(dumped)
            this_model = md['model_info']['type']
            assert this_model == mtype or mtype is None, 'Model mismatch: {} vs {}'.format(mtype, this_model)
            mtype = this_model
            if mtype == 'classifier':
                feat = extract_class_probabilities(feat, md, 'training_info' in md)
            elif mtype == 'regressor':
                raise NotImplementedError
            else:
                raise ValueError('Unkown model type: {}'.format(mtype))
    return feat


def main():
    """
    :return:
    """
    args = parse_command_line()
    feat = col.defaultdict(list)
    feat = collect_values(args.runmdfile, feat)
    assert feat, 'No features derived from input files'
    feat = [(k, st.mean(v)) for k, v in feat.items()]
    df = pd.DataFrame(feat, columns=['name', args.featname])
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('ftdrv', df, format='fixed')
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
