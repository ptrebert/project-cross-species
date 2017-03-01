#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as js
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--train-file', '-tf', type=str, dest='trainfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    args = parse_command_line()
    load_file = args.trainfile
    if load_file.endswith('.pck'):
        load_file = load_file.replace('.pck', '.json')
    class_priors = []
    with open(load_file, 'r') as dumped:
        train_infos = js.load(dumped)
        sample_names = train_infos['sample_info']['names']
        true_class_probs = train_infos['training_info']['true_class_prob']
        true_class = train_infos['sample_info']['targets']
    for n, p, c in zip(sample_names, true_class_probs, true_class):
        # probabilities * 100:
        # expressed in percent just by convention of other feature definitions
        if c == 1:
            class_priors.append([n, p * 100])
        else:
            class_priors.append([n, (1 - p) * 100])
    df = pd.DataFrame(class_priors, columns=['name', 'ftprior_pct_active'])

    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('prior', df, format='fixed')

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
