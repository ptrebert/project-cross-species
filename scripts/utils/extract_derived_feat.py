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
import functools as fnt


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--run-md-file', '-tf', type=str, nargs='+', dest='runmdfile', required=True)
    parser.add_argument('--output', '-o', type=str, dest='outputfile', required=True)
    parser.add_argument('--feat-name', '-ftn', type=str, dest='featname', required=True)
    parser.add_argument('--label-name', '-lbn', type=str, dest='labelname', default='class_label')
    args = parser.parse_args()
    return args


def extract_class_probabilities(feat, labels, md, training):
    """
    :param feat:
    :param md:
    :param training:
    :return:
    """
    sample_names = md['sample_info']['names']
    if training:
        class_probs = md['attribute_info']['oob_decision_function_']
        class_label = md['sample_info']['targets']
        for n, (p0, p1), l in zip(sample_names, class_probs, class_label):
            # the derived feature in the classification scenario gives
            # the probability that this gene is active
            if l == 1:
                feat[n].append(p1 * 100.)
            else:
                feat[n].append((1 - p1) * 100.)
            # the oob predicted label is recorded for later (easier) subset
            # selection of the training data set
            pred_label = 1 if p1 > p0 else 0
            labels[n].append(pred_label)
    else:
        class_probs = md['testing_info']['probabilities']['pred']
        class_label = md['testing_info']['targets']['pred']
        for n, p, l in zip(sample_names, class_probs, class_label):
            if l == 1:
                feat[n].append(p * 100.)
            else:
                feat[n].append((1 - p) * 100.)
            labels[n].append(l)
    return feat, labels


def extract_regression_estimates(feat, md, training):
    """
    :param feat:
    :param md:
    :param training:
    :return:
    """
    sample_names = md['sample_info']['names']
    if training:
        est_value = md['attribute_info']['oob_prediction_']
        for n, v in zip(sample_names, est_value):
            feat[n].append(v)
    else:
        est_value = md['testing_info']['targets']['pred']
        for n, v in zip(sample_names, est_value):
            feat[n].append(v)
    return feat


def merge_feat_labels(features, labels, featname, labelname):
    """
    :param features:
    :param labels:
    :return:
    """
    records = []
    for k, labs in labels.items():
        feat = features[k]
        if len(set(labs)) == 1:
            records.append([k, st.mean(feat), labs.pop(0)])
        else:
            select = sorted([(f, l) for f, l in zip(feat, labs)], reverse=True)
            records.append([k, st.mean(feat), select[0][1]])
    df = pd.DataFrame(records, columns=['name', featname, labelname])
    assert not df.empty, 'Created empty dataframe when merging features/labels'
    return df


def process_est(vtype, values):
    """
    :param vtype:
    :param values:
    :return:
    """
    if vtype.startswith('int'):
        est = int(round(st.median(values), ndigits=0))
    elif vtype.startswith('float'):
        est = float(st.mean(values))
    else:
        raise ValueError('Unexpected type of target variable: {}'.format(vtype))
    return est


def merge_feat(features, featname, vartype):
    """
    :param features:
    :param featname:
    :return:
    """
    records = []
    make_est = fnt.partial(process_est, *(vartype, ))
    for k, vals in features.items():
        records.append([k, make_est(vals)])
    df = pd.DataFrame(records, columns=['name', featname])
    assert not df.empty, 'Created empty dataframe when merging features'
    return df


def collect_values(mdfile, featname, labelname):
    """
    :param mdfile:
    :param featname:
    :param labelname:
    :return:
    """
    feat = col.defaultdict(list)
    labels = col.defaultdict(list)
    mtype = None
    df = None
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
                feat, labels = extract_class_probabilities(feat, labels, md, 'training_info' in md)
                df = merge_feat_labels(feat, labels, featname, labelname)
            elif mtype == 'regressor':
                feat = extract_regression_estimates(feat, md, 'training_info' in md)
                df = merge_feat(feat, featname, md['dataset_info']['target_type'])
            else:
                raise ValueError('Unkown model type: {}'.format(mtype))
    return df


def main():
    """
    :return:
    """
    args = parse_command_line()
    df = collect_values(args.runmdfile, args.featname, args.labelname)
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
