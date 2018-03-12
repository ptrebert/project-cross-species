#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as js
import numpy as np
import pandas as pd
import collections as col
import gzip as gz
from sklearn.metrics import r2_score


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='inputfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--model-query', '-mq', type=str, dest='query')
    parser.add_argument('--model-type', '-mt', type=str,choices=['full', 'sig', 'seq'], dest='mtype')
    parser.add_argument('--test-type', '-tt', type=str, dest='test')
    parser.add_argument('--spec-match', '-sm', type=str, dest='spec')
    parser.add_argument('--setting', '-s', type=int, choices=[1, 2, 3], dest='setting')
    parser.add_argument('--load-genes', '-g', type=str, dest='genes')
    parser.add_argument('--load-column', '-col', type=int, dest='column', default=1)
    parser.add_argument('--missing-as-zero', action='store_true', default=False, dest='nazero')
    args = parser.parse_args()
    assert os.path.isfile(args.inputfile), 'Invalid path to input file: {}'.format(args.inputfile)
    return args


def init_dataframe(args):
    """
    :param args:
    :return:
    """
    col_idx = args.column
    gene_ids = []
    with gz.open(args.genes, 'rt') as genefile:
        for line in genefile:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.strip().split()
            gene_ids.append(cols[col_idx])
    gene_ids = sorted(set(gene_ids))
    df = pd.DataFrame([], index=[], columns=gene_ids, dtype='float64')
    return df


def load_class_scores(perfdata, sample_index):
    """
    :param perfdata:
    :return:
    """
    lolim, hilim = 0.4, 0.6
    sample_names = perfdata['sample_info']['names']
    # list of lists (2-tuple)
    class_probs = perfdata['testing_info']['probabilities']['all']
    assert perfdata['testing_info']['targets']['order'] == sorted(perfdata['testing_info']['targets']['order']),\
        'Handling unsorted class labels not supported'
    true_labels = perfdata['testing_info']['targets']['true']
    pred_labels = perfdata['testing_info']['targets']['pred']
    perf_score = perfdata['testing_info']['performance']
    assert len(true_labels) == len(class_probs) == len(sample_names),\
        'Missing information: {} vs {} vs {}'.format(len(true_labels), len(class_probs), len(sample_names))
    counts = []
    true_preds = set()
    for truelab, predlab, prob, name in zip(true_labels, pred_labels, class_probs, sample_names):
        if lolim < prob[truelab] < hilim:
            counts.append(0)
        elif truelab != predlab:
            counts.append(-1)
        else:
            counts.append(1)
            true_preds.add(name)
    sample = pd.DataFrame([counts], index=[sample_index], columns=sample_names, dtype='float64')
    return sample, true_preds, perf_score


def load_reg_scores(perfdata, sample_index):
    """
    :param perfdata:
    :param sample_index:
    :return:
    """
    sample_names = perfdata['sample_info']['names']
    true_values = np.expm1(perfdata['testing_info']['targets']['true'])
    pred_values = np.expm1(perfdata['testing_info']['targets']['pred'])
    r2_realspace = r2_score(true_values, pred_values)
    assert len(true_values) == len(pred_values) == len(sample_names),\
        'Missing information: {} vs {} vs {}'.format(len(true_values), len(pred_values), len(sample_names))
    scores = []
    true_est = set()
    for trueval, predval, name in zip(true_values, pred_values, sample_names):
        true_int = round(trueval)
        pred_int = round(predval)
        lo, hi = round(trueval * 0.8), round(trueval * 1.2)
        if lo <= round(predval) <= hi:
            true_est.add(name)
        res = true_int - pred_int
        scores.append(res)
    sample = pd.DataFrame([scores], index=[sample_index], columns=sample_names, dtype='float64')
    return sample, true_est, r2_realspace


def select_data_loader(setnum):
    """
    :param setnum:
    :return:
    """
    options = {1: load_class_scores, 2: load_reg_scores, 3: load_reg_scores}
    return options[setnum]


def load_run_data(fpath, loader, index, df):
    """
    :param fpath:
    :param df:
    :return:
    """
    with open(fpath, 'r') as rundata:
        data = js.load(rundata)
        scores, true_est, perf = loader(data, index)
    df = pd.concat([df, scores], axis=0)
    return df, true_est, perf, data['run_info']['data_file'], data['run_info']['model_file']


def main():
    """
    :return:
    """
    args = parse_command_line()
    select_pair = args.spec, args.test
    df = init_dataframe(args)
    data_loader = select_data_loader(args.setting)
    count_cons = col.Counter()
    run_counts = 0
    with open(args.inputfile, 'r') as md:
        metadata = js.load(md)
        assert_data_exists(metadata['run_pairs'], select_pair)
        out_buffer = []
        for groupid, runs in metadata['test_runs'].items():
            for run in runs:
                if (args.setting == run['setting_number'] and
                        select_pair == (run['run_spec_match'], run['run_test_type']) and
                        args.query == run['model_metadata']['model_query'] and
                        args.mtype in run['model_metadata']['model_spec']):
                    sample_index = groupid + '_' + str(args.setting) + str(run['run_number'])
                    df, true_est, perf, datafile, modelfile = load_run_data(run['run_file'], data_loader, sample_index, df)
                    out_buffer.append((datafile, modelfile, np.round(perf, 4)))
                    count_cons.update(true_est)
                    run_counts += 1
        out_buffer = sorted(out_buffer)
        for item in out_buffer:
            print(item[0], ' - ', item[1], ':  ', item[2])
    cons_threshold = round(run_counts * 0.2)
    df.dropna(axis=1, how='all', inplace=True)
    if args.nazero:
        df.fillna(value=0., inplace=True)
    total_sum = df.sum(axis=0)
    assert len(total_sum) == df.shape[1], 'Wrong dimensionality: {} vs {}'.format(len(total_sum), df.shape[1])
    #avg_per_run = df.mean(axis=0)
    #assert len(avg_per_run) == df.shape[1], 'wrong dimension: {} vs {}'.format(len(avg_per_run), df.shape[1])
    skipped = 0
    # with open(args.outputfile.replace('.h5', '.tsv'), 'w') as out:
    #     for a, n in zip(total_sum, df.columns):
    #         if np.isnan(a):
    #             skipped += 1
    #             continue
    #         out_a = str(int(a))
    #         out_n = n.split('.')[0]
    #         _ = out.write(out_n + '\t' + out_a + '\n')
    with open(args.outputfile.replace('.h5', '.tsv'), 'w') as out:
        mc = count_cons.most_common()
        for name, count in mc:
            if count < cons_threshold:
                break
            out_n = name.split('.')[0]
            #_ = out.write(out_n + '\t' + str(count) + '\n')
            _ = out.write(out_n + '\n')
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('/counts', df)
    return


def assert_data_exists(runpairs, selectpair):
    """
    :param runpairs:
    :param selectpair:
    :return:
    """
    assert any([selectpair == tuple(t[0]) for t in runpairs]), \
        'Specified pairing does not exist in dataset: {}'.format(selectpair)
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
