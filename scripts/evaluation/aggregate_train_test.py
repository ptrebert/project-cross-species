#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import fnmatch as fnm
import json as js
import csv as csv


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--train-root', '-tr', type=str, dest='trainroot')
    parser.add_argument('--test-root', '-ts', type=str, dest='testroot')
    parser.add_argument('--target', '-trg', type=str, dest='target')

    parser.add_argument('--output', '-o', type=str, dest='outputfile')

    args = parser.parse_args()
    return args


def extract_training_infos(md):
    """
    :param md:
    :return:
    """
    train_infos = dict()
    train_infos['num_features'] = md['dataset_info']['num_features']
    train_infos['num_samples'] = md['dataset_info']['num_samples']
    rule = md['dataset_info']['derive_target']
    if rule is None or not rule:
        train_infos['target_rule'] = 'raw'
    else:
        train_infos['target_rule'] = '(' + rule + ')'
    train_infos['target_var'] = md['dataset_info']['target_var']
    train_infos['score'] = md['training_info']['best_score']
    subset = md['dataset_info']['subset']
    if subset is None or not subset:
        train_infos['subset'] = 'full'
    else:
        train_infos['subset'] = '(' + subset + ')'
    train_infos['scoring'] = md['training_info']['scoring']
    train_infos['model_type'] = md['model_info']['type']
    train_infos['model_file'] = md['run_info']['model_file']
    return train_infos


def extract_testing_infos(md):
    """
    :param md:
    :return:
    """
    test_infos = dict()
    test_infos['num_features'] = md['dataset_info']['num_features']
    test_infos['num_samples'] = md['dataset_info']['num_samples']
    rule = md['dataset_info']['derive_target']
    if rule is None or not rule:
        test_infos['target_rule'] = 'raw'
    else:
        test_infos['target_rule'] = '(' + rule + ')'
    test_infos['target_var'] = md['dataset_info']['target_var']
    test_infos['score'] = md['testing_info']['performance']
    test_infos['scoring'] = md['testing_info']['scoring']
    subset = md['dataset_info']['subset']
    if subset is None or not subset:
        test_infos['subset'] = 'full'
    else:
        test_infos['subset'] = '(' + subset + ')'
    test_infos['model_type'] = md['model_info']['type']
    test_infos['model_file'] = md['run_info']['model_file']
    return test_infos


def collect_training_stats(rootpath, target):
    """
    :param rootpath:
    :return:
    """
    run_infos = []
    for root, dirs, files in os.walk(rootpath):
        if files:
            subdir = os.path.split(root)[1]
            if not subdir.startswith(target):
                continue
            trg, qry = subdir.split('_to_')
            assert trg == target, 'File walker is lost and confused: {} and {}'.format(root, subdir)
            models = fnm.filter(files, '*.pck')
            for mf in models:
                groupid = mf.split('_')[0]
                # manual fix due to lacking kidney data in hg19
                if (trg == 'hg19' or qry == 'hg19') and groupid.startswith('G09'):
                    continue
                num_hist = groupid[3]
                num_dnase = groupid[4]
                mdfile = os.path.join(root, mf.replace('.pck', '.json'))
                with open(mdfile, 'r') as dump:
                    md = js.load(dump)
                    run = {'target': trg, 'query': qry, 'groupid': groupid,
                           'num_histone': num_hist, 'num_dnase': num_dnase,
                           'setting': 'training'}
                    traininfos = extract_training_infos(md)
                    run.update(traininfos)
                    run_infos.append(run)
    return run_infos


def collect_testing_stats(testroot, target):
    """
    :param testroot:
    :param target:
    :return:
    """
    run_infos = []
    match_target = '*_trg-{}_qry*.json'.format(target)
    for root, dirs, files in os.walk(testroot):
        if files:
            outputs = fnm.filter(files, match_target)
            for o in outputs:
                groupid, trg, qry, _ = o.split('_', 3)
                trg = trg.replace('trg-', '')
                qry = qry.replace('qry-', '')
                assert trg == target, 'File walker is lost and confused: {} and {}'.format(root, o)
                # same fix as for training data
                if (trg == 'hg19' or qry == 'hg19') and groupid.startswith('G09'):
                    continue
                num_hist = groupid[3]
                num_dnase = groupid[4]
                with open(os.path.join(root, o), 'r') as dump:
                    md = js.load(dump)
                    run = {'target': trg, 'query': qry, 'groupid': groupid,
                           'setting': 'testing', 'num_histone': num_hist,
                           'num_dnase': num_dnase}
                    testinfos = extract_testing_infos(md)
                    run.update(testinfos)
                    run_infos.append(run)
    return run_infos


def main():
    """
    :return:
    """
    args = parse_command_line()
    train_runs = collect_training_stats(args.trainroot, args.target)
    test_runs = collect_testing_stats(args.testroot, args.target)

    fields = sorted(train_runs[0].keys())
    with open(args.outputfile, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=fields,
                                delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted(train_runs, key=lambda d: d['groupid']))
        writer.writerows(sorted(test_runs, key=lambda d: d['groupid']))

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
