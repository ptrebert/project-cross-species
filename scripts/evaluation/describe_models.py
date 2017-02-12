#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as js
import csv as csv

import pandas as pd


# Note to self: model_query is vital for access, since DataFrame.query is a method of the object
__model_columns__ = ['groupid', 'group', 'model_name', 'model_task', 'model_type', 'n_estimators', 'min_samples_split',
                     'model_target', 'model_query', 'eid', 'epi_project', 'tid', 'trans_project', 'cell',
                     'num_hist', 'num_dnase', 'cv_scoring', 'cv_perf', 'model_file']


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--datasets', '-d', type=str,
                        default='/home/pebert/work/code/mpggit/crossspecies/annotation/exec/datasets.tsv',
                        dest='datasetfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def load_dataset_annotation(fpath):
    """
    :param fpath:
    :return:
    """
    datasets = dict()
    with open(fpath, 'r') as ann:
        rows = csv.DictReader(ann, delimiter='\t')
        for r in rows:
            datasets[r['id']] = r
    return datasets


def read_model_information(metadata, datasets, index):
    """
    :param metadata:
    :return:
    """
    infos = metadata['model_info']
    row = {'model_name': infos['name'], 'n_estimators': infos['params']['n_estimators'],
           'min_samples_split': infos['params']['min_samples_split'], 'model_task': infos['type']}
    train = metadata['training_info']
    row['cv_perf'] = train['best_score']
    row['cv_scoring'] = train['scoring']
    run = metadata['run_info']
    modelparts = run['model_file'].split('.')
    modeltype = modelparts[4]
    sampleparts = modelparts[0].split('_')
    trg, qry, epi, trans, cell = sampleparts[3], modelparts[2], sampleparts[1], sampleparts[2], sampleparts[4]
    epi_proj, trans_proj = datasets[epi]['project'], datasets[trans]['project']
    if modeltype == 'seq':
        group, groupid, hist, dnase = sampleparts[0][:3], sampleparts[0], 0, 0
    else:
        group, groupid, hist, dnase = sampleparts[0][:3], sampleparts[0], int(sampleparts[0][3]), int(sampleparts[0][4])
    tmp = {'model_file': run['model_file'], 'model_type': modeltype, 'model_target': trg, 'model_query': qry,
           'group': group, 'groupid': groupid, 'num_hist': hist, 'num_dnase': dnase, 'eid': epi, 'tid': trans,
           'epi_project': epi_proj, 'trans_project': trans_proj, 'cell': cell}
    row.update(tmp)
    check = [(k, k in __model_columns__) for k in row.keys()]
    assert all([t[1] for t in check]), 'Missing information: {}'.format(check)
    df_row = pd.DataFrame.from_records([row], index=[index])
    #row.index = index
    return df_row


def read_feature_information(metadata, index):
    """
    :param metadata:
    :param index:
    :return:
    """
    imps = metadata['attribute_info']['feature_importances_']
    names = metadata['feature_info']['order']
    row = pd.DataFrame([imps], index=[index], columns=names, dtype='float64')
    return row


def init_dataframes(fpath):
    """
    :param fpath:
    :return:
    """

    with open(fpath, 'r') as mdfile:
        md = js.load(mdfile)
        feat_names = sorted(md['feature_info']['order'])
        feat_desc = pd.DataFrame([], index=[], columns=feat_names, dtype='float64')
        model_desc = pd.DataFrame([], index=[], columns=__model_columns__, dtype='object')
    return model_desc, feat_desc


def main():
    """
    :return:
    """
    args = parse_command_line()
    model_desc = None
    feat_desc = None
    datasets = load_dataset_annotation(args.datasetfile)
    for idx, fpath in enumerate(sorted(args.inputfiles)):
        load_path = fpath.replace('.pck', '.json')
        if model_desc is None:
            model_desc, feat_desc = init_dataframes(load_path)
        with open(load_path, 'r') as mdf:
            metadata = js.load(mdf)
            this_model = read_model_information(metadata, datasets, idx)
            model_desc = pd.concat([model_desc, this_model], ignore_index=False)
            this_feat = read_feature_information(metadata, idx)
            feat_desc = pd.concat([feat_desc, this_feat], ignore_index=False)

    with pd.HDFStore(args.outputfile, 'w', complevel=9, complib='blosc') as hdf:
        hdf.put('/models', model_desc, format='fixed')
        hdf.put('/features', feat_desc, format='fixed')
        hdf.flush()
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
