#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import fnmatch as fnm
import json as js
import collections as col
import csv as csv

import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support,\
    accuracy_score, roc_auc_score, roc_curve

__DATASET_FILE__ = '/home/pebert/work/code/mpggit/crossspecies/annotation/exec/datasets.tsv'
__MATCHTYPES__ = '/home/pebert/work/code/mpggit/crossspecies/annotation/exec/cellmatches_ro.json'
__APPLY_ROOT__ = '/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/processing/norm/task_applymodel_exp/sub_status'
__ORTHOLOGS__ = '/TL/deep/fhgfs/projects/pebert/thesis/refdata/orthologs/hdf/odb9_6species.h5'
__AGG_EXPRESSION__ = '/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/rawdata/conv/agg'
__CONSERVATION__ = '/TL/deep/fhgfs/projects/pebert/thesis/refdata/conservation/genes'
__TAB_ASSEMBLY__ = '/home/pebert/work/code/mpggit/refdata/annotation/assemblies.tsv'


TISSUE_MAP = {('hepa', 'liver'): 1,
              ('GM12878', 'CH12'): 1,
              ('CH12', 'GM12878'): 1,
              ('K562', 'MEL'): 1,
              ('MEL', 'K562'): 1,
              ('liver', 'hepa'): 1,
              ('ESE14', 'H1hESC'): 1,
              ('H1hESC', 'ESE14'): 1,
              ('blood', 'ncd4'): 1,
              ('ncd4', 'blood'): 1}

CONS_MAP = {'hg19': os.path.join(__CONSERVATION__, 'hg19_pc-genes_phylop_tsswin.h5'),
            'mm9': os.path.join(__CONSERVATION__, 'mm9_pc-genes_phylop_tsswin.h5')}


def sample_type(sample):
    """
    :param sample:
    :return:
    """
    imm_lut = {'MEL', 'K562', 'GM12878', 'CH12'}
    stem_lut = {'ESE14', 'H1hESC'}
    if sample in imm_lut:
        return 'immortal'
    elif sample in stem_lut:
        return 'stem'
    else:
        return 'primary'


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()

    parser.add_argument('--datasets', '-ds', type=str, default=__DATASET_FILE__, dest='datasets')
    parser.add_argument('--match-types', '-mt', type=str, default=__MATCHTYPES__, dest='matchtypes')
    parser.add_argument('--base-dir', '-bd', type=str, default=__APPLY_ROOT__, dest='applyroot')

    parser.add_argument('--orthologs', '-orth', type=str, default=__ORTHOLOGS__, dest='orthologs')
    parser.add_argument('--assemblies', '-assm', type=str, default=__TAB_ASSEMBLY__, dest='assemblies')
    parser.add_argument('--cons-dir', '-cons', type=str, default=__CONSERVATION__, dest='consdir')
    parser.add_argument('--select-cons', '-slc', type=str, dest='selectcons')

    parser.add_argument('--expression', '-exp', type=str, default=__AGG_EXPRESSION__, dest='expression')

    parser.add_argument('--output', '-o', type=str, dest='outputfile')

    args = parser.parse_args()
    return args


def annotate_test_output(base_path, ds_ids, match_types, req_models):
    """
    :param base_path:
    :param ds_ids:
    :param match_types:
    :return:
    """
    collector = col.defaultdict(list)
    match_count = col.defaultdict(col.Counter)
    found_models = {'hg19': set(), 'mm9': set()}
    for root, dirs, files in os.walk(base_path):
        if files:
            test_out = fnm.filter(files, '*.json')
            for f in test_out:
                fp = os.path.join(root, f)
                with open(fp, 'r') as infile:
                    dump = js.load(infile)
                testfile = dump['run_info']['data_file']
                _, test_epi, test_trans, _, _ = testfile.split('_')
                test_epi_sample = ds_ids[test_epi]['sample']
                test_trans_sample = ds_ids[test_trans]['sample']
                test_key = test_epi_sample + '-' + test_trans_sample
                try:
                    test_type = match_types[test_key]
                except KeyError:
                    # sys.stderr.write('\nAssuming negative test: {}\n'.format(test_key))
                    test_type = 'neg'
                match_count[(test_trans, test_trans_sample)][test_type] += 1
                modelfile = dump['run_info']['model_file']
                trainfile, _, qry, modelinfo = modelfile.split('.', 3)
                gid, epigenome, train_trans, trg, cell = trainfile.split('_')
                assert trg in ['hg19', 'mm9'], 'Unexpected target species: {}'.format(trg)
                assert ds_ids[epigenome]['sample'] == ds_ids[train_trans]['sample'] == cell, \
                    'Sample mismatch for training data: {}'.format(trainfile)
                assert epigenome == test_epi, 'Epigenome mismatch: {}'.format(f)
                found_models[trg].add((epigenome, train_trans))
                modeltype = modelinfo.split('.')[1]
                infos = {'group': gid, 'train_epi': epigenome, 'train_trans': train_trans,
                         'train_cell': cell, 'target_assm': trg, 'query_assm': qry, 'model_type': modeltype,
                         'test_dataset': testfile, 'train_dataset': trainfile,
                         'metadata': fp, 'test_type': test_type, 'test_epi': test_epi,
                         'test_trans': test_trans, 'test_epi_sample': test_epi_sample,
                         'test_trans_sample': test_trans_sample, 'featcomp': gid[-2:]}
                collector[(trg, qry, test_epi, test_trans)].append(infos)
    for assm, models in req_models.items():
        models_trained = found_models[assm]
        for m in models:
            if m not in models_trained:
                sys.stderr.write('\nRequired model missing: {} - {} - {}\n'.format(assm, ds_ids[m[0]]['sample'], m))
    return collector


def identify_missing_tests(assm_trans, models, tests):
    """
    :param assm_trans:
    :param models:
    :param tests:
    :return:
    """
    for trg, train_pairs in models.items():

        for qry, trans in assm_trans.items():
            if trg == qry:
                continue
            if trg == 'mm9' and qry == 'canFam3':
                continue
            qry_trans = [t[0] for t in trans]
            qry_types = dict((t[0], sample_type(t[1])) for t in trans)
            train_epi = [t[0] for t in train_pairs]
            done = set()
            for te in train_epi:
                runs = [k for k in tests.keys() if k[0] == trg and k[1] == qry and k[2] == te]
                [done.add(r[3]) for r in runs]
                delta = set(qry_trans) - done
                delta = set([(d, qry_types[d]) for d in delta if qry_types[d] in ['primary', 'stem']])
                print('======')
                print(trg, ' --> ', qry, ': ', te)
                print('Missing: ', delta)
    return 0


def extract_switching_genes(basepath):
    """
    :param basepath:
    :return:
    """
    expfiles = os.listdir(basepath)
    switch_genes = dict()
    for ef in expfiles:
        assm = ef.split('_')[0]
        fpath = os.path.join(basepath, ef)
        to_drop = []
        with pd.HDFStore(fpath, 'r') as hdf:
            data = hdf['/norm/tpm']
            for c in data.columns:
                cell = c.split('_')[2]
                st = sample_type(cell)
                if st == 'immortal':
                    to_drop.append(c)
            if len(to_drop) > 0:
                data.drop(to_drop, axis=1, inplace=True)
            always_on = (data >= 1).all(axis=1)
            always_off = (data < 1).all(axis=1)
            subset = data.loc[~ np.logical_or(always_on, always_off), :]
            gene_names = subset.index.tolist()
            switch_genes[assm] = gene_names
    return switch_genes


def read_dataset_ids(fpath):
    """
    :param fpath:
    :return:
    """
    dsids = dict()
    train_pairs = {'hg19': col.defaultdict(set),
                   'mm9': col.defaultdict(set)}
    assm_trans = col.defaultdict(list)
    with open(fpath, 'r') as infile:
        rows = csv.DictReader(infile, delimiter='\t')
        for r in rows:
            dsids[r['id']] = {'assembly': r['assembly'], 'sample': r['biosample']}
            if r['assembly'] in ['hg19', 'mm9']:
                if sample_type(r['biosample']) in ['primary', 'stem']:
                    train_pairs[r['assembly']][r['biosample']].add(r['id'])
            if r['id'].startswith('T'):
                assm_trans[r['assembly']].append((r['id'], r['biosample']))
    req_models = {'hg19': set(), 'mm9': set()}
    for assm, samples in train_pairs.items():
        for sample, datasets in samples.items():
            epigenomes = [d for d in datasets if d.startswith('E')]
            transcriptomes = [d for d in datasets if d.startswith('T')]
            for e in epigenomes:
                for t in transcriptomes:
                    req_models[assm].add((e, t))
    return dsids, req_models, assm_trans


def load_gene_orthologs(fpath, spec_a, spec_b):
    """
    :param fpath:
    :param spec_a:
    :param spec_b:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        grp = hdf['/auto/groups']
        grp = grp.loc[:, ['{}_name'.format(spec_a), '{}_name'.format(spec_b)]]
        pairs = hdf['/auto/pairs/{}/{}'.format(spec_a, spec_b)]
        pairs = pairs.loc[:, ['{}_name'.format(spec_a), '{}_name'.format(spec_b)]]
    return pairs, grp


def load_conservation_scores(fpath, select):
    """
    :param fpath:
    :param select:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        data = hdf['genes']
        data.index = data['name']
        cons_cols = [c for c in data.columns if c.endswith(select)]
        data = data.loc[:, cons_cols]
    new_names = []
    for c in data.columns:
        if c.startswith('grp'):
            new_names.append('cons_level')
        elif c.startswith('rnk'):
            new_names.append('cons_rank')
        elif c.startswith('ext'):
            new_names.append('cons_score')
        else:
            raise ValueError('Unexpected column in cons dataframe: {}'.format(c))
    data.columns = new_names
    return data


def collect_perf_metrics(data, metadata, row_name):
    """
    :param data:
    :param metadata:
    :return:
    """
    # save as string to avoid serialization problems due to
    # mixed-integer dtype
    # metadata[prefix + 'num_active'] = str((data['true_class'] == 1).sum())
    # metadata[prefix + 'num_inactive'] = str((data['true_class'] == 0).sum())
    # metadata[prefix + 'tp'] = str(data['tp'].sum())
    # metadata[prefix + 'tn'] = str(data['tn'].sum())
    # metadata[prefix + 'true'] = str(int(metadata[prefix + 'tp']) + int(metadata[prefix + 'tn']))
    # metadata[prefix + 'fp'] = str(data['fp'].sum())
    # metadata[prefix + 'fn'] = str(data['fn'].sum())
    # metadata[prefix + 'false'] = str(int(metadata[prefix + 'fp']) + int(metadata[prefix + 'fn']))
    # metadata[prefix + 'p060'] = str(data['true_prob_060'].sum())
    # metadata[prefix + 'p070'] = str(data['true_prob_070'].sum())
    # metadata[prefix + 'p080'] = str(data['true_prob_080'].sum())
    # metadata[prefix + 'p090'] = str(data['true_prob_090'].sum())
    # metadata[prefix + 'p100'] = str(data['true_prob_100'].sum())
    col_sums = data.sum(axis=0, numeric_only=True)
    labels, values = [], []
    for col, val in col_sums.iteritems():
        if col == 'true_class':
            pos_class = val
            neg_class = data.shape[0] - pos_class
            labels.extend(['pos_class', 'neg_class'])
            values.extend([pos_class, neg_class])
        elif col.startswith('true_prob') or col.startswith('pred_prob'):
            labels.append(col)
            values.append(val)
        elif col == 'pred_class':
            continue
        elif '_name' in col:
            continue
        elif col.startswith('cons'):
            continue
        elif col.startswith('true_class_prob') or col.startswith('pred_class_prob') or col.startswith('pos_class_prob'):
            continue
        else:
            labels.append(col)
            values.append(val)
    positives = data['tp'].sum() + data['tn'].sum()
    negatives = data['fp'].sum() + data['fn'].sum()
    labels.extend(['positives', 'negatives'])
    values.extend([positives, negatives])
    prec, recall, f1score, _ = precision_recall_fscore_support(data['true_class'], data['pred_class'],
                                                               beta=1.0, pos_label=1, average='macro')
    acc = accuracy_score(data['true_class'], data['pred_class'])
    sens_tpr = data['tp'].sum() / (data['tp'].sum() + data['fn'].sum())
    spec_tnr = data['tn'].sum() / (data['tn'].sum() + data['fp'].sum())
    ppv = data['tp'].sum() / (data['tp'].sum() + data['fp'].sum())
    npv = data['tn'].sum() / (data['tn'].sum() + data['fn'].sum())
    fpr = data['fp'].sum() / (data['fp'].sum() + data['tn'].sum())
    fnr = data['fn'].sum() / (data['tp'].sum() + data['fn'].sum())
    fdr = data['fp'].sum() / (data['tp'].sum() + data['fp'].sum())
    labels.extend(['precision', 'recall', 'f1score', 'accuracy', 'sensitivity',
                   'true_pos_rate', 'specificity', 'true_neg_rate', 'pos_pred_value',
                   'neg_pred_value', 'false_pos_rate', 'false_neg_rate',
                   'false_discovery_rate', 'auc_roc'])
    auc = roc_auc_score(data['true_class'], data['pos_class_prob'], average='macro')
    values.extend([prec, recall, f1score, acc, sens_tpr, sens_tpr,
                   spec_tnr, spec_tnr, ppv, npv, fpr, fnr, fdr, auc])
    if metadata is None:
        metadata = pd.DataFrame([values], index=[row_name], columns=labels, dtype=np.float32)
    else:
        new_row = pd.DataFrame([values], index=[row_name], columns=labels, dtype=np.float32)
        metadata = pd.concat([metadata, new_row], ignore_index=False, axis=0)
    return metadata


def record_cons_performance(dataset):
    """
    :param dataset:
    :return:
    """
    cons_levels = sorted(dataset['cons_level'].unique())
    metrics = ['selected', 'relevant', 'pos_class', 'neg_class',
               'positives', 'negatives', 'precision', 'recall', 'f1score',
               'accuracy', 'sensitivity', 'true_pos_rate',
               'specificity', 'true_neg_rate', 'pos_pred_value',
               'neg_pred_value', 'false_pos_rate', 'false_neg_rate',
               'false_discovery_rate', 'auc_roc']
    cons_scoring = pd.DataFrame(np.zeros((len(metrics), len(cons_levels)), dtype=np.float32),
                                index=metrics, columns=cons_levels)
    for lvl in cons_levels:
        lvl_sub = dataset.loc[dataset['cons_level'] >= lvl, :]
        lvl_metrics = [lvl_sub.shape[0], lvl_sub.shape[0]]
        lvl_posclass = lvl_sub['true_class'].sum()
        lvl_negclass = lvl_sub.shape[0] - lvl_posclass
        lvl_pos = (lvl_sub['pred_class'] == lvl_sub['true_class']).sum()
        lvl_neg = lvl_sub.shape[0] - lvl_pos
        lvl_prec, lvl_recall, lvl_f1, _ = precision_recall_fscore_support(lvl_sub['true_class'],
                                                                          lvl_sub['pred_class'],
                                                                          beta=1.0, pos_label=1,
                                                                          average='macro')
        lvl_acc = accuracy_score(dataset['true_class'], dataset['pred_class'])
        lvl_sens_tpr = dataset['tp'].sum() / (dataset['tp'].sum() + dataset['fn'].sum())
        lvl_spec_tnr = dataset['tn'].sum() / (dataset['tn'].sum() + dataset['fp'].sum())
        lvl_ppv = dataset['tp'].sum() / (dataset['tp'].sum() + dataset['fp'].sum())
        lvl_npv = dataset['tn'].sum() / (dataset['tn'].sum() + dataset['fn'].sum())
        lvl_fpr = dataset['fp'].sum() / (dataset['fp'].sum() + dataset['tn'].sum())
        lvl_fnr = dataset['fn'].sum() / (dataset['tp'].sum() + dataset['fn'].sum())
        lvl_fdr = dataset['fp'].sum() / (dataset['tp'].sum() + dataset['fp'].sum())
        lvl_auc = roc_auc_score(dataset['true_class'], dataset['pos_class_prob'], average='macro')
        lvl_metrics.extend([lvl_posclass, lvl_negclass, lvl_pos, lvl_neg,
                            lvl_prec, lvl_recall, lvl_f1,
                            lvl_acc, lvl_sens_tpr, lvl_sens_tpr, lvl_spec_tnr, lvl_spec_tnr,
                            lvl_ppv, lvl_npv, lvl_fpr, lvl_fnr, lvl_fdr, lvl_auc])
        cons_scoring[lvl] = lvl_metrics
    return cons_scoring


def compute_pos_class_prob(dataset):
    """
    :param dataset:
    :return:
    """
    true_is_one = np.array(dataset['true_class'] == 1, dtype=np.bool)
    true_is_zero = np.array(dataset['true_class'] == 0, dtype=np.bool)
    pos_class_prob = np.zeros(true_is_one.size, dtype=np.float32)
    pos_class_prob[true_is_one] = dataset['true_class_prob'][true_is_one]
    pos_class_prob[true_is_zero] = 1 - dataset['true_class_prob'][true_is_zero]
    assert pos_class_prob.max() <= 1, 'Max pos. class prob. false: {}'.format(pos_class_prob.max())
    assert pos_class_prob.min() >= 0, 'Min pos. class prob. false: {}'.format(pos_class_prob.min())
    return pos_class_prob


def process_run_metadata(collected_infos, genes_switching, assm_map, gene_orthologs, cons_select, outputfile):
    """
    :param collected_infos:
    :param genes_switching:
    :param gene_orthologs:
    :param cons_select:
    :param outputfile:
    :return:
    """
    current_pairing = None, None
    pair_ortho, group_ortho = None, None
    filemode = 'w'
    for (trg, qry, epi, exp) in sorted(collected_infos.keys()):
        if (trg, qry) != current_pairing:
            trg_spec, qry_spec = assm_map[trg], assm_map[qry]
            pair_ortho, group_ortho = load_gene_orthologs(gene_orthologs, trg_spec, qry_spec)
            if (trg, qry) == ('hg19', 'mm9'):
                conservation = load_conservation_scores(CONS_MAP['mm9'], cons_select)
            elif (trg, qry) == ('mm9', 'hg19'):
                conservation = load_conservation_scores(CONS_MAP['hg19'], cons_select)
            else:
                conservation = None
        for run in collected_infos[(trg, qry, epi, exp)]:
            run['target_spec'] = trg_spec
            run['query_spec'] = qry_spec
            if run['model_type'] == 'seq':
                outpath = '/'.join(['', '{test_type}', '{model_type}', '{target_spec}', '{query_spec}',
                                    '{train_trans}_{target_assm}_{train_cell}',
                                    '{test_trans}_{query_assm}_{test_trans_sample}'])
            else:
                outpath = '/'.join(['', '{test_type}', '{model_type}', '{target_spec}', '{query_spec}',
                                    '{train_epi}_{train_trans}_{target_assm}_{train_cell}',
                                    '{test_epi}_{test_trans}_{query_assm}_{test_trans_sample}'])
            outpath = outpath.format(**run)
            gene_col = '{}_name'.format(qry_spec)
            with open(run['metadata'], 'r') as run_info:
                dump = js.load(run_info)
            run['metadata'] = os.path.basename(run['metadata'])
            gene_names = dump['sample_info']['names']
            true_class = np.array(dump['testing_info']['targets']['true'], dtype=np.int8)
            pred_class = np.array(dump['testing_info']['targets']['pred'], dtype=np.int8)
            # predicted probabilities for true class labels
            # if p(label) > 0.5, the correct class has been predicted
            true_class_prob = np.array(dump['testing_info']['probabilities']['true'], dtype=np.float32)
            pred_class_prob = np.array(dump['testing_info']['probabilities']['pred'], dtype=np.float32)
            prob_steps = list(range(50, 100, 5))
            prob_steps.extend([99, 100])
            prob_labels = list(map(str, prob_steps))
            prob_steps = np.array(prob_steps, dtype=np.float32)
            prob_steps /= 100.
            df = pd.DataFrame(None, index=gene_names, dtype=np.int8)
            df.index.name = gene_col
            df['true_class'] = true_class
            df['pred_class'] = pred_class
            df['true_class_prob'] = true_class_prob
            df['pred_class_prob'] = pred_class_prob
            df['pos_class_prob'] = compute_pos_class_prob(df)
            df['tp'] = np.array(np.logical_and(df['true_class'] == 1, df['pred_class'] == 1), dtype=np.int8)
            df['tn'] = np.array(np.logical_and(df['true_class'] == 0, df['pred_class'] == 0), dtype=np.int8)
            df['fp'] = np.array(np.logical_and(df['true_class'] == 0, df['pred_class'] == 1), dtype=np.int8)
            df['fn'] = np.array(np.logical_and(df['true_class'] == 1, df['pred_class'] == 0), dtype=np.int8)
            if conservation is not None:
                assert df.shape[0] == conservation.shape[0], \
                    'Differing number of genes: {} vs {}'.format(df.shape[0], conservation.shape[0])
                df = df.join(conservation, on=None, how='left')
                cons_scoring = record_cons_performance(df)
            else:
                cons_scoring = None
            # add indicator columns for membership of genes to certain probability intervals
            # for true predictions
            for idx, (s, l) in enumerate(zip(prob_steps[1:-1], prob_labels[1:-1]), start=1):
                indicator = np.array(np.logical_and(prob_steps[idx-1] < true_class_prob,
                                                    true_class_prob <= s), dtype=np.int8)
                df['true_prob_iv_{}-{}'.format(prob_labels[idx-1], l)] = indicator
                indicator = np.array(prob_steps[idx-1] < true_class_prob, dtype=np.int8)
                df['true_prob_lo_{}'.format(prob_labels[idx-1])] = indicator

            # for predicted class - also captures false predictions
            for idx, (s, l) in enumerate(zip(prob_steps[1:-1], prob_labels[1:-1]), start=1):
                indicator = np.array(np.logical_and(prob_steps[idx-1] < pred_class_prob,
                                                    pred_class_prob <= s), dtype=np.int8)
                df['pred_prob_iv_{}-{}'.format(prob_labels[idx-1], l)] = indicator
                indicator = np.array(prob_steps[idx-1] < pred_class_prob, dtype=np.int8)
                df['pred_prob_lo_{}'.format(prob_labels[idx-1])] = indicator

            df['group_ortho'] = np.array(df.index.isin(group_ortho[gene_col]), dtype=np.int8)
            df['pair_ortho'] = np.array((df.index.isin(pair_ortho[gene_col])), dtype=np.int8)
            df['switching'] = np.array(df.index.isin(genes_switching[qry]), dtype=np.int8)
            perf = collect_perf_metrics(df, None, 'perf_wg')
            perf = collect_perf_metrics(df.loc[df['switching'] == 1, :].copy(), perf, 'perf_switch')
            perf = collect_perf_metrics(df.loc[df['group_ortho'] == 1, :].copy(), perf, 'perf_group')
            perf = collect_perf_metrics(df.loc[df['pair_ortho'] == 1, :].copy(), perf, 'perf_pair')
            mdf = pd.DataFrame.from_dict(run, orient='index', dtype='object')
            mdf.index.name = 'key'
            mdf.columns = ['value']
            stored_keys = set()
            with pd.HDFStore(outputfile, filemode, complib='blosc', complevel=9) as hdf:
                data_path = outpath + '/data'
                assert data_path not in stored_keys, 'Data path duplicate: {}'.format(data_path)
                hdf.put(data_path, df, format='fixed')
                stored_keys.add(data_path)
                md_path = outpath + '/metadata'
                assert md_path not in stored_keys, 'Metadata path duplicate: {}'.format(md_path)
                hdf.put(md_path, mdf, format='table')
                stored_keys.add(md_path)
                perf_path = outpath + '/perf'
                assert perf_path not in stored_keys, 'Performance path duplicate: {}'.format(perf_path)
                hdf.put(perf_path, perf)
                # record values for ROC plotting
                fpr, tpr, thresholds = roc_curve(df['true_class'], df['pos_class_prob'], pos_label=1)
                df_roc = pd.DataFrame([fpr, tpr, thresholds], index=['fpr', 'tpr', 'thresholds'], dtype=np.float32)
                roc_path = outpath + '/roc'
                assert roc_path not in stored_keys, 'ROC curve path duplicate: {}'.format(roc_path)
                hdf.put(roc_path, df_roc, format='fixed')
                if cons_scoring is not None:
                    cons_path = outpath + '/cons'
                    hdf.put(cons_path, cons_scoring, format='fixed')
                hdf.flush()
            filemode = filemode.replace('w', 'a')
    return 0


def main():
    """
    :return:
    """
    args = parse_command_line()
    switch = extract_switching_genes(args.expression)
    dsids, req_models, assm_trans = read_dataset_ids(args.datasets)
    mtypes = js.load(open(args.matchtypes, 'r'))['matchtypes']
    collect_test = annotate_test_output(args.applyroot, dsids, mtypes, req_models)
    lut_assm = dict()
    with open(args.assemblies, 'r') as tab:
        rows = csv.DictReader(tab, delimiter='\t')
        for row in rows:
            lut_assm[row['assembly']] = row['common_name']
    # no more missing tests, skip that step
    # _ = identify_missing_tests(assm_trans, req_models, collect_test)
    _ = process_run_metadata(collect_test, switch, lut_assm, args.orthologs, args.selectcons, args.outputfile)
    return 0


if __name__ == '__main__':
    try:
        _ = main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
