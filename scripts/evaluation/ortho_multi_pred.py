#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp

import pandas as pd
import numpy as np
from scipy.stats import kendalltau as ktau
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import r2_score as r2s


SPECIES_MAP = {'human': 'hg19',
               'mouse': 'mm9',
               'dog': 'canFam3',
               'cow': 'bosTau7',
               'chicken': 'galGal3',
               'pig': 'susScr2'}

TISSUE_MAP = {('hepa', 'liver'): 1,
              ('GM12878', 'CH12'): 1,
              ('K562', 'MEL'): 1,
              ('liver', 'hepa'): 1,
              ('ESE14', 'H1hESC'): 1,
              ('H1hESC', 'ESE14'): 1}


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()

    parser.add_argument('--exp-files', '-exp', type=str, nargs='+', dest='expfiles')
    parser.add_argument('--ortho-file', '-orth', type=str, dest='orthofile')
    parser.add_argument('--cons-file', '-cons', type=str, nargs='+', dest='consfiles')
    parser.add_argument('--select-cons', '-slc', type=str, default='', dest='selectcons')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--group-root', '-gr', type=str, dest='grouproot', default='/auto/pairs')
    parser.add_argument('--tpm-threshold', '-tpm', type=int, dest='threshold', default=1)

    args = parser.parse_args()
    return args


def kendall_tau_scorer(x, y):
    """
    :param x:
    :param y:
    :return:
    """
    statistic, p_value = ktau(x, y, nan_policy='raise')
    return statistic


def normalize_col_names(cnames, species):
    """
    :param cnames:
    :return:
    """
    # TS18_susScr2_kidney_mRNA
    norm = []
    datacols = []
    for c in cnames:
        if c == 'gene_name':
            norm.append(c)
        else:
            comp = c.split('_')
            new_c = '_'.join([comp[0], species, comp[2]])
            norm.append(new_c)
            datacols.append(new_c)
    return norm, datacols


def load_conservation_scores(consfiles, species, select):
    """
    :param consfiles:
    :param species:
    :return:
    """
    assm = SPECIES_MAP[species]
    to_load = [f for f in consfiles if os.path.basename(f).startswith(assm)]
    assert len(to_load) == 1, 'Cannot load conservation data for assembly {}: {}'.format(assm, consfiles)
    with pd.HDFStore(to_load[0], 'r') as hdf:
        data = hdf['genes']
        cons_cols = [c for c in data.columns if c.endswith(select)]
        data = data[['name'] + cons_cols]
    new_names = []
    for c in data.columns:
        if c.startswith('grp'):
            new_names.append('cons_level')
        elif c.startswith('rnk'):
            new_names.append('cons_rank')
        elif c.startswith('ext'):
            new_names.append('cons_score')
        elif c == 'name':
            new_names.append('{}_name'.format(species))
        else:
            raise ValueError('Unexpected column in cons dataframe: {}'.format(c))
    data.columns = new_names
    assert not data.empty, 'Conservation score DF is empty'
    assert data.notnull().all(axis=1).all(), 'Loaded conservation scores contain NaN'
    return data


def identify_species_pairs(fpath, grouproot):
    """
    :param fpath:
    :return:
    """
    if not grouproot.startswith('/'):
        grouproot = '/' + grouproot
    with pd.HDFStore(fpath, 'r') as hdf:
        load_groups = [k for k in hdf.keys() if k.startswith(grouproot)]
        assert load_groups, 'No data groups to load from HDF for root: {}'.format(grouproot)
    species = []
    for path in load_groups:
        tmp = path.replace(grouproot, '')
        a, b = os.path.split(tmp)
        assert a != b, 'Invalid species pair: {} and {}'.format(a, b)
        species.append((a.strip('/'), b.strip('/'), path))
    return species


def load_expression_data(species, expfiles):
    """
    :param species:
    :param expfiles:
    :return:
    """
    assm = SPECIES_MAP[species.strip('/')]
    expf = [fp for fp in expfiles if os.path.basename(fp).startswith(assm)]
    assert len(expf) == 1, 'Could not identify expression file: {}'.format(expf)
    with pd.HDFStore(expf[0], 'r') as hdf:
        tpm = hdf['/norm/tpm']
        ranks = hdf['/norm/ranks']
    assert sorted(tpm.columns) == sorted(ranks.columns), 'Incompatible column names in expression data'
    return tpm, ranks


def tissue_match(a, b):
    """
    :param a:
    :param b:
    :return:
    """
    m = (a == b) or ((a, b) in TISSUE_MAP) or ((b, a) in TISSUE_MAP)
    return m


def make_ortholog_pred(species_a, tpm_a, ranks_a,
                       species_b, tpm_b, ranks_b,
                       orthologs, threshold, tsscons,
                       comp, outputfile):
    """

    Assumed directionality: from species A to species B - i.e. labels of
    species B are "true", labels of species A are "pred"

    :param tpm_a:
    :param ranks_a:
    :param tpm_b:
    :param ranks_b:
    :param orthologs:
    :param outputfile:
    :return:
    """
    assert pd.notnull(tpm_a).all(axis=0).all(), 'Data for A {} contain NaN/NULL values'.format(species_a)
    assert pd.notnull(tpm_b).all(axis=0).all(), 'Data for B {} contain NaN/NULL values'.format(species_b)
    name_a, name_b = '{}_name'.format(species_a), '{}_name'.format(species_b)
    # for HCOP, this may happen
    if name_a not in orthologs.columns or name_b not in orthologs.columns:
        return
    tpm_threshold = threshold
    ortho_rel = orthologs.loc[:, [name_a, name_b]].copy()
    select_a = tpm_a.index.isin(orthologs.loc[:, name_a])
    select_b = tpm_b.index.isin(orthologs.loc[:, name_b])

    sub_tpm_a = tpm_a.loc[select_a, :].copy()
    sub_tpm_a[name_a] = sub_tpm_a.index

    sub_tpm_b = tpm_b.loc[select_b, :].copy()
    sub_tpm_b[name_b] = sub_tpm_b.index

    tpm = ortho_rel.merge(sub_tpm_a, how='outer', on=name_a, copy=True)
    tpm = tpm.merge(sub_tpm_b, how='outer', on=name_b, copy=True)
    assert pd.notnull(tpm).all(axis=0).all(),\
        'Merged TPM dataframe contains NULL/NaN: {} and {}'.format(species_a, species_b)

    sub_ranks_a = ranks_a.loc[select_a, :].copy()
    sub_ranks_a[name_a] = sub_ranks_a.index
    sub_ranks_b = ranks_b.loc[select_b, :].copy()
    sub_ranks_b[name_b] = sub_ranks_b.index
    ranks = ortho_rel.merge(sub_ranks_a, how='outer', on=name_a, copy=True)
    ranks = ranks.merge(sub_ranks_b, how='outer', on=name_b, copy=True)
    new_cols = []
    for c in ranks.columns:
        if c.endswith('name'):
            new_cols.append(c)
        else:
            c = c.replace('mRNA', 'rank')
            new_cols.append(c)
    ranks.columns = new_cols

    num_orthologs = tpm.shape[0]
    for ca in sub_tpm_a.columns:
        if ca.endswith('_name'):
            continue
        norm_a = ca.rsplit('_', 1)[0]
        rca = norm_a + '_rank'
        tissue_a = norm_a.rsplit('_', 1)[-1]
        spec_a_labels = tpm[ca] >= tpm_threshold
        spec_a_active = spec_a_labels.sum()
        spec_a_inactive = num_orthologs - spec_a_active
        for cb in sub_tpm_b.columns:
            if cb.endswith('_name'):
                continue
            norm_b = cb.rsplit('_', 1)[0]
            rcb = norm_b + '_rank'
            tissue_b = norm_b.rsplit('_', 1)[-1]
            spec_b_labels = tpm[cb] >= tpm_threshold
            spec_b_active = spec_b_labels.sum()
            spec_b_inactive = num_orthologs - spec_b_active

            # zero_wt = 0.5 / spec_b_inactive
            # one_wt = 0.5 / spec_b_active
            # wt_vec = np.array([zero_wt if val < tpm_threshold else one_wt for val in tpm[cb]], dtype=np.float64)
            # assert np.isclose(wt_vec.sum(), 1), 'Weight vector does not sum to 1: {}...'.format(wt_vec[:5])
            # # record performance metrics for whole dataset
            # acc_score = acc(spec_b_labels, spec_a_labels, sample_weight=wt_vec)

            # average='macro' - both classes have the same weight, irrespective of support
            prec, recall, f1, _ = precision_recall_fscore_support(spec_b_labels, spec_a_labels,
                                                                  beta=1.0, pos_label=1, average='macro')
            r2_score_all = r2s(tpm[cb], tpm[ca])
            ktau_score_all = kendall_tau_scorer(tpm[cb], tpm[ca])

            conservation = None
            if tsscons is not None:
                conservation = tsscons.loc[tsscons[name_b].isin(tpm[name_b]), :].copy()
                cons_levels = sorted(conservation['cons_level'].unique())
                metrics = ['selected', 'relevant', 'pos_class', 'neg_class',
                           'positives', 'negatives', 'precision', 'recall', 'f1score']
                cons_scoring = pd.DataFrame(np.zeros((len(metrics), len(cons_levels)), dtype=np.float32),
                                            index=metrics, columns=cons_levels)
                for lvl in cons_levels:
                    lvl_genes = conservation.loc[conservation['cons_level'] >= lvl, name_b]
                    lvl_metrics = [lvl_genes.shape[0], sub_tpm_b[name_b].isin(lvl_genes).sum()]
                    lvl_sub = tpm.loc[tpm[name_b].isin(lvl_genes), :]
                    lvl_labels_b = lvl_sub[cb] >= tpm_threshold
                    lvl_labels_a = lvl_sub[ca] >= tpm_threshold
                    lvl_posclass = lvl_labels_b.sum()
                    lvl_negclass = lvl_sub.shape[0] - lvl_posclass
                    lvl_pos = (lvl_labels_a == lvl_labels_b).sum()
                    lvl_neg = lvl_sub.shape[0] - lvl_pos
                    lvl_prec, lvl_recall, lvl_f1, _ = precision_recall_fscore_support(lvl_labels_b, lvl_labels_a,
                                                                                      beta=1.0, pos_label=1,
                                                                                      average='macro')
                    lvl_metrics.extend([lvl_posclass, lvl_negclass, lvl_pos, lvl_neg,
                                        lvl_prec, lvl_recall, lvl_f1])
                    cons_scoring[lvl] = lvl_metrics

            ranked = ranks.loc[:, [rca, rcb]].copy().rank(axis=0, method='dense', pct=True)
            deltas = np.abs(ranked[rca] - ranked[rcb]) < 0.05
            num_consistent_all_5 = deltas.sum()
            num_inconsistent_all_5 = deltas.size - num_consistent_all_5

            deltas = np.abs(ranked[rca] - ranked[rcb]) < 0.1
            num_consistent_all_10 = deltas.sum()
            num_inconsistent_all_10 = deltas.size - num_consistent_all_10

            # record performance metrics for subset predicted as active
            active_subset = tpm.loc[spec_a_labels, :].copy()
            assert active_subset.shape[0] == spec_a_active, \
                'Selecting active subset failed: should be {}, is {}'.format(spec_a_active, active_subset.shape[0])
            r2_score_act = r2s(active_subset[cb], active_subset[ca])
            ktau_score_act = kendall_tau_scorer(active_subset[cb], active_subset[ca])

            active_rank_subset = ranks.loc[spec_a_labels, [rca, rcb]].copy().rank(axis=0, method='dense', pct=True)
            deltas = np.abs(active_rank_subset[rca] - active_rank_subset[rcb]) < 0.05
            num_consistent_act_5 = deltas.sum()
            num_inconsistent_act_5 = deltas.size - num_consistent_act_5

            deltas = np.abs(active_rank_subset[rca] - active_rank_subset[rcb]) < 0.1
            num_consistent_act_10 = deltas.sum()
            num_inconsistent_act_10 = deltas.size - num_consistent_act_10

            idx_tp = np.array(np.logical_and(tpm[ca] >= tpm_threshold, tpm[cb] >= tpm_threshold), dtype=np.int8)
            assert np.isfinite(idx_tp).all(), 'True positive label not finite'
            idx_fp = np.array(np.logical_and(tpm[ca] >= tpm_threshold, tpm[cb] < tpm_threshold), dtype=np.int8)
            assert np.isfinite(idx_fp).all(), 'False positive label not finite'
            idx_tn = np.array(np.logical_and(tpm[ca] < tpm_threshold, tpm[cb] < tpm_threshold), dtype=np.int8)
            assert np.isfinite(idx_tn).all(), 'True negative label not finite'
            idx_fn = np.array(np.logical_and(tpm[ca] < tpm_threshold, tpm[cb] >= tpm_threshold), dtype=np.int8)
            assert np.isfinite(idx_fn).all(), 'False negative label not finite'

            num_tp = idx_tp.sum()
            num_fp = idx_fp.sum()
            num_tn = idx_tn.sum()
            num_fn = idx_fn.sum()

            dataset = tpm.loc[:, (name_a, ca, name_b, cb)].copy()
            dataset.columns = [name_a, norm_a, name_b, norm_b]
            dataset = dataset.merge(ranks.loc[:, (name_a, norm_a + '_rank')], on=name_a, how='outer', copy=True)
            dataset = dataset.merge(ranks.loc[:, (name_b, norm_b + '_rank')], on=name_b, how='outer', copy=True)
            assert dataset.notnull().all(axis=1).all(), 'Dataset contains NaN after merging ranks'
            dataset['true_pos'] = idx_tp
            dataset['false_pos'] = idx_fp
            dataset['true_neg'] = idx_tn
            dataset['false_neg'] = idx_fn
            if tsscons is not None:
                dataset = dataset.merge(conservation, on=name_b, how='outer', copy=True)
                assert dataset.notnull().all(axis=1).all(), 'Dataset contains NaN after merging conservation scores'

            metadata = {'num_orthologs': num_orthologs,
                        name_a + '_active': spec_a_active, name_a + '_inactive': spec_a_inactive,
                        name_b + '_inactive': spec_b_inactive, name_b + '_inactive': spec_b_inactive,
                        'perf_prec_all': prec, 'perf_recall_all': recall, 'perf_f1score_all': f1,
                        'perf_r2_all': r2_score_all, 'perf_ktau_all': ktau_score_all,
                        'perf_r2_active': r2_score_act, 'perf_ktau_active': ktau_score_act,
                        'perf_num_tp': num_tp, 'perf_num_fp': num_fp, 'perf_num_tn': num_tn, 'perf_num_fn': num_fn,
                        'perf_num_consistent_all_5': num_consistent_all_5, 'perf_num_inconsistent_all_5': num_inconsistent_all_5,
                        'perf_num_consistent_all_10': num_consistent_all_10, 'perf_num_inconsistent_all_10': num_inconsistent_all_10,
                        'perf_num_consistent_act_5': num_consistent_act_5, 'perf_num_inconsistent_act_5': num_inconsistent_act_5,
                        'perf_num_consistent_act_10': num_consistent_act_10, 'perf_num_inconsistent_act_10': num_inconsistent_act_10}
            metadata = pd.DataFrame.from_dict(metadata, orient='index')

            if tissue_match(tissue_a, tissue_b):
                base_group = '/pos/{}/{}/{}/{}/{}'.format(comp, species_a, species_b, norm_a, norm_b)
                with pd.HDFStore(outputfile, 'a', complib='blosc', complevel=9) as hdf:
                    hdf.put(os.path.join(base_group, 'data'), dataset, format='fixed')
                    hdf.put(os.path.join(base_group, 'metadata'), metadata, format='fixed')
                    if tsscons is not None:
                        hdf.put(os.path.join(base_group, 'cons'), cons_scoring, format='fixed')
                    hdf.flush()
            else:
                base_group = '/neg/{}/{}/{}/{}/{}'.format(comp, species_a, species_b, norm_a, norm_b)
                with pd.HDFStore(outputfile, 'a', complib='blosc', complevel=9) as hdf:
                    hdf.put(os.path.join(base_group, 'data'), dataset, format='fixed')
                    hdf.put(os.path.join(base_group, 'metadata'), metadata, format='fixed')
                    if tsscons is not None:
                        hdf.put(os.path.join(base_group, 'cons'), cons_scoring, format='fixed')
                    hdf.flush()
    return


def make_group_matrix(tpm_data, ortho_group):
    """
    :param tpm_data:
    :param ortho_group:
    :return:
    """

    for spec, tpm in tpm_data:
        shared_keys = set(tpm.columns).intersection(ortho_group.columns)
        ortho_group = ortho_group.merge(tpm, how='outer', on=list(shared_keys),
                                        suffixes=('', ''), indicator=False, copy=True)
    # normalize columns
    new_cols = []
    to_drop = []
    for c in ortho_group.columns:
        if c.endswith('_name') or c.endswith('_symbol'):
            new_cols.append(c)
        elif c.endswith('_mRNA'):
            new_c = c.replace('_mRNA', '')
            new_cols.append(new_c)
        else:
            to_drop.append(c)
            new_cols.append(c)
    ortho_group.columns = new_cols
    if to_drop:
        ortho_group.drop(to_drop, axis=1, inplace=True)
    ortho_group.dropna(how='any', axis=0, inplace=True)
    return ortho_group


def main():
    """
    :return:
    """
    args = parse_command_line()

    species_pairs = identify_species_pairs(args.orthofile, args.grouproot)
    all_tpms = []
    spec_done = set()
    for spec_a, spec_b, group in species_pairs:
        tpm_a, ranks_a = load_expression_data(spec_a, args.expfiles)
        tpm_b, ranks_b = load_expression_data(spec_b, args.expfiles)
        if (spec_a, spec_b) in [('human', 'mouse'), ('mouse', 'human')]:
            cons_scores = load_conservation_scores(args.consfiles, spec_b, args.selectcons)
        else:
            cons_scores = None
        with pd.HDFStore(args.orthofile, 'r') as hdf:
            orthologs = hdf[group]
        make_ortholog_pred(spec_a, tpm_a, ranks_a,
                           spec_b, tpm_b, ranks_b,
                           orthologs, args.threshold,
                           cons_scores, 'pair', args.outputfile)
        with pd.HDFStore(args.orthofile, 'r') as hdf:
            orthologs = hdf['/auto/groups']
        make_ortholog_pred(spec_a, tpm_a, ranks_a,
                           spec_b, tpm_b, ranks_b,
                           orthologs, args.threshold,
                           None, 'group', args.outputfile)
        if spec_a not in spec_done:
            tpm_a['{}_name'.format(spec_a)] = tpm_a.index
            tpm_a = tpm_a.reset_index(drop=True, inplace=False)
            all_tpms.append((spec_a, tpm_a))
            spec_done.add(spec_a)
        if spec_b not in spec_done:
            tpm_b['{}_name'.format(spec_b)] = tpm_b.index
            tpm_b = tpm_b.reset_index(drop=True, inplace=False)
            all_tpms.append((spec_b, tpm_b))
            spec_done.add(spec_b)

    with pd.HDFStore(args.orthofile, 'r') as hdf:
        group_genes = hdf['/auto/groups']
    mat = make_group_matrix(all_tpms, group_genes)
    with pd.HDFStore(args.outputfile, 'a', complevel=9, complib='blosc') as hdf:
        hdf.put('/matrix/auto/data', mat, format='table')
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
