#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp

import pandas as pd
import numpy as np
from scipy.stats import kendalltau as ktau
from sklearn.metrics import accuracy_score as acc
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
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--group-root', '-gr', type=str, dest='grouproot', default='/auto/pairs')

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


def process_species_data(species, tpm, rank, orthologs, metadata):
    """
    :param tpm:
    :param rank:
    :param orthologs:
    :return:
    """
    new_cols, data_cols = normalize_col_names(tpm.columns, species)
    select_cols = data_cols + ['og_id']
    tpm.columns = new_cols
    tpm_orth = orthologs.merge(tpm, on='gene_name', copy=True)
    tpm_strict = tpm_orth.loc[tpm_orth['group_size'] == 5, select_cols]
    tpm_strict = tpm_strict.groupby('og_id').mean()
    tpm_balanced = tpm_orth.loc[tpm_orth['group_balanced'] == 1, select_cols + ['species', 'gene_name']]
    for s in tpm_balanced['species'].unique():
        gbal = tpm_balanced.loc[tpm_balanced['species'] == s, 'gene_name'].unique().size
        metadata['balanced_genes_{}'.format(s)] = int(gbal)
    tpm_balanced.drop(['species', 'gene_name'], axis=1, inplace=True)
    tpm_balanced = tpm_balanced.groupby('og_id').mean()
    tpm_all = tpm_orth.loc[:, select_cols]
    tpm_all = tpm_all.groupby('og_id').mean()
    # same for ranks
    new_cols, data_cols = normalize_col_names(rank.columns, species)
    select_cols = data_cols + ['og_id']
    rank.columns = new_cols
    rank_orth = orthologs.merge(rank, on='gene_name', copy=True)
    rank_strict = rank_orth.loc[rank_orth['group_size'] == 5, select_cols]
    rank_strict = rank_strict.groupby('og_id').median()
    rank_strict = rank_strict.round(decimals=0)
    rank_balanced = rank_orth.loc[tpm_orth['group_balanced'] == 1, select_cols]
    rank_balanced = rank_balanced.groupby('og_id').median()
    rank_balanced = rank_balanced.round(decimals=0)
    rank_all = rank_orth.loc[:, select_cols]
    rank_all = rank_all.groupby('og_id').median()
    rank_all = rank_all.round(decimals=0)
    tpm_data = {'strict': tpm_strict, 'balanced': tpm_balanced, 'all': tpm_all}
    rank_data = {'strict': rank_strict, 'balanced': rank_balanced, 'all': rank_all}
    return tpm_data, rank_data, metadata


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
                       orthologs, outputfile):
    """
    :param tpm_a:
    :param ranks_a:
    :param tpm_b:
    :param ranks_b:
    :param orthologs:
    :param outputfile:
    :return:
    """
    name_a, name_b = '{}_name'.format(species_a), '{}_name'.format(species_b)
    ortho_rel = orthologs.loc[:, [name_a, name_b]].copy()
    select_a = tpm_a.index.isin(orthologs.loc[:, name_a])
    select_b = tpm_b.index.isin(orthologs.loc[:, name_b])

    sub_tpm_a = tpm_a.loc[select_a, :].copy()
    sub_tpm_a[name_a] = sub_tpm_a.index

    sub_tpm_b = tpm_b.loc[select_b, :].copy()
    sub_tpm_b[name_b] = sub_tpm_b.index
    tpm = ortho_rel.merge(sub_tpm_a, how='outer', on=name_a, copy=True)
    tpm = tpm.merge(sub_tpm_b, how='outer', on=name_b, copy=True)

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
        tissue_a = norm_a.rsplit('_', 1)[-1]
        spec_a_labels = tpm[ca] >= 1
        spec_a_active = spec_a_labels.sum()
        spec_a_inactive = num_orthologs - spec_a_active
        for cb in sub_tpm_b.columns:
            if cb.endswith('_name'):
                continue
            norm_b = cb.rsplit('_', 1)[0]
            tissue_b = norm_b.rsplit('_', 1)[-1]
            spec_b_labels = tpm[cb] >= 1
            spec_b_active = spec_b_labels.sum()
            spec_b_inactive = num_orthologs - spec_b_active
            zero_wt = 0.5 / spec_b_inactive
            one_wt = 0.5 / spec_b_active
            wt_vec = [zero_wt if val < 1 else one_wt for val in tpm[cb]]
            # record performance metrics for whole dataset
            acc_score = acc(spec_b_labels, spec_a_labels, sample_weight=wt_vec)
            r2_score_all = r2s(tpm[cb], tpm[ca])
            ktau_score_all = kendall_tau_scorer(tpm[cb], tpm[ca])

            # record performance metrics for subset predicted as active
            active_subset = tpm.loc[spec_a_labels, :].copy()
            assert active_subset.shape[0] == spec_a_active, \
                'Selecting active subset failed: should be {}, is {}'.format(spec_a_active, active_subset.shape[0])
            r2_score_act = r2s(active_subset[cb], active_subset[ca])
            ktau_score_act = kendall_tau_scorer(active_subset[cb], active_subset[ca])

            num_tp = (np.logical_and(tpm[ca] >= 1, tpm[cb] >= 1)).sum()
            num_fp = (np.logical_and(tpm[ca] >= 1, tpm[cb] < 1)).sum()
            num_tn = (np.logical_and(tpm[ca] < 1, tpm[cb] < 1)).sum()
            num_fn = (np.logical_and(tpm[ca] < 1, tpm[cb] >= 1)).sum()

            dataset = tpm.loc[:, (name_a, ca, name_b, cb)].copy()
            dataset.columns = [name_a, norm_a, name_b, norm_b]
            dataset['{}_weights'.format(species_b)] = wt_vec
            dataset = dataset.merge(ranks.loc[:, (name_a, norm_a + '_rank')], on=name_a, how='outer', copy=True)
            dataset = dataset.merge(ranks.loc[:, (name_b, norm_b + '_rank')], on=name_b, how='outer', copy=True)

            metadata = {'num_orthologs': num_orthologs,
                        name_a + '_active': spec_a_active, name_a + '_inactive': spec_a_inactive,
                        name_b + '_inactive': spec_b_inactive, name_b + '_inactive': spec_b_inactive,
                        'perf_wt_acc': acc_score, 'perf_r2_all': r2_score_all, 'perf_ktau_all': ktau_score_all,
                        'perf_r2_active': r2_score_act, 'perf_ktau_active': ktau_score_act,
                        'perf_num_tp': num_tp, 'perf_num_fp': num_fp, 'perf_num_tn': num_tn, 'perf_num_fn': num_fn}
            metadata = pd.DataFrame.from_dict(metadata, orient='index')

            if tissue_match(tissue_a, tissue_b):
                base_group = '/pos/{}/{}/{}/{}'.format(species_a, species_b, norm_a, norm_b)
                with pd.HDFStore(outputfile, 'a', complib='blosc', complevel=9) as hdf:
                    hdf.put(os.path.join(base_group, 'data'), dataset, format='fixed')
                    hdf.put(os.path.join(base_group, 'metadata'), metadata, format='fixed')
                    hdf.flush()
            else:
                base_group = '/neg/{}/{}/{}/{}'.format(species_a, species_b, norm_a, norm_b)
                with pd.HDFStore(outputfile, 'a', complib='blosc', complevel=9) as hdf:
                    hdf.put(os.path.join(base_group, 'data'), dataset, format='fixed')
                    hdf.put(os.path.join(base_group, 'metadata'), metadata, format='fixed')
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
        with pd.HDFStore(args.orthofile, 'r') as hdf:
            orthologs = hdf[group]
        make_ortholog_pred(spec_a, tpm_a, ranks_a,
                           spec_b, tpm_b, ranks_b,
                           orthologs, args.outputfile)
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
