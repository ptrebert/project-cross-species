#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import argparse as argp
import traceback as trb
import collections as col
import logging as log
import json as js

import numpy as np
import pandas as pd

import scipy.stats as stats
import sklearn.metrics as sklm


def parse_command_line():
    """
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--assembly-a', '-asma', type=str, dest='assemblya')
    parser.add_argument('--assembly-b', '-asmb', type=str, dest='assemblyb')
    parser.add_argument('--exptab-a', '-expa', type=str, dest='exptaba')
    parser.add_argument('--exptab-b', '-expb', type=str, dest='exptabb')
    parser.add_argument('--exp-path', '-expp', type=str, default='/norm/tpm', dest='exppath')
    parser.add_argument('--orthologs', '-orth', type=str, dest='orthologs')
    parser.add_argument('--ortho-path', '-orthp', type=str, dest='orthopath',
                        default='/orthologs/subset/proteincoding')
    parser.add_argument('--debug', '-dbg', action='store_true', default=False, dest='debug')
    parser.add_argument('--json-out', '-jo', type=str, dest='jsonout')
    parser.add_argument('--hdf-out', '-ho', type=str, dest='hdfout')
    args = parser.parse_args()
    return args


def init_logger(debug):
    """
    :param debug:
    :return:
    """
    if debug:
        log.basicConfig(stream=sys.stderr, level=log.DEBUG, format='%(asctime)s - %(message)s')
    else:
        log.basicConfig(stream=sys.stderr, level=log.WARNING, format='%(asctime)s - %(message)s')
    logger = log.getLogger()
    logger.debug('Logger initialized')
    return logger


def load_orthologs(fpath, loadpath):
    """
    """

    with pd.HDFStore(fpath, 'r') as hdf:
        data = hdf[loadpath]
        md = hdf['/metadata']
        species1 = (md.loc[(md['key'] == 'species1'), 'value']).iloc[0]
        species2 = (md.loc[(md['key'] == 'species2'), 'value']).iloc[0]
        assert any([c.startswith(species1) for c in data.columns]),\
            'Species1 {} not in dataset: {}'.format(species1, data.columns)
        assert any([c.startswith(species2) for c in data.columns]),\
            'Species2 {} not in dataset: {}'.format(species2, data.columns)
    return data, species1, species2


def get_biosample_map():
    """
    :return:
    """
    samples = "K562 MEL GM12878 CH12 H1hESC ESE14 liver hepa".split()
    bsmap = dict()
    for s in range(0, len(samples), 2):
        assert samples[s] not in bsmap, 'Conflict: {}'.format(samples)
        bsmap[samples[s]] = samples[s+1]
        assert samples[s+1] not in bsmap, 'Conflict: {}'.format(samples)
        bsmap[samples[s+1]] = samples[s]
    return bsmap


def get_matching(mapping, col1, col2):
    """
    :param mapping:
    :param col1:
    :param col2:
    :return:
    """
    bio1 = col1.split('_')[2]
    bio2 = col2.split('_')[2]
    partner = mapping.get(bio1, 'n/a')
    return partner == bio2


def kendall_tau_scorer(x, y):
    """
    :param x:
    :param y:
    :return:
    """
    statistic, p_value = stats.kendalltau(x, y, nan_policy='raise')
    return statistic


def collect_metrics(expvals):
    """
    :param expvals:
    :return:
    """
    metrics = dict()
    expvals.insert(0, 'true_class', (expvals['true'] >= 1).astype(np.int8))
    expvals.insert(0, 'orth_class', (expvals['orth'] >= 1).astype(np.int8))
    class_counts = expvals['true_class'].value_counts(sort=True)  # assert that class 0 is at index 0
    weight_zero = 0.5 / class_counts.ix[0]
    weight_one = 0.5 / class_counts.ix[1]
    sample_weights = np.array([weight_zero if item < 1 else weight_one for item in expvals['true_class']],
                              dtype=np.float32)
    assert np.isclose(sample_weights.sum(), 1), 'Sample weights do not sum to 1: {}'.format(sample_weights.sum())
    metrics['accuracy'] = round(sklm.accuracy_score(expvals['true_class'].values,
                                                    expvals['orth_class'].values,
                                                    sample_weight=sample_weights), 5)
    metrics['ktau_all'] = round(kendall_tau_scorer(expvals['true'].values, expvals['orth']), 5)
    metrics['r2_all'] = round(sklm.r2_score(expvals['true'].values,
                                            expvals['orth'].values,
                                            sample_weight=sample_weights))
    subset = expvals.loc[expvals['orth_class'] == 1, :]
    metrics['ktau_sub'] = round(kendall_tau_scorer(subset['true'].values, subset['orth'].values), 5)
    metrics['r2_sub'] = round(sklm.r2_score(subset['true'].values, subset['orth'].values), 5)
    metadata = dict()
    metadata['true_active'] = expvals.loc[expvals['true_class'] == 1, :].index.tolist()
    metadata['num_true_active'] = len(metadata['true_active'])
    metadata['orth_active'] = subset.index.tolist()
    metadata['num_orth_active'] = len(metadata['orth_active'])
    return metadata, metrics


def evaluate_performance(expvals1, expvals2, select1, select2, orthmap, bsmap):
    """
    :param expvals1:
    :param expvals2:
    :param bsmap:
    :return:
    """
    collect_infos = {'pos': col.defaultdict(dict), 'neg': col.defaultdict(dict)}
    pos_pairs = None
    for c1 in sorted(expvals1.columns):
        for c2 in sorted(expvals2.columns):
            analog = get_matching(bsmap, c1, c2)
            group = 'pos' if analog else 'neg'
            source_vals = pd.DataFrame(expvals1[c1].values, columns=['values'], index=expvals1.index)
            source_vals.insert(0, select1, source_vals.index)
            source_vals = source_vals.reset_index(drop=True)
            modmap = orthmap.merge(source_vals, on=select1, how='outer', copy=True)
            orth_vals = modmap.groupby(select2).mean()
            true_vals = expvals2[c2]
            # make sure that genes/rows are properly aligned
            joined = pd.concat([true_vals, orth_vals], ignore_index=False, axis=1)
            assert joined.shape[1] == 2, 'Expected two columns after ortholog mapping: {}'.format(joined.shape)
            joined.columns = ['true', 'orth']
            metadata, metrics = collect_metrics(joined.copy())
            joined.columns = ['true_{}'.format(c2), 'orth_{}'.format(c1)]
            if group == 'pos':
                if pos_pairs is None:
                    pos_pairs = joined.copy()
                else:
                    pos_pairs = pd.concat([pos_pairs, joined], ignore_index=False, axis=1)
            else:
                # negative pairings will not be used to compare to, so save some space
                # by not keeping all names of active genes
                metadata['true_active'] = []
                metadata['orth_active'] = []
            metadata['select1'] = select1
            metadata['select2'] = select2
            key = c1 + '-' + c2
            assert key not in collect_infos, 'Duplicate detected: {}'.format(key)
            collect_infos[group][key]['metrics'] = metrics
            collect_infos[group][key]['metadata'] = metadata
    # remove redundant data - avoid having these in the DF should be implemented at some point...
    pos_pairs = pos_pairs.transpose()
    pos_pairs = pos_pairs.drop_duplicates(keep='first', inplace=False)
    pos_pairs = pos_pairs.transpose()
    return collect_infos, pos_pairs


def select_gene_subset(orthos, loadpath, selecta, exptaba, selectb, exptabb):
    """
    :param orthos:
    :param exptaba:
    :param exptabb:
    :return:
    """
    with pd.HDFStore(exptaba, 'r') as hdf:
        data_a = hdf[loadpath]
    with pd.HDFStore(exptabb, 'r') as hdf:
        data_b = hdf[loadpath]
    ortho_subset = orthos.loc[(orthos[selecta].isin(data_a.index) & orthos[selectb].isin(data_b.index)), :]
    assert (ortho_subset[selecta].isin(data_a.index)).all(), 'Some genes of species A not present in subset'
    assert (ortho_subset[selectb].isin(data_b.index)).all(), 'Some genes of species B not present in subset'
    data_a = data_a.loc[data_a.index.isin(ortho_subset[selecta]), :]
    data_b = data_b.loc[data_b.index.isin(ortho_subset[selectb]), :]
    return ortho_subset, data_a, data_b


def compare_groups(perfinfos):
    """
    :param perfinfos:
    :return:
    """
    all_pos = col.defaultdict(list)
    all_neg = col.defaultdict(list)
    for _, infos in perfinfos['pos'].items():
        for m, value in infos['metrics'].items():
            all_pos[m].append(value)
    for _, infos in perfinfos['neg'].items():
        for m, value in infos['metrics'].items():
            all_neg[m].append(value)
    group_data = {'test': 'KS_2samp', 'alternative': '2_sided'}
    for k, pos in all_pos.items():
        neg = all_neg[k]
        pv = float(stats.ks_2samp(pos, neg)[1])
        group_data[k] = pv
        group_data['num_pos'] = len(pos)
        group_data['num_neg'] = len(neg)
    return group_data


def main():
    """
    :return:
    """
    args = parse_command_line()
    logger = init_logger(args.debug)
    bsmap = get_biosample_map()
    orthos, spec1, spec2 = load_orthologs(args.orthologs, args.orthopath)
    logger.debug('Orthologs loaded: {}'.format(orthos.shape))
    logger.debug('Matching {} to {} and {} to {}'.format(spec1, args.assemblya, spec2, args.assemblyb))
    selecta = '{}_gene_name'.format(spec1)
    selectb = '{}_gene_name'.format(spec2)
    orthos, expa, expb = select_gene_subset(orthos, args.exppath, selecta, args.exptaba, selectb, args.exptabb)
    logger.debug('Ortholog subset selected: {}'.format(orthos.shape))
    logger.debug('Orthologs species A: {}'.format(expa.shape))
    logger.debug('Orthologs species B: {}'.format(expb.shape))
    perf_AB, exp_AB = evaluate_performance(expa, expb, selecta, selectb, orthos[[selecta, selectb]], bsmap)
    logger.debug('Performance evaluation complete: from {} to {}'.format(spec1, spec2))
    perf_BA, exp_BA = evaluate_performance(expb, expa, selectb, selecta, orthos[[selecta, selectb]], bsmap)
    logger.debug('Performance evaluation complete: from {} to {}'.format(spec2, spec1))
    group_AB = compare_groups(perf_AB)
    group_BA = compare_groups(perf_BA)
    logger.debug('Group comparisons complete')
    dump_object = {'species_A': spec1, 'species_B': spec2, 'assembly_A': args.assemblya, 'assembly_B': args.assemblyb,
                   'group_comparison': {'AB': group_AB, 'BA': group_BA}, 'performance': {'AB': perf_AB, 'BA': perf_BA}}
    with open(args.jsonout, 'w') as dump:
        js.dump(dump_object, dump, sort_keys=True, ensure_ascii=True, indent=1)
    logger.debug('JSON dumped')
    with pd.HDFStore(args.hdfout, 'w', complevel=9, complib='blosc') as hdf:
        hdf.put('AB/{}/{}'.format(args.assemblya, args.assemblyb), exp_AB, format='fixed')
        hdf.put('BA/{}/{}'.format(args.assemblyb, args.assemblya), exp_BA, format='fixed')
        hdf.flush()
    logger.debug('HDF dumped')
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(err))
        sys.exit(1)
    else:
        sys.exit(0)
