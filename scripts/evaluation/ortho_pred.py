#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import argparse as argp
import traceback as trb
import collections as col
import functools as fnt
import logging as log

import pandas as pd

import scipy.stats as stats
from sklearn.metrics import accuracy_score, r2_score, make_scorer


def parse_command_line():
    """
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--assembly-a', '-spa', type=str, dest='assemblya')
    parser.add_argument('--assembly-b', '-spb', type=str, dest='assemblyb')
    parser.add_argument('--exptab-a', '-expa', type=str, dest='exptaba')
    parser.add_argument('--exptab-b', '-expb', type=str, dest='exptabb')
    parser.add_argument('--exp-path', '-expp', type=str, default='/norm/tpm', dest='exppath')
    parser.add_argument('--orthologs', '-orth', type=str, dest='orthologs')
    parser.add_argument('--ortho-path', '-orthp', type=str, dest='orthopath',
                        default='/orthologs/subset/proteincoding')
    parser.add_argument('--bio-samples', '-bsam', type=str, dest='biosamples')
    parser.add_argument('--debug', '-dbg', action='store_true', default=False, dest='debug')
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


def read_biosample_map(fpath):
    """
    :param fpath:
    :return:
    """
    bsmap = dict()
    with open(fpath, 'r') as infile:
        for line in infile:
            if not line.strip():
                continue
            p1, p2 = line.strip().split()
            assert p1 not in bsmap, 'Sample {} already in mapping'.format(p1)
            bsmap[p1] = p2
            assert p2 not in bsmap, 'Sample {} already in mapping'.format(p2)
            bsmap[p2] = p1
    return bsmap


def eval_performance(truevals, orthmap, mode):
    """
    """
    true_genes = set(truevals['name'].tolist())
    orth_genes = set(orthmap['name'].tolist())
    genes = orth_genes.intersection(true_genes)
    subset_true = truevals.loc[truevals['name'].isin(genes), :]
    # subset normalized
    subset_orth = orthmap.loc[orthmap['name'].isin(genes), :]
    joined = pd.merge(subset_true, subset_orth, on='name', how='outer', copy=True)
    true_cols = sorted([c for c in subset_true.columns if c not in ['name', 'T010_MEL']])
    orth_cols = sorted([c for c in subset_orth.columns if c not in ['name', 'T010_MEL']])
    assert true_cols != orth_cols, 'Wrong columns selected'
    collect_class = col.defaultdict(list)
    collect_reg = col.defaultdict(list)
    collect_exp = col.defaultdict(list)
    use_comparison = 1 if mode == 'pos' else 0
    for row in joined.itertuples():
        for tc in true_cols:
            for oc in orth_cols:
                if use_comparison == 1:
                    use = get_compat_info(tc, oc)
                    if use != use_comparison:
                        continue
                key = tc + '-' + oc
                tc_val_raw = getattr(row, tc)
                tc_val = 0 if tc_val_raw < 1 else round(tc_val_raw)
                oc_val_raw = getattr(row, oc)
                oc_val = 0 if oc_val_raw < 1 else round(oc_val_raw)
                true_label = 0 if tc_val < 1 else 1
                pred_label = 0 if oc_val < 1 else 1
                collect_class[key].append((row.name, true_label, pred_label))
                collect_reg[key].append((row.name, tc_val_raw, oc_val_raw))
                lo, hi = round(tc_val * 0.8), round(tc_val * 1.2)
                if ((oc_val >= 1 and tc_val >= 1) and
                        (lo <= oc_val <= hi)):
                    collect_exp[key].append((row.name, tc_val_raw, oc_val_raw))
    return collect_class, collect_reg, collect_exp
      

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


def evaluate_performance(expvals1, expvals2, select1, select2, orthmap, bsmap):
    """
    :param expvals1:
    :param expvals2:
    :param bsmap:
    :return:
    """
    collect_r2 = {'pos': col.defaultdict(list), 'neg': col.defaultdict(list)}
    collect_kt = {'pos': col.defaultdict(list), 'neg': col.defaultdict(list)}
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
            # make sure that are arranged properly
            joined = pd.concat([true_vals, orth_vals], ignore_index=False, axis=1)
            joined.columns = ['true', 'orth']
            assert joined.shape[1] == 2, 'Expected two columns after ortholog mapping: {}'.format(joined.shape)
            r2_stat = r2_score(joined['true'].values, joined['orth'].values)
            kt_stat = kendall_tau_scorer(joined['true'].values, joined['orth'].values)
            collect_kt[group][(c1, c2)].append(kt_stat)
            collect_r2[group][(c1, c2)].append(r2_stat)
    all_pos = []
    all_neg = []
    for _, vals in collect_kt['pos'].items():
        all_pos.extend(vals)
    for _, vals in collect_kt['neg'].items():
        all_neg.extend(vals)
    print(len(all_pos))
    print(len(all_neg))
    print(stats.ks_2samp(all_neg, all_pos))
    return


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


def main():
    """
    :return:
    """
    args = parse_command_line()
    logger = init_logger(args.debug)
    bsmap = read_biosample_map(args.biosamples)
    orthos, spec1, spec2 = load_orthologs(args.orthologs, args.orthopath)
    logger.debug('Orthologs loaded: {}'.format(orthos.shape))
    logger.debug('Matching {} to {} and {} to {}'.format(spec1, args.assemblya, spec2, args.assemblyb))
    selecta = '{}_gene_name'.format(spec1)
    selectb = '{}_gene_name'.format(spec2)
    orthos, expa, expb = select_gene_subset(orthos, args.exppath, selecta, args.exptaba, selectb, args.exptabb)
    logger.debug('Ortholog subset selected: {}'.format(orthos.shape))
    logger.debug('Orthologs species A: {}'.format(expa.shape))
    logger.debug('Orthologs species B: {}'.format(expb.shape))
    evaluate_performance(expa, expb, selecta, selectb, orthos[[selecta, selectb]], bsmap)

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



"""
if __name__ == '__main__':
    orth = load_orthologs()
    run_mode = sys.argv[1]
    src_species, dest_species = sys.argv[2:]
    print(src_species, dest_species)
    classperf, regperf, consperf = find_ortho_cons_genes(src_species, dest_species, orth, run_mode)
    cons_genes = col.Counter()
    for key, values in consperf.items():
        this_genes = set([t[0] for t in values])
        print(len(this_genes))
        cons_genes.update(this_genes)
    for key, values in classperf.items():
        true = [t[1] for t in values]
        pred = [t[2] for t in values]
        print('Acc. {}: {}'.format(key, accuracy_score(true, pred)))
    for key, values in regperf.items():
        true = [t[1] for t in values]
        pred = [t[2] for t in values]
        print('R2 {}: {}'.format(key, r2_score(true, pred)))
    outpath = os.path.join(outdir, '{}_cons_orth_pm20pct.txt'.format(dest_species))
    with open(outpath, 'w') as outf:
        for gene, count in cons_genes.most_common():
            if count > 1:
                _ = outf.write(gene + '\n')
    #dump_cons_genes(expperf, outpath)
    sys.exit(0)
"""
