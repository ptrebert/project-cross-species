#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import pandas as pd
import argparse as argp


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()

    parser.add_argument('--exp-files', '-ef', type=str, nargs='+', dest='expfiles')
    parser.add_argument('--ortho-file', '-orth', type=str, dest='orthofile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')

    args = parser.parse_args()
    return args


def load_orthologs(fpath):
    """
    :param fpath:
    :return:
    """

    with pd.HDFStore(fpath, 'r') as hdf:
        dataset = hdf['/shared/subset']
    indexer = dataset.loc[dataset['chrom'].isin(['chrX', 'chrY']), 'og_id']
    dataset = dataset.loc[~dataset['og_id'].isin(indexer), :]
    return dataset


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


def load_expression_data(fpath, orthologs):
    """
    :param fpath:
    :param orthologs:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        tpm_data = hdf['/norm/tpm']
        rank_data = hdf['/norm/ranks']
    indexer = orthologs['gene_name'].isin(tpm_data.index)
    spec_orth = orthologs.loc[indexer, :]
    species = spec_orth['species'].unique()
    assert species.size == 1, 'More than one species: {}'.format(species)
    species = species[0]
    indexer = tpm_data.index.isin(spec_orth['gene_name'])
    tpm_data = tpm_data.loc[indexer, :]
    rank_data = rank_data.loc[indexer, :]
    tpm_data['gene_name'] = tpm_data.index
    tpm_data = tpm_data.reset_index(drop=True)
    rank_data['gene_name'] = rank_data.index
    rank_data = rank_data.reset_index(drop=True)
    return species, tpm_data, rank_data, spec_orth


def process_species_data(species, tpm, rank, orthologs):
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
    tpm_balanced = tpm_orth.loc[tpm_orth['group_balanced'] == 1, select_cols]
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
    return {'strict': tpm_strict, 'balanced': tpm_balanced, 'all': tpm_all},\
           {'strict': rank_strict, 'balanced': rank_balanced, 'all': rank_all}


def main():
    """
    :return:
    """
    args = parse_command_line()
    orths = load_orthologs(args.orthofile)
    df_tpm = {'strict': None, 'balanced': None, 'all': None}
    df_rank = {'strict': None, 'balanced': None, 'all': None}
    for ef in args.expfiles:
        spec, tpm, rank, spec_orth = load_expression_data(ef, orths)
        tpms, ranks = process_species_data(spec, tpm, rank, spec_orth)
        for k, df in df_tpm.items():
            dset_df = tpms[k]
            if df is None:
                df = dset_df
            else:
                df = pd.concat([df, dset_df], axis=1, join='outer', ignore_index=False)
            df_tpm[k] = df
        for k, df in df_rank.items():
            dset_df = ranks[k]
            if df is None:
                df = dset_df
            else:
                df = pd.concat([df, dset_df], axis=1, join='outer', ignore_index=False)
            df_rank[k] = df
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        for k, df in df_tpm.items():
            hdf.put('/tpm/{}'.format(k), df, format='fixed')
        hdf.flush()
        for k, df in df_rank.items():
            hdf.put('/rank/{}'.format(k), df, format='fixed')
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
