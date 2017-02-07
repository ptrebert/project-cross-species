#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp

import numpy as np
import scipy.stats as stats
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles', required=True)
    parser.add_argument('--agg-output', '-ao', type=str, dest='aggoutput')
    parser.add_argument('--split-output', '-so', type=str, dest='splitoutput')
    parser.add_argument('--gene-model', '-gm', type=str, dest='genemodel')
    args = parser.parse_args()
    return args


def nonzero_qnorm(mat):
    """
    :param mat:
    :return:
    """
    add_zero_columns = False
    col_idx = mat.sum(axis=0) > 0
    if np.sum(col_idx) != mat.shape[1]:
        # at least one all-zero column
        add_zero_columns = True
        mat = mat[:, col_idx]
    ranks = np.zeros_like(mat)
    for row_idx in range(mat.shape[0]):
        ranks[row_idx, :] = stats.rankdata(mat[row_idx, :], method='dense')
    mat.sort(axis=1)
    col_means = np.unique(mat.mean(axis=0))
    mean_ranks = stats.rankdata(col_means, method='dense')
    for row_idx in range(ranks.shape[0]):
        indices = np.digitize(ranks[row_idx, :], mean_ranks, right=True)
        replace_vals = col_means[indices]
        ranks[row_idx, :] = replace_vals
    norm_mat = ranks
    if add_zero_columns:
        add_zeros = np.zeros((ranks.shape[0], col_idx.size))
        col_indices = np.arange(col_idx.size)[col_idx]
        add_zeros[:, col_indices] = ranks[:]
        norm_mat = add_zeros
    return norm_mat


def normalize_tpm_data(dataset):
    """
    :param dataset:
    :return:
    """
    dataset = dataset.transpose()
    gene_columns = dataset.columns
    sample_rows = dataset.index
    norm_dataset = nonzero_qnorm(dataset.values)
    norm_dataset = pd.DataFrame(norm_dataset, columns=gene_columns, index=sample_rows, dtype=np.int32)
    norm_dataset = norm_dataset.transpose()
    return norm_dataset


def read_salmon_files(quantfiles):
    """
    :param quantfiles:
    :return:
    """
    col_names = ['name', 'length', 'effective_length', 'tpm', 'num_reads']
    use_cols = [0, 3]
    data_types = {'name': str, 'tpm': np.float32}
    collect_data = []
    num_rows = 0
    for qf in quantfiles:
        sample_id = os.path.split(os.path.split(qf)[0])[1]
        assert '_mRNA' in sample_id and sample_id.startswith('T'), \
               'Could not extract sample ID: {}'.format(sample_id)
        this_dataset = pd.read_csv(qf, header=0, delimiter='\t', names=col_names,
                                   usecols=use_cols, dtype=data_types, index_col=0,
                                   skip_blank_lines=True)
        if num_rows == 0:
            num_rows = this_dataset.shape[0]
        else:
            assert num_rows == this_dataset.shape[0], \
                    'Different number of genes in samples - ' \
                    'expected {}, got {} from {}'.format(num_rows, this_dataset.shape[0], sample_id)
        collect_data.append((sample_id, this_dataset))
    collect_data = sorted(collect_data, key=lambda x: x[0])
    sample_names = [t[0] for t in collect_data]
    full_dataset = pd.concat([t[1] for t in collect_data], ignore_index=False, axis=1)
    full_dataset.columns = sample_names
    assert full_dataset.shape[0] == num_rows,\
        'Full dataset contains {} rows, expected {}'.format(full_dataset.shape[0], num_rows)
    full_dataset.sort_index(axis=0, inplace=True)
    return full_dataset


def main():
    """
    :return:
    """
    args = parse_command_line()
    raw_dataset = read_salmon_files(args.inputfiles)
    norm_dataset = normalize_tpm_data(raw_dataset.copy(deep=True))
    assert raw_dataset.shape == norm_dataset.shape,\
        'Dataset size differs after normalization: before {} - after {}'.format(raw_dataset.shape, norm_dataset.shape)
    assert raw_dataset.columns == norm_dataset.columns,\
        'Column index mismatch: {} vs {}'.format(raw_dataset.columns, norm_dataset.columns)
    raw_ranks = raw_dataset.rank(axis=0, method='dense')
    norm_ranks = norm_dataset.rank(axis=0, method='dense')
    with pd.HDFStore(args.aggoutput, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('/raw/tpm', raw_dataset, format='fixed')
        hdf.put('/norm/tpm', norm_dataset, format='fixed')
        hdf.put('/raw/ranks', raw_ranks, format='fixed')
        hdf.put('/norm/ranks', norm_ranks, format='fixed')
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
