#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import gzip as gz
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
    parser.add_argument('--regime-output', '-ro', type=str, dest='regimeoutput')
    parser.add_argument('--model-assembly', '-ma', type=str, dest='modelassembly')
    parser.add_argument('--model-dir', '-md', type=str, dest='modeldir')
    args = parser.parse_args()
    return args


def select_gene_model(folder, assembly):
    """
    :param folder:
    :param assembly:
    :return:
    """
    all_models = os.listdir(folder)
    select_model = [f for f in all_models if assembly in f and f.endswith('.genes.bed.gz')]
    assert len(select_model) == 1, 'No/ambig. gene model selected: {}'.format(select_model)
    return os.path.join(folder, select_model[0])


def nonzero_qnorm(mat):
    """
    :param mat:
    :return:
    """
    assert np.isfinite(mat).all(), 'Matrix contains invalid (NaN, INF) values'
    mat = mat.round(decimals=3)
    add_zero_columns = False
    col_idx = np.array(mat.sum(axis=0) > 0, dtype=np.bool)
    non_zero_cols = np.sum(col_idx)
    if non_zero_cols == 0:
        return mat
    orig_shape = mat.shape
    if non_zero_cols != mat.shape[1]:
        # at least one all-zero column
        add_zero_columns = True
        mat = mat[:, col_idx]
    # determine data ranks
    data_ranks = pd.DataFrame(mat).rank(axis=1, method='min')
    data_ranks = data_ranks.astype(np.int32)
    uniq_data_ranks = np.unique(data_ranks.values.flatten())

    # get column means
    mat.sort(axis=1)
    col_means = mat.mean(axis=0)
    col_means = col_means.round(decimals=3).astype(np.float32)

    # ranks for column means
    col_mean_ranks = stats.rankdata(col_means, method='min')
    col_mean_ranks = col_mean_ranks.astype(np.int32)

    # lookup table for rank -> mean replacement
    replace_means = pd.Series(col_means, index=col_mean_ranks,
                              dtype=np.float32)
    # reduce to unique entries (i.e., group by index value = rank of mean value)
    replace_means = replace_means.groupby(by=replace_means.index).first()

    # Check if there are ranks missing
    # There are cases where data_ranks contains values that do
    # not appear in col_mean_ranks - apparently, the ranking of
    # floating point numbers is an issue here
    missing_ranks = set(uniq_data_ranks).difference(set(replace_means.index))
    if len(missing_ranks) > 0:
        missing_means = []
        missing_positions = []
        for mr in sorted(missing_ranks):
            lo = replace_means[replace_means.index < mr].idxmax()
            hi = replace_means[replace_means.index > mr].idxmin()
            lo_val = replace_means[lo]
            hi_val = replace_means[hi]
            new_val = np.round((lo_val + hi_val) / 2., 3)
            missing_means.append(new_val)
            missing_positions.append(mr)

        additional_means = pd.Series(missing_means, index=missing_positions, dtype=np.float32)
        replace_means = replace_means.append(additional_means)

    for row_idx in range(mat.shape[0]):
        # by constructions, this must give a valid replacement vector
        mat[row_idx, :] = (replace_means[data_ranks.iloc[row_idx, :]]).values

    if add_zero_columns:
        add_zeros = np.zeros(orig_shape, dtype=np.float32)
        col_indices = np.arange(col_idx.size)[col_idx]
        add_zeros[:, col_indices] = mat[:]
        mat = add_zeros
    # Note that this is ok since all zero matrices are returned
    # immediately in the beginning
    assert mat.max() > 0, 'Normalization failed: ends with all zero data matrix'
    assert np.isfinite(mat).all(), 'Result matrix contains invalid (NaN, INF) values'
    return mat


def normalize_tpm_data(dataset):
    """
    :param dataset:
    :return:
    """
    dataset = dataset.transpose()
    gene_columns = dataset.columns
    sample_rows = dataset.index
    norm_dataset = nonzero_qnorm(dataset.values)
    norm_dataset = pd.DataFrame(norm_dataset, columns=gene_columns, index=sample_rows, dtype=np.float32)
    norm_dataset = norm_dataset.transpose()
    return norm_dataset


def read_salmon_files(quantfiles, genes):
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
        new_dataset = pd.read_csv(qf, header=0, delimiter='\t', names=col_names,
                                  usecols=use_cols, dtype=data_types, index_col=0,
                                  skip_blank_lines=True)
        row_indexer = new_dataset.index.isin(genes)
        this_dataset = pd.DataFrame(new_dataset.loc[row_indexer, :])
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


def read_gene_annotation(fpath):
    """
    :param fpath:
    :return:
    """
    with gz.open(fpath, 'rt') as ann:
        header = ann.readline().strip().split('\t')
        header[0] = header[0].strip('#')
    genes = pd.read_csv(fpath, header=0, sep='\t', names=header,
                        skip_blank_lines=True, comment=None)
    select_chroms = '^chr[0-9A-F]+'  # not interested in gonosomes, and samples are of mixed karyotype
    select_rows = genes['chrom'].str.match(select_chroms, as_indexer=True)
    auto_genes = pd.DataFrame(genes.loc[select_rows, :])
    return auto_genes


def make_split_output(rawdat, rawrk, normdat, normrk, genes, outdir):
    """
    :param rawdat:
    :param rawrk:
    :param normdat:
    :param normrk:
    :param genes:
    :param outdir:
    :return:
    """
    sample_names = rawdat.columns
    header = genes.columns.tolist()
    header += ['tpm_raw', 'rank_raw', 'tpm_norm', 'rank_norm']
    write_cols = list(header)
    header[0] = '#' + header[0]
    out_shape = None
    for sn in sample_names:
        outfile = os.path.join(outdir, sn + '.genes.bed.gz')
        assert (rawdat[sn].index == normrk[sn].index).all(),\
            'Index mismatch (1)'
        assert (normdat[sn].index == rawrk[sn].index).all(),\
            'Index mismatch (2)'
        out_df = pd.concat([rawdat[sn], rawrk[sn], normdat[sn], normrk[sn]],
                           ignore_index=False, copy=True, axis=1)
        if out_shape is None:
            out_shape = out_df.shape
        else:
            assert out_shape == out_df.shape,\
                'Created datasets of different sizes for sample {}. ' \
                'Expected {} - got {}'.format(sn, out_shape, out_df.shape)
        out_df.columns = ['tpm_raw', 'rank_raw', 'tpm_norm', 'rank_norm']
        out_df.insert(0, 'name', out_df.index)
        out_df = out_df.merge(genes, how='outer', on='name')
        out_df = out_df.sort_values(by=['chrom', 'start', 'end'], axis=0,
                                    inplace=False, na_position='first')
        out_df.to_csv(outfile, sep='\t', header=header, columns=write_cols,
                      compression='gzip', line_terminator='\n', na_rep='N/A',
                      index=False)
    return


def make_regime_output(expdata, genes, assm, outdir):
    """
    :param expdata:
    :param genes:
    :param assm:
    :param outdir:
    :return:
    """
    # on/off regime
    stable_off = (expdata < 1).all(axis=1)
    stable_on = (expdata >= 1).all(axis=1)

    stable_on_off = np.logical_or(stable_off, stable_on)

    indicator = [stable_off, stable_on, stable_on_off]

    stable_low = (expdata < 5).all(axis=1)
    stable_high = (expdata >= 5).all(axis=1)

    stable_low_high = np.logical_or(stable_low, stable_high)

    indicator.extend([stable_low, stable_high, stable_low_high])

    unstable = np.logical_not(np.logical_or(stable_on_off, stable_low_high))

    indicator.append(unstable)

    labels = ['stable_off', 'stable_on', 'stable_on_off',
              'stable_low', 'stable_high', 'stable_low_high',
              'unstable']
    for idx, lab in zip(indicator, labels):
        selected_genes = expdata.loc[idx, :].index.tolist()
        if len(selected_genes) > 0:
            with open(os.path.join(outdir, '{}_{}.txt'.format(assm, lab)), 'w') as dump:
                _ = dump.write('\n'.join(selected_genes) + '\n')
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    gene_model = select_gene_model(args.modeldir, args.modelassembly)
    genes = read_gene_annotation(gene_model)
    raw_dataset = read_salmon_files(args.inputfiles, genes['name'])
    norm_dataset = normalize_tpm_data(raw_dataset.copy(deep=True))
    assert raw_dataset.shape == norm_dataset.shape,\
        'Dataset size differs after normalization: before {} - after {}'.format(raw_dataset.shape, norm_dataset.shape)
    assert (raw_dataset.columns == norm_dataset.columns).all(),\
        'Column index mismatch: {} vs {}'.format(raw_dataset.columns, norm_dataset.columns)
    # due to the changed q-norm implementation involving rounding and lower precision
    # floating point arithmetic, the following ranking can actually show
    # some minor variation (in an ideal scenario, they would be consistent)
    raw_ranks = raw_dataset.rank(axis=0, method='dense').astype(np.int32)
    norm_ranks = norm_dataset.rank(axis=0, method='dense').astype(np.int32)
    with pd.HDFStore(args.aggoutput, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('/raw/tpm', raw_dataset, format='fixed')
        hdf.put('/norm/tpm', norm_dataset, format='fixed')
        hdf.put('/raw/ranks', raw_ranks, format='fixed')
        hdf.put('/norm/ranks', norm_ranks, format='fixed')
        hdf.flush()
    make_split_output(raw_dataset, raw_ranks, norm_dataset,
                      norm_ranks, genes, args.splitoutput)
    make_regime_output(norm_dataset, genes, args.modelassembly, args.regimeoutput)
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
