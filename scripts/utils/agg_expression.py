#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp

import numpy as np
import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, nargs='+', dest='inputfiles')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    col_names = ['chrom', 'start', 'end', 'name', 'tpm', 'strand', 'symbol']
    args = parse_command_line()
    agg_df = None
    species = ''
    for inpf in args.inputfiles:
        df = pd.read_csv(inpf, header=0, delimiter='\t', names=col_names,
                         usecols=[3, 4, 6], dtype={'name': str, 'tpm': np.float64, 'symbol': str})
        df = df.set_index(['name', 'symbol'], drop=True)
        assert df.shape[1] == 1, 'Re-creating index failed: {}'.format(df.shape)
        ds_name = os.path.basename(inpf).split('_')
        assert not species or ds_name[1] == species, 'Mismatch between datasets: {}'.format(inpf)
        species = ds_name[1]
        ds_name = '_'.join([ds_name[0], ds_name[2]])
        df.columns = [ds_name]
        if agg_df is None:
            agg_df = df
        else:
            agg_df = pd.concat([agg_df, df], axis=1, join='inner')
    with pd.HDFStore(args.outputfile, 'w', complib='blosc', complevel=9) as hdf:
        hdf.put('/{}/TPM'.format(species), agg_df, format='fixed')
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
