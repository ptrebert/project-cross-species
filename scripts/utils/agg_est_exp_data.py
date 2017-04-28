#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp

import pandas as pd


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, dest='inputfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    parser.add_argument('--queries', '-q', type=str, nargs='+',
                        default=['human', 'mouse', 'cow'], dest='queries')

    args = parser.parse_args()
    return args


def collect_query_data(fpath, queries):
    """
    :param fpath:
    :param queries:
    :return:
    """
    # ('/', 'testing', 'level', 'all', trg_spec, qry_spec, gid, datafile, modelfile)
    data_collect = dict()
    with pd.HDFStore(fpath, 'r') as hdf:
        for k in hdf.keys():
            if k.startswith('/testing') and k.endswith('/data'):
                components = k.split('/')
                out_type = components[2]
                query = components[5]
                if out_type == 'level' and query in queries:
                    subset = components[3]
                    target = components[4]
                    data = hdf[k]
                    data.drop([c for c in data.columns if c.startswith('ortho_')], axis=1, inplace=True)
                    datafile_parts = components[7].split('_')
                    exp_name = datafile_parts[2] + '_' + datafile_parts[4]
                    modelfile_parts = components[8].split('_')
                    mod_name = exp_name + '_' + '_'.join(modelfile_parts[:3])
                    if data.columns[0] == 'true':
                        data.columns = [exp_name, mod_name]
                    else:
                        data.columns = [mod_name, exp_name]
                    key = (target, query, subset)
                    if key in data_collect:
                        collected = data_collect[key]
                        if exp_name in collected:
                            data.drop(exp_name, axis=1, inplace=True)
                        assert mod_name not in collected.columns, 'Model duplicate under path {}'.format(k)
                        collected = pd.concat([collected, data], axis=1, join='outer', ignore_index=False)
                        data_collect[key] = collected
                    else:
                        data_collect[key] = data
    for k, v in data_collect.items():
        if k[2] == 'active':
            # since different genes are predicted to
            # be active, this can contain NaN
            v.dropna(axis=0, how='any', inplace=True)
            assert not v.empty, 'Dataframe empty after dropping NA rows: {}'.format(k)
    return data_collect


def main():
    """
    :return:
    """
    args = parse_command_line()
    query_data = collect_query_data(args.inputfile, args.queries)

    with pd.HDFStore(args.outputfile, 'w', complevel=9, complib='blosc') as hdf:
        for k, v in query_data.items():
            hdf.put(os.path.join(k[2], k[0], k[1]), v, format='fixed')
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
