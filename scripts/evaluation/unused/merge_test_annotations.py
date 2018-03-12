#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import collections as col
import json as js


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
    args = parse_command_line()
    run_collector = col.defaultdict(list)
    run_pairs = col.Counter()
    for fp in args.inputfiles:
        new_pairs = col.Counter()
        data = js.load(open(fp, 'r'))
        for gid, runs in data.items():
            for idx, run in enumerate(sorted(runs, key=lambda d: d['run_file']), start=1):
                new_pairs[(run['run_spec_match'], run['run_test_type'])] += 1
                run['run_number'] = idx
                if 'Status' in run['setting']:
                    run['setting_number'] = 1
                elif 'TPM' in run['setting']:
                    run['setting_number'] = 2
                elif 'Est.' in run['setting']:
                    run['setting_number'] = 3
                else:
                    raise ValueError('Unexpected setting of run {}: {}'.format(run['run_file']), run['setting'])
                run_collector[gid].append(run)
        if run_pairs:
            assert run_pairs == new_pairs, 'Different run types recorded: {} - {}'.format(run_pairs, new_pairs)
        else:
            run_pairs = new_pairs
    serialize = sorted(run_pairs.most_common(), key=lambda x: x[1], reverse=True)
    with open(args.outputfile, 'w') as outf:
        js.dump({'run_pairs': serialize, 'test_runs': run_collector}, outf, indent=1, sort_keys=True)
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
