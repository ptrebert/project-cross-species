#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as js


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--train-file', '-ff', type=str, dest='trainfile')
    parser.add_argument('--output', '-o', type=str, dest='outputfile')
    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    args = parse_command_line()
    load_file = args.trainfile
    if load_file.endswith('.pck'):
        load_file = load_file.replace('.pck', '.json')
    class_priors = []
    with open(load_file, 'r') as dumped:
        train_infos = js.load(dumped)
        sample_names = train_infos['sample_info']['names']


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
