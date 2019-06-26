#!/usr/bin/env python3
# coding=utf-8

__author__ = "Peter Ebert"
__copyright__ = "Copyright (C) 2019 Peter Ebert"
__license__ = "GPLv3"

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as json
import re as re
import gzip as gzip
import io as io
import logging as logging
import logging.config as logconf

import pandas as pd


logger = logging.getLogger()


def parse_command_line():
    """
    :return:
    """
    script_full_path = os.path.realpath(__file__)
    script_dir = os.path.dirname(script_full_path)
    log_config_default_path = os.path.join(script_dir, 'configs', 'log_config.json')
    if not os.path.isfile(log_config_default_path):
        script_root = os.path.split(script_dir)[0]
        log_config_default_path = os.path.join(script_root, 'configs', 'log_config.json')
        if not os.path.isfile(log_config_default_path):
            log_config_default_path = ''

    parser = argp.ArgumentParser(add_help=True, allow_abbrev=False)
    parser.add_argument('--debug', '-d', action='store_true', dest='debug',
                        help='Print progress messages (by default: to stderr).')

    parser.add_argument('--use-logger', '-ul', type=str,
                        default='default', dest='use_logger',
                        help='Name of logger to use (default).')
    parser.add_argument('--log-config', '-lc', type=str,
                        default=log_config_default_path, dest='log_config',
                        help='Full path to JSON file containing '
                             'configuration parameters for the '
                             'loggers to use. A logger named "debug" '
                             'must be present in the configuration file.')
    parser.add_argument('--fasta-in', '-fin', type=str, dest='fastain')
    parser.add_argument('--transcripts', '-t', type=str, dest='transcripts')
    parser.add_argument('--fasta-compare', '-fc', type=str, dest='compare',
                        default='')
    parser.add_argument('--fasta-out', '-fout', type=str, dest='fastaout')

    args = parser.parse_args()
    return args


def init_logger(cli_args):
    """
    :param cli_args:
    :return:
    """
    if not os.path.isfile(cli_args.log_config):
        return
    with open(cli_args.log_config, 'r') as log_config:
        config = json.load(log_config)
    if 'debug' not in config['loggers']:
        raise ValueError('Logger named "debug" must be present '
                         'in log config JSON: {}'.format(cli_args.logconfig))
    logconf.dictConfig(config)
    global logger
    if cli_args.debug:
        logger = logging.getLogger('debug')
    else:
        logger = logging.getLogger(cli_args.use_logger)
    logger.debug('Logger initialized')
    return


def read_fasta_sequences(fasta_path):
    """
    :param fasta_path:
    :return:
    """
    transcript = re.compile('^>(ENS[A-Z0-9]+)')
    assert fasta_path.endswith('.gz'), 'Currently, only gzipped files supported'

    transcript_sequences = dict()
    current_transcript = None
    seq_buffer = io.StringIO()
    with gzip.open(fasta_path, 'rt') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if current_transcript is not None:
                    transcript_sequences[current_transcript] = seq_buffer
                    current_transcript = None
                    seq_buffer = io.StringIO

                mobj = transcript.search(line)
                if mobj is None:
                    if line.startswith('>ENS'):
                        raise ValueError('Potential ID miss: {}'.format(line.strip()))
                    continue
                else:
                    current_transcript = mobj.group(0).strip('>')
                    seq_buffer = io.StringIO()
            elif current_transcript is None:
                continue
            else:
                seq_buffer.write(line)
    if current_transcript is not None:
        transcript_sequences[current_transcript] = seq_buffer
    logger.debug('Read {} FASTA sequences from file {}'.format(len(transcript_sequences), fasta_path))
    return transcript_sequences


def compare_sequence_sets(input_seq, compare_seq):
    """
    :param input_seq:
    :param compare_seq:
    :return:
    """
    shared = set(input_seq.keys()).intersection(set(compare_seq.keys()))
    logger.debug('Shared transcripts between input and comparison set: {}'.format(len(shared)))
    identical = 0
    mismatch_length = 0
    mismatch_other = 0
    for s in shared:
        a = input_seq[s].getvalue()
        b = compare_seq[s].getvalue()
        if a == b:
            identical += 1
        elif len(a) != len(b):
            mismatch_length += 1
        else:
            mismatch_other += 1
    logger.debug('Identical sequences: {}'.format(identical))
    logger.debug('Length mismatch sequences: {}'.format(mismatch_length))
    logger.debug('Other mismatches: {}'.format(mismatch_other))
    return


def dump_fasta_sequences(sequences, transcripts, output_file):
    """
    :param sequences:
    :param transcripts:
    :param output_file:
    :return:
    """
    output_dir = os.path.abspath(os.path.dirname(output_file))
    os.makedirs(output_dir, exist_ok=True)

    with open(output_file, 'w') as fasta:
        for transcript in transcripts:
            _ = fasta.write('>' + transcript + '\n')
            _ = fasta.write(sequences[transcript].getvalue())
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    init_logger(args)
    logger.debug('Reading input FASTA file')
    input_seq = read_fasta_sequences(args.fastain)
    transcripts = pd.read_csv(args.transcripts, sep='\t')
    logger.debug('Read {} transcripts from file {}'.format(transcripts.shape[0], args.transcripts))
    transcripts = set(transcripts['name'].values)
    shared = transcripts.intersection(set(input_seq.keys()))
    logger.debug('{} transcripts shared between FASTA and annotation file'.format(len(shared)))
    # this should never be the case - all Ensembl
    assert len(shared) == len(transcripts), \
        'Transcript ID in annotation but not in FASTA file'
    if args.compare:
        logger.debug('Reading sequences for comparison')
        compare_seq = read_fasta_sequences(args.compare)
        compare_sequence_sets(input_seq, compare_seq)
    logger.debug('Writing transcript sequences...')
    dump_fasta_sequences(input_seq, shared, args.fastaout)
    logger.debug('Writing complete')

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
