#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as json
import logging as logging
import logging.config as logconf
import gzip as gzip
import io as io
import re as re
import itertools as itt
import collections as col

import pandas as pd


logger = logging.getLogger()


NUCLEOTIDES = list('ACGTNacgtn')


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
    parser.add_argument('--fasta-in', '-fin', type=str, dest='input', required=True)
    parser.add_argument('--fasta-out', '-fout', type=str, dest='output', required=True)
    parser.add_argument('--chromosome-sizes', '-csz', type=str, dest='chromsizes', required=True)
    parser.add_argument('--chromosome-regions', '-cre', type=str, dest='chromregions', required=True)
    parser.add_argument('--autosome-sizes', '-asz', type=str, dest='autosomes', default='')
    parser.add_argument('--autosome-regions', '-are', type=str, dest='autoregions', default='')
    parser.add_argument('--genome-metrics', '-met', type=str, dest='metrics', default='')

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


def dump_genome_sequence(chrom_order, chrom_sequence, output_file):
    """
    :param chrom_order:
    :param chrom_sequence:
    :param output_file:
    :return:
    """
    os.makedirs(os.path.abspath(os.path.dirname(output_file)), exist_ok=True)
    with open(output_file, 'w') as fasta:
        for chrom_name in chrom_order:
            _ = fasta.write('>' + chrom_name + '\n')
            _ = fasta.write(chrom_sequence[chrom_name].getvalue())
    return


def dump_chromosome_size_table(chrom_order, chrom_sizes, output_file, regions=False):
    """
    :param chrom_order:
    :param chrom_sizes:
    :param output_file:
    :param regions:
    :return:
    """
    os.makedirs(os.path.abspath(os.path.dirname(output_file)), exist_ok=True)
    with open(output_file, 'w') as table:
        for chrom_name in chrom_order:
            if regions:
                _ = table.write(chrom_name + '\t0\t' + str(chrom_sizes[chrom_name]) + '\n')
            else:
                _ = table.write(chrom_name + '\t' + str(chrom_sizes[chrom_name]) + '\n')
    return


def dump_metrics_table(chrom_order, chrom_metrics, output_file):
    """
    :param chrom_order:
    :param chrom_metrics:
    :param output_file:
    :return:
    """
    os.makedirs(os.path.abspath(os.path.dirname(output_file)), exist_ok=True)
    with open(output_file, 'w') as table:
        metrics = pd.DataFrame.from_records([chrom_metrics[chrom] for chrom in chrom_order],
                                            columns=NUCLEOTIDES + ['other'], index=chrom_order)
        metrics.fillna(0, inplace=True)
        metrics = metrics.astype('int64')
        gw_row = pd.DataFrame([metrics.sum(axis=0)],
                              columns=NUCLEOTIDES + ['other'],
                              index=['genome'])

        metrics = metrics.append(gw_row, ignore_index=False, sort=False)
        metrics['effective_size'] = metrics.sum(axis=1).astype('int64')
        metrics.loc[:, 'effective_size'] -= (metrics.loc[:, 'N'] + metrics.loc[:, 'n'] + metrics.loc[:, 'other'])
        metrics.to_csv(table, sep='\t', index=True, index_label='region')
    return


def read_gzipped_fasta(input_path):
    """
    :param input_path:
    :return:
    """
    chromosome_sizes = dict()
    chromosome_sequences = dict()
    current_chrom = None
    chrom_size = 0
    known_size = 0
    seq_buffer = io.StringIO()
    skip = False
    with gzip.open(input_path, 'rt') as fasta:
        for line in fasta:
            if line.startswith('>'):
                parts = line.strip().split()
                if parts[-1] != 'REF':
                    logger.debug('Skipping record {}'.format(parts[0]))
                    skip = True
                    continue
                skip = False
                if current_chrom is not None:
                    try:
                        assert known_size == chrom_size, \
                            'Chromosome size mismatch: {} vs {}'.format(chrom_size, known_size)
                    except AssertionError:
                        if not (current_chrom == 'chrY' and known_size == 56887902):
                            raise
                    chromosome_sizes[current_chrom] = chrom_size
                    chromosome_sequences[current_chrom] = seq_buffer

                current_chrom = parts[0].strip('>')
                seq_buffer = io.StringIO()
                chrom_size = 0
                known_size = int(parts[2].split(':')[-2])
                if 'chromosome' in parts[1]:
                    current_chrom = 'chr' + current_chrom
                elif 'primary_assembly' in parts[1]:
                    if 'EquCab' in parts[2]:
                        if not parts[0].startswith('>PJAA'):
                            current_chrom = 'chr' + current_chrom
                    elif 'ARS-UCD' in parts[2]:
                        if not parts[0].startswith('>NKL'):
                            current_chrom = 'chr' + current_chrom
                    else:
                        pass
                else:
                    pass
                logger.debug('Reading chromosome {}'.format(current_chrom))
            elif line.strip() and not skip:
                seq_buffer.write(line)
                chrom_size += len(line.strip())
            else:
                continue
    chromosome_sizes[current_chrom] = chrom_size
    chromosome_sequences[current_chrom] = seq_buffer
    return chromosome_sizes, chromosome_sequences


def compute_sequence_metrics(chrom_sequences):
    """
    :param chrom_sequences:
    :return:
    """
    counts_per_chrom = dict()
    for chrom, seqbuffer in chrom_sequences.items():
        clean_count = col.Counter()
        base_count = col.Counter(seqbuffer.getvalue())
        for base, abundance in base_count.items():
            if base in NUCLEOTIDES:
                clean_count[base] = abundance
            elif base in ['\n', '\t', '.', '-', '_']:
                continue
            else:
                clean_count['other'] += abundance
        counts_per_chrom[chrom] = clean_count
    return counts_per_chrom


def main():
    """
    :return:
    """
    args = parse_command_line()
    init_logger(args)
    logger.debug('Reading input FASTA file')
    sizes, sequences = read_gzipped_fasta(args.input)
    if args.metrics:
        logger.debug('Computing genome metrics')
        metrics = compute_sequence_metrics(sequences)
    else:
        metrics = None

    sort_order_main = []
    other = []
    main_ordering = {'X': 100., 'Y': 200., 'MT': 300.,
                     'Z': 100., 'W': 200.,
                     '2A': 2.3, '2B': 2.6,
                     '2a': 2.3, '2b': 2.6}
    # this is for cat and platypus
    for pos, (letter, number) in enumerate(itt.product(['A', 'B', 'C', 'D', 'E', 'F', 'X'], range(1, 7), repeat=1)):
        main_ordering[letter + str(number)] = float(pos)
    for chrom_name in sizes.keys():
        try:
            order_number = float(chrom_name.strip('chr'))
            sort_order_main.append((order_number, chrom_name))
        except ValueError:
            main_chromosome = re.match('^chr[0-9A-Z]+$', chrom_name) is not None
            if main_chromosome and not any([x in chrom_name for x in ['Un', 'random']]):
                order_number = main_ordering[chrom_name.strip('chr')]
                sort_order_main.append((order_number, chrom_name))
            else:
                other.append(chrom_name)
    logger.debug('Established chromosome sort order')
    
    sort_order = [c for n, c in sorted(sort_order_main)] + sorted(other)
    logger.debug('Writing genome sequence')
    dump_genome_sequence(sort_order, sequences, args.output)
    logger.debug('Writing chromosome size table')
    dump_chromosome_size_table(sort_order, sizes, args.chromsizes)
    logger.debug('Writing chromosome region BED file')
    dump_chromosome_size_table(sort_order, sizes, args.chromregions, True)

    if args.autosomes or args.autoregions:
        autosome_select = re.compile('chr[A-F0-9]{1,2}$')
        autosomes = []
        for c in sort_order:
            if autosome_select.match(c) is not None:
                autosomes.append(c)
        logger.debug('Selected {} autosomes'.format(len(autosomes)))
        if args.autosomes:
            logger.debug('Writing autosome size table')
            dump_chromosome_size_table(autosomes, sizes, args.autosomes)
        if args.autoregions:
            logger.debug('Writing autosome region table')
            dump_chromosome_size_table(autosomes, sizes, args.autoregions, True)

    if metrics is not None:
        logger.debug('Writing metrics file')
        dump_metrics_table(sort_order, metrics, args.metrics)
    logger.debug('Done')
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
