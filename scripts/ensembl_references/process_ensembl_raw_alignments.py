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
import logging as logging
import logging.config as logconf
import gzip as gzip
import collections as col
import multiprocessing as mp


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
    parser.add_argument('--maf-folder', '-mf', type=str, dest='maffolder')
    parser.add_argument('--reference-chromosomes', '-rc', type=str, dest='refchroms')
    parser.add_argument('--target-chromosomes', '-tc', type=str, dest='trgchroms')
    parser.add_argument('--num-cpu', '-n', type=int, default=1, dest='numcpu')
    parser.add_argument('--out-block-counts', '-oc', type=str, dest='outcounts',
                        help='Path to output file holding counts of MAF alignments'
                             ' between chromosomes (block level).')

    parser.add_argument('--out-block-splits', '-os', type=str, dest='outsplits',
                        help='Path to output file holding coordinates of non-gapped'
                             ' alignments between target and reference')
    parser.add_argument('--keep-chromosome-splits', '-ks', action='store_true',
                        default=False, dest='keepsplits')

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


def read_chromosome_file(fpath):
    """
    :param fpath:
    :return:
    """
    chromosomes = []
    with open(fpath, 'r') as table:
        for line in table:
            chromosomes.append(line.split()[0])
    species = os.path.basename(fpath).split('.')[0]
    return species, chromosomes


def determine_maf_reference_chromosome(maf_files, keep_chromosomes):
    """
    Due to irregular file naming, safest option is to
    read first alignment block and determine reference
    chromosome - discard file if chromosome is not needed

    :param maf_files:
    :param keep_chromosomes:
    :return:
    """
    maf_chrom_sets = col.defaultdict(list)
    for maf in maf_files:
        filename = os.path.basename(maf)
        skip_file = False
        if 'readme' in maf.lower() or 'md5sum' in maf.lower():
            logger.debug('Skipping metadata file {}'.format(filename))
        elif maf.endswith('.maf') or maf.endswith('.maf.gz'):
            if maf.endswith('.maf.gz'):
                open_file = gzip.open
                open_mode = 'rt'
            else:
                open_file = open
                open_mode = 'r'
            ref_chrom = None
            with open_file(maf, open_mode) as alignments:
                for line in alignments:
                    if line.startswith('s '):
                        # ref species comes first
                        ref_chrom = line.split()[1].split('.')[1]
                        if not ref_chrom.startswith('chr'):
                            ref_chrom = 'chr' + ref_chrom
                        if ref_chrom not in keep_chromosomes:
                            skip_file = True
                        break
            if ref_chrom is None:
                logger.warning('No reference chromosome identified for file {} - skipping'.format(filename))
                continue
            if skip_file:
                logger.debug('Skipping file {} - reference chromosome {}'.format(filename, ref_chrom))
                continue
            logger.debug('Linked file {} to reference chromosome {}'.format(filename, ref_chrom))
            maf_chrom_sets[ref_chrom].append(maf)
        else:
            logger.warning('Skipping unrecognized file {}'.format(filename))
            continue
    return maf_chrom_sets


def split_alignment_block(ref_parts, trg_parts):
    """
    By manual examination, it seems like Ensembl
    MAF files list blocks/sequences actually as
    0-based; this function follows this observation,
    somewhat high risk for off-by-one

    :param ref_parts:
    :param trg_parts:
    :return:
    """
    ref_chrom = ref_parts[1].split('.')[1]
    trg_chrom = trg_parts[1].split('.')[1]

    ref_frw_start = int(ref_parts[2])
    ref_strand = ref_parts[4]
    assert ref_strand == '+', 'Inverse reference: {}'.format(ref_parts)

    trg_strand = trg_parts[4]
    if trg_strand == '+':
        trg_frw_start = int(trg_parts[2])
    else:
        # chrom_size - rev_start - length (#chars)
        trg_frw_start = int(trg_parts[5]) - int(trg_parts[2]) - int(trg_parts[3])

    ref_seq = ref_parts[-1].strip()
    assert ref_seq[0] != '-', 'Gap at beginning of sequence: {}'.format(ref_parts)
    trg_seq = trg_parts[-1].strip()
    assert trg_seq[0] != '-', 'Gap at beginning of sequence: {}'.format(trg_parts)

    ref_frw_end = ref_frw_start
    trg_frw_end = trg_frw_start

    block_splits = []

    for r, t in zip(ref_seq, trg_seq):
        if r == '-' or t == '-':
            if ref_frw_end > ref_frw_start:
                # we have a proper split
                assert trg_frw_end > trg_frw_start, 'Invalid split: {} - {}'.format(ref_frw_start, ref_parts)
                block_splits.append((ref_chrom, str(ref_frw_start), str(ref_frw_end), '+',
                                     trg_chrom, str(trg_frw_start), str(trg_frw_end), trg_strand))
                ref_frw_start = ref_frw_end
                trg_frw_start = trg_frw_end
            if r == '-' and t == '-':
                # nothing, just iterate over the gaps
                pass
            elif r == '-':  # and not t
                trg_frw_start += 1
                trg_frw_end = trg_frw_start  # no usable split
            elif t == '-':  # and not r
                ref_frw_start += 1
                ref_frw_end = ref_frw_start  # no usable split
            else:
                raise RuntimeError('Impossible condition: {} - {}'.format(ref_parts, trg_parts))
        else:
            ref_frw_end += 1
            trg_frw_end += 1

    if ref_frw_end > ref_frw_start:
        assert trg_frw_end > trg_frw_start, 'Invalid split: {} - {}'.format(ref_frw_start, ref_parts)
        block_splits.append((ref_chrom, str(ref_frw_start), str(ref_frw_end), '+',
                             trg_chrom, str(trg_frw_start), str(trg_frw_end), trg_strand))

    return block_splits


def check_alignment_block(ref_species, ref_chrom, ref_record, trg_species, trg_chroms, trg_record):
    """
    :param ref_species:
    :param ref_chrom:
    :param ref_record:
    :param trg_species:
    :param trg_chroms:
    :param trg_record:
    :return:
    """
    trg_parts = trg_record.split()
    trg_chrom = trg_parts[1].split('.')[1]
    if not trg_chrom.startswith('chr'):
        trg_chrom = 'chr' + trg_chrom
    if trg_chrom not in trg_chroms:
        return 0, None, None, None
    trg_parts[1] = trg_species + '.' + trg_chrom
    new_trg_record = ' '.join(trg_parts)

    ref_parts = ref_record.split()
    ref_parts[1] = ref_species + '.' + ref_chrom
    block_start = int(ref_parts[2])
    new_ref_record = ' '.join(ref_parts)
    new_block = 'a\n{}\n{}\n'.format(new_ref_record, new_trg_record)

    splits = split_alignment_block(ref_parts, trg_parts)

    return block_start, new_block, trg_chrom, splits


def process_maf_chromosome_set(params):
    """
    :param params:
    :return:
    """
    ref_species, ref_chrom, trg_species, trg_chromosomes, maf_files, maf_output, split_output = params

    alignment_blocks = []
    alignment_splits = []

    block_counter = col.Counter()

    open_file = open
    open_mode = 'r'
    if maf_files[0].endswith('.maf.gz'):
        open_file = gzip.open
        open_mode = 'rt'

    for maf in maf_files:

        ref = ''
        trg = ''
        with open_file(maf, open_mode) as alignments:
            for line in alignments:
                if line.startswith('s '):
                    if ref:
                        trg = line.strip()
                        assert ref and trg, 'Corrupt block - last line in file {}: {}'.format(maf, line.strip())
                        block_start, new_block, trg_chrom, block_splits = check_alignment_block(
                            ref_species, ref_chrom, ref,
                            trg_species, trg_chromosomes, trg
                        )
                        if block_start > 0:
                            block_counter[(ref_chrom, trg_chrom)] += 1
                            alignment_blocks.append((block_start, new_block))
                            alignment_splits.extend(block_splits)
                        ref = ''
                        trg = ''
                    else:
                        ref = line.strip()

    alignment_blocks = sorted(alignment_blocks)

    with gzip.open(maf_output, 'wt') as dump:
        for _, block in alignment_blocks:
            _ = dump.write(block + '\n')

    alignment_splits = sorted(alignment_splits, key=lambda x: int(x[1]))

    with gzip.open(split_output, 'wt') as dump:
        for pos, split in enumerate(alignment_splits, 1):
            split_name = '\tR-{}_{}_T-{}\t'.format(split[0], pos, split[4])
            named_split = '\t'.join(split[:4]) + split_name + '\t'.join(split[4:]) + '\n'
            _ = dump.write(named_split)

    return maf_output, split_output, block_counter, ref_chrom


def main():
    """
    :return:
    """
    args = parse_command_line()
    init_logger(args)
    ref_species, ref_chroms = read_chromosome_file(args.refchroms)
    trg_species, trg_chroms = read_chromosome_file(args.trgchroms)
    logger.debug('Loading MAF files from folder {}'.format(args.maffolder))
    all_files = sorted([os.path.join(args.maffolder, f) for f in os.listdir(args.maffolder)])
    maf_chrom_sets = determine_maf_reference_chromosome(all_files, ref_chroms)

    logger.debug('Building parameter list for processing')
    job_parameters = []
    output_dir = os.path.abspath(os.path.dirname(args.outcounts))
    os.makedirs(output_dir, exist_ok=True)
    for ref_chrom, maf_files in maf_chrom_sets.items():
        maf_output_file = '{}_vs_{}.{}.maf.gz'.format(ref_species, trg_species, ref_chrom)
        maf_output_path = os.path.join(output_dir, maf_output_file)
        split_output_file = '{}_vs_{}.{}.tsv.gz'.format(ref_species, trg_species, ref_chrom)
        split_output_path = os.path.join(output_dir, split_output_file)
        job_parameters.append((ref_species, ref_chrom, trg_species, trg_chroms,
                               sorted(maf_files), maf_output_path, split_output_path))

    logger.debug('Created parameter list of size {}'.format(len(job_parameters)))

    merged_counts = col.Counter()
    created_split_files = []

    with mp.Pool(min(args.numcpu, len(job_parameters))) as pool:
        logger.debug('Initialized worker pool')
        resit = pool.imap(process_maf_chromosome_set, job_parameters)
        for maf_out, split_out, block_counts, ref_chrom in resit:
            logger.debug('Finished preprocessing for: {}'.format(os.path.basename(maf_out)))
            merged_counts.update(block_counts)
            created_split_files.append((int(ref_chrom.strip('chr')), split_out))

    stat_file = args.outcounts
    logger.debug('Writing block counts to file {}'.format(stat_file))
    with open(stat_file, 'w') as table:
        _ = table.write('{}_vs_{}\t'.format(ref_species, trg_species))
        _ = table.write('\t'.join(trg_chroms) + '\n')
        for r in ref_chroms:
            _ = table.write(r)
            for t in trg_chroms:
                _ = table.write('\t' + str(merged_counts[(r, t)]))
            _ = table.write('\n')

    logger.debug('Merging split files')

    seen_ids = set()
    with gzip.open(args.outsplits, 'wt') as dump:
        for _, split_file in sorted(created_split_files):
            with gzip.open(split_file, 'rt') as infile:
                for line in infile:
                    split_id = line.split()[4]
                    if split_id in seen_ids:
                        raise ValueError('Duplicate ID: {} ({})'.format(split_id, split_file))
                    seen_ids.add(split_id)
                    _ = dump.write(line)
    logger.debug('Done merging split files')
    if not args.keepsplits:
        for _, split_file in created_split_files:
            try:
                os.unlink(split_file)
            except (OSError, IOError) as ose:
                logger.warning('Could not remove split file: {} - {}'.format(split_file, str(ose)))
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
