#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import json as json
import logging as logging
import logging.config as logconf
import re as re
import collections as col


import gffutils as gffutils
import pandas as pd


# This is the info as stated on
# the GENCODE website for the
# FASTA sequence files for download
GENCODE_PROTEIN_CODING_BIOTYPES = [
    'protein_coding',
    'nonsense_mediated_decay',
    'non_stop_decay',
    'IG_[A-Z]+_gene',
    'TR_[A-Z]+_gene',
    'polymorphic_pseudogene'
]

IS_PCG = re.compile('(' + '|'.join(GENCODE_PROTEIN_CODING_BIOTYPES) + ')')


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
    parser.add_argument('--gff-gene-model', '-gff', type=str, dest='genemodel',
                        required=True)
    parser.add_argument('--db-gene-model', '-db', type=str, dest='database',
                        required=True)
    parser.add_argument('--chromosomes', '-chrom', type=str, dest='chromosomes',
                        required=True)
    parser.add_argument('--gene-transcript-map', '-map', type=str, dest='mapping',
                        required=True)
    parser.add_argument('--tsv-genes', '-genes', type=str, dest='genes',
                        required=True)
    parser.add_argument('--tsv-transcripts', '-trans', type=str, dest='transcripts',
                        required=True)
    parser.add_argument('--bed-gene-promoters', '-prom', type=str, dest='promoters',
                        required=True)
    parser.add_argument('--bed-gene-bodies', '-body', type=str, dest='bodies',
                        required=True)
    parser.add_argument('--enforce-basic', '-fb', action='store_true',
                        default=False, dest='force_basic',
                        help='If no "tag" attribute is present, assume that a transcript '
                             'is not basic by default.')

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


def dump_mapping_table(transcripts, out_path):
    """
    :param transcripts:
    :param out_path:
    :return:
    """
    output_dir = os.path.abspath(os.path.dirname(out_path))
    os.makedirs(output_dir, exist_ok=True)

    transcripts.to_csv(out_path, sep='\t', columns=['name', 'gene'],
                       header=False, index=False)
    return


def dump_tabular_output(data_table, out_path, bed_style=False):
    """
    :param data_table:
    :param out_path:
    :param bed_style:
    :return:
    """
    output_dir = os.path.abspath(os.path.dirname(out_path))
    os.makedirs(output_dir, exist_ok=True)
    if bed_style:
        with open(out_path, 'w') as bed_file:
            _ = bed_file.write('#')
            data_table.to_csv(bed_file, sep='\t', header=True, index=False)
    else:
        data_table.to_csv(out_path, sep='\t', header=True, index=False)

    return


def load_gene_model_annotation(args):
    """
    :param args:
    :return:
    """
    if os.path.isfile(args.database):
        logger.debug('Loading existing database from {}'.format(args.database))
        db = gffutils.FeatureDB(args.database)
    else:
        os.makedirs(os.path.dirname(os.path.abspath(args.database)), exist_ok=True)
        db = gffutils.create_db(args.genemodel, args.database, force=True,
                                merge_strategy='create_unique', id_spec=['ID', 'Name'],
                                sort_attribute_values=False, keep_order=False,
                                checklines=50)
    return db


def load_chromosome_annotation(args):
    """
    :param args:
    :return:
    """
    chroms = dict()
    with open(args.chromosomes, 'r') as table:
        for line in table:
            if line.startswith('#'):
                continue
            name, size = line.strip().split()[:2]
            size = int(size)
            assert size > 0, 'Zero-length chromosome: {}'.format(line.strip())
            chroms[name] = size
    return chroms


def assemble_gene_info(gene, gene_chrom):
    """
    :param gene:
    :param gene_chrom:
    :return:
    """
    chrom = gene_chrom
    start = int(gene.start)
    end = int(gene.end)
    name = gene.attributes['gene_id'][0]
    assert name.startswith('ENS'), \
        'Non-Ensembl name/id: {} / {}'.format(gene, gene.attributes)
    try:
        symbol = gene.attributes['Name'][0]
    except KeyError:
        symbol = 'NoGeneSymbol'
    assert gene.strand in ['+', '-'], \
        'No strand for gene: {} / {}'.format(gene, gene.attributes)
    strand = -1 if gene.strand == '-' else 1
    length = end - start
    return chrom, start, end, name, length, strand, symbol


def assemble_transcript_info(transcript, gene_parent, gene_chrom):
    """
    :param transcript:
    :param gene_parent:
    :param gene_chrom:
    :return:
    """
    chrom = transcript.seqid
    assert chrom == gene_chrom or 'chr' + chrom == gene_chrom, \
        'Location mismatch transcript - gene: {} / {}'.format(transcript, gene_chrom)
    chrom = gene_chrom
    start = int(transcript.start)
    end = int(transcript.end)
    name = transcript.attributes['transcript_id'][0]
    assert name.startswith('ENS'), \
        'Non-Ensembl name/id: {} / {}'.format(transcript, transcript.attributes)
    try:
        symbol = transcript.attributes['Name'][0]
    except KeyError:
        symbol = 'NoTranscriptSymbol'
    assert transcript.strand in ['+', '-'], \
        'No strand for transcript: {} / {}'.format(transcript, transcript.attributes)
    strand = -1 if transcript.strand == '-' else 1
    length = end - start
    parent = gene_parent
    return chrom, start, end, name, length, strand, symbol, parent


def derive_gene_promoters(genes):
    """
    :param genes:
    :return:
    """
    genes['promoter_start'] = 0
    genes['promoter_end'] = 0

    # plus strand genes
    genes.loc[genes['strand'] > 0, 'promoter_start'] = genes.loc[genes['strand'] > 0, 'start'] - 1000
    genes.loc[genes['strand'] > 0, 'promoter_end'] = genes.loc[genes['strand'] > 0, 'start'] + 500

    # minus strand genes
    # make sure that start < end always holds
    genes.loc[genes['strand'] < 0, 'promoter_end'] = genes.loc[genes['strand'] < 0, 'end'] + 1000
    genes.loc[genes['strand'] < 0, 'promoter_start'] = genes.loc[genes['strand'] < 0, 'end'] - 500

    genes['start'] = genes['promoter_start']
    genes['end'] = genes['promoter_end']
    genes.drop(['promoter_start', 'promoter_end'], inplace=True, axis=1)
    genes['length'] = genes['end'] - genes['start']

    assert (genes['length'] > 0).all(), 'Zero length: deriving promoter coordinates failed'
    assert (genes['start'] < genes['end']).all(), 'End > start: deriving promoter coordinates failed'

    return genes


def build_annotation_tables(genes, transcripts, transcript_count):
    """
    :param genes:
    :param transcripts:
    :param transcript_count:
    :return:
    """
    gene_header = ['chrom', 'start', 'end', 'name', 'length', 'strand', 'symbol']
    genes = pd.DataFrame(genes, columns=gene_header)
    logger.debug('Built gene table with dimension: {}'.format(genes.shape))
    genes.drop_duplicates(['name'], inplace=True, keep='first')
    logger.debug('Dropped duplicate genes - new dimension: {}'.format(genes.shape))

    genes_with_transcript = set(transcript_count.keys())
    genes = genes.loc[genes['name'].isin(genes_with_transcript), :].copy()
    logger.debug('Reduced set to genes with selected transcripts - new dimension: {}'.format(genes.shape))

    genes['start'] = genes['start'].astype('int32')
    genes['end'] = genes['end'].astype('int32')
    genes['length'] = genes['length'].astype('int32')
    genes['strand'] = genes['strand'].astype('int8')
    genes.sort_values(['chrom', 'start', 'end', 'name'], inplace=True, ascending=True)

    transcript_header = gene_header + ['gene']
    transcripts = pd.DataFrame(transcripts, columns=transcript_header)
    logger.debug('Built transcript table with dimension: {}'.format(transcripts.shape))
    transcripts.drop_duplicates(['name'], inplace=True, keep='first')
    logger.debug('Dropped duplicate transcripts - new dimension: {}'.format(transcripts.shape))

    transcripts['start'] = transcripts['start'].astype('int32')
    transcripts['end'] = transcripts['end'].astype('int32')
    transcripts['length'] = transcripts['length'].astype('int32')
    transcripts['strand'] = transcripts['strand'].astype('int8')
    transcripts.sort_values(['chrom', 'start', 'end', 'name'], inplace=True, ascending=True)

    assert genes.shape[0] > 0, 'No genes retained'
    assert transcripts.shape[0] > 0, 'No transcripts retained'

    return genes, transcripts


def load_gene_transcript_information(database, chromosomes, force_basic):
    """
    :param database:
    :param chromosomes:
    :param force_basic:
    :return:
    """
    # Apparently, in the Ensembl annotation, not everything
    # is just a "transcript" with varying biotype, so need
    # to check for all in order not to miss IG/TR genes
    transcript_types = [
        'mRNA',
        'transcript',
        'C_gene_segment',
        'D_gene_segment',
        'J_gene_segment',
        'V_gene_segment'
    ]
    check_types = []
    for feat_type in database.featuretypes():
        if feat_type in transcript_types:
            check_types.append(feat_type)
    assert len(check_types) > 0, 'No feature type selected'
    logger.debug('Checking for the following transcript types: {}'.format(sorted(check_types)))

    exclude_gene = col.Counter()
    excluded_genes = set()
    exclude_transcript = col.Counter()
    excluded_transcripts = set()
    transcript_count = col.Counter()
    selected_genes = []
    selected_transcripts = []
    for gene in database.features_of_type('gene'):
        chrom = gene.seqid
        if chrom in chromosomes:
            pass
        elif 'chr' + chrom in chromosomes:
            chrom = 'chr' + chrom
        else:
            exclude_gene['location'] += 1
            excluded_genes.add(gene.id + '@Location@' + gene.seqid)
            continue
        biotypes = gene.attributes['biotype']
        assert len(biotypes) == 1, 'Multi-biotype: {} / {}'.format(gene, gene.attributes)
        biotype = biotypes[0]
        if IS_PCG.match(biotype) is None:
            exclude_gene['non_coding'] += 1
            excluded_genes.add(gene.id + '@Biotype@' + biotype)
            continue
        gene_record = assemble_gene_info(gene, chrom)

        # discard genes that are so close to the
        # chromosome border that the promoter
        # does not fit
        chrom_boundary = chromosomes[chrom]
        if (gene_record[1] - 1500) < 0 or (gene_record[2] + 1500) > chrom_boundary:
            exclude_gene['location'] += 1
            excluded_genes.add(gene.id + '@Boundary@' + chrom)
            continue
        selected_genes.append(gene_record)

        for ftype in check_types:
            for transcript in database.children(gene, featuretype=ftype):
                biotypes = transcript.attributes['biotype']
                assert len(biotypes) == 1, 'Multi-biotype: {} / {}'.format(transcript, transcript.attributes)
                biotype = biotypes[0]
                if IS_PCG.match(biotype) is None:
                    exclude_transcript['non_coding'] += 1
                    excluded_transcripts.add(transcript.id + '@Biotype@' + biotype)
                    continue
                if 'tag' in transcript.attributes:
                    if 'basic' not in transcript.attributes['tag']:
                        exclude_transcript['non_basic'] += 1
                        excluded_transcripts.add(transcript.id + '@NonBasic@' + gene.id)
                        continue
                elif force_basic:
                    exclude_transcript['non_basic'] += 1
                    excluded_transcripts.add(transcript.id + '@NonBasic@' + gene.id)
                    continue
                else:
                    pass
                transcript_count[gene_record[3]] += 1
                transcript_record = assemble_transcript_info(transcript, gene_record[3], chrom)
                selected_transcripts.append(transcript_record)

    logger.debug('Selected genes: {}'.format(len(selected_genes)))
    logger.debug('Selected transcripts: {}'.format(len(selected_transcripts)))
    logger.debug('Excluded genes based on genomic location: {}'.format(exclude_gene['location']))
    logger.debug('Excluded genes based on biotype: {}'.format(exclude_gene['non_coding']))
    logger.debug('Excluded transcripts based on biotype: {}'.format(exclude_transcript['non_coding']))
    logger.debug('Excluded non-basic transcripts: {}'.format(exclude_transcript['non_basic']))
    logger.debug('\n===================\nList of excluded genes\n=========================\n')
    logger.debug('{}'.format('\n'.join(sorted(excluded_genes))))
    logger.debug('\n===================\nList of excluded transcripts\n========================\n')
    logger.debug('{}'.format('\n'.join(sorted(excluded_transcripts))))

    return selected_genes, selected_transcripts, transcript_count


def main():
    """
    :return:
    """
    args = parse_command_line()
    init_logger(args)
    logger.debug('Loading GFF annotation from {}'.format(args.genemodel))
    db = load_gene_model_annotation(args)
    logger.debug('Loading chromosome list from {}'.format(args.chromosomes))
    chromosomes = load_chromosome_annotation(args)
    logger.debug('Loading gene and transcript information')
    genes, transcripts, transcript_count = load_gene_transcript_information(db, chromosomes, args.force_basic)
    genes, transcripts = build_annotation_tables(genes, transcripts, transcript_count)
    logger.debug('Writing transcript to gene mapping table...')
    dump_mapping_table(transcripts, args.mapping)
    logger.debug('Writing gene and transcript tables....')
    dump_tabular_output(genes, args.genes)
    dump_tabular_output(transcripts, args.transcripts)
    logger.debug('Writing BED output...')
    dump_tabular_output(genes, args.bodies, True)
    genes = derive_gene_promoters(genes)
    dump_tabular_output(genes, args.promoters, True)
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
