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

import pandas as pd

logger = logging.getLogger()

# task specific information
NAME_MATCH_FIELDS = ['aliases', 'common_name', 'display_name', 'name']

# exclude species based on prior knowledge
IS_AQUATIC_SPECIES = ['tursiops_truncatus']
HAS_NO_CHROMOSOME_ASSEMBLY = ['callithrix_jacchus', 'cavia_porcellus',
                              'erinaceus_europaeus', 'meriones_unguiculatus',
                              'mesocricetus_auratus', 'mustela_putorius_furo',
                              'sarcophilus_harrisii', 'tupaia_belangeri',
                              'ursus_americanus', 'heterocephalus_glaber_female',
                              'heterocephalus_glaber_male']
HAS_NO_ASSEMBLY = ['suncus_murinus', 'chlorocebus_aethiops']
TARGET_WITHOUT_TRANSCRIPTOME = ['chlorocebus_sabaeus']

EXCLUDE_SPECIES = IS_AQUATIC_SPECIES + HAS_NO_CHROMOSOME_ASSEMBLY + HAS_NO_ASSEMBLY + TARGET_WITHOUT_TRANSCRIPTOME


def parse_command_line():
    """
    :return:
    """
    script_full_path = os.path.realpath(__file__)
    script_dir = os.path.dirname(script_full_path)
    log_config_default_path = os.path.join(script_dir, 'configs', 'log_config.json')

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
    parser.add_argument('--json-dump', '-jd', type=str, required=True, dest='dump',
                        help='Full path to JSON dump with Ensembl species information.')
    parser.add_argument('--species-table', '-st', type=str, required=True, dest='table',
                        help='Full path to TSV table with "species" column.')
    parser.add_argument('--output', '-o', type=str, required=True, dest='output',
                        help='Path to output TSV table with selected species information.')
    parser.add_argument('--create-info-files', '-info', action='store_true', default=False,
                        dest='info_files', help='Also create one TSV file per species.')

    args = parser.parse_args()
    return args


def init_logger(cli_args):
    """
    :param cli_args:
    :return:
    """
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


def load_species_dump(dump_path):
    """
    :param dump_path:
    :return:
    """
    with open(dump_path, 'r') as dump:
        content = json.load(dump)['species']
    logger.debug('Loaded Ensembl annotation for {} species'.format(len(content)))
    return content


def load_species_list(table_path):
    """
    :param table_path:
    :return:
    """
    table = pd.read_csv(table_path, sep='\t')
    species = table['species'].unique().tolist()
    species = [s.strip('"').lower() for s in species]
    logger.debug('Loaded {} species to be possibly selected from Ensembl annotation'.format(len(species)))
    return species


def filter_species_annotation(select_species, ensembl_species):
    """
    :param select_species:
    :param ensembl_species:
    :return:
    """

    species_annotation = col.defaultdict(list)

    for species_record in ensembl_species:
        scientific_name = species_record['name']
        if scientific_name in EXCLUDE_SPECIES:
            continue
        for name_field in NAME_MATCH_FIELDS:
            species_names = species_record[name_field]
            if not isinstance(species_names, list):
                species_names = [species_names]
            for species_id in species_names:
                any_match = [s for s in select_species if re.match(species_id, s, re.IGNORECASE) is not None]
                if any_match:
                    species_info = {'scientific_name': scientific_name,
                                    'assembly': species_record['assembly'],
                                    'accession': species_record['accession'],
                                    'taxon_id': species_record['taxon_id'],
                                    'ensembl_release': species_record['release'],
                                    'common_name': species_record['common_name'].replace(' ', '_').lower()}
                    for match in any_match:
                        species_annotation[match].append(species_info)
    logger.debug('Selected {} species from Ensembl annotation'.format(len(species_annotation)))
    unique_species = []
    for short_name, species_info in species_annotation.items():
        if not len(set(si['scientific_name'] for si in species_info)) == 1:
            raise ValueError('Non-unique species match for {}: {}'.format(short_name, species_info))
        record = species_info[0]
        record['name'] = short_name
        unique_species.append(record)
    return unique_species


def dump_selected_species(species_info, output_path, info_files):
    """
    :param species_info:
    :param output_path:
    :param info_files:
    :return:
    """
    header = ['name', 'scientific_name', 'taxon_id', 'assembly', 'accession', 'common_name', 'ensembl_release']
    table = pd.DataFrame.from_records(species_info, columns=header)
    table.sort_values(['name'], inplace=True)
    os.makedirs(os.path.abspath(os.path.dirname(output_path)), exist_ok=True)
    logger.debug('Writing species info to table {}'.format(output_path))
    table.to_csv(output_path, sep='\t', header=True, index=False)
    if info_files:
        output_folder = os.path.dirname(output_path)
        for row in table.itertuples():
            file_name = row.name + '.info'
            file_path = os.path.join(output_folder, file_name)
            with open(file_path, 'w') as info_file:
                _ = info_file.write('\t'.join([str(getattr(row, h)) for h in header]) + '\n')
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    init_logger(args)
    species_annotation = load_species_dump(args.dump)
    select_species = load_species_list(args.table)
    species_infos = filter_species_annotation(select_species, species_annotation)
    dump_selected_species(species_infos, args.output, args.info_files)
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
