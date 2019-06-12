#!/usr/bin/env python3
# coding=utf-8

import os
import sys
import traceback as trb
import argparse as argp
import json as json
import logging
import logging.config as logconf
import math
import time

import requests


logger = logging.getLogger()


class EnsemblRestService:

    def __init__(self, scale_rps=1):
        self._querying_period = 3600  # sec => 1 hour
        self._total_queries = 55000  # per period => ~15 / sec
        self._requests_per_second = int(math.floor(self._total_queries / self._querying_period))
        if scale_rps < 1:
            self._requests_per_second = int(math.floor(self._requests_per_second * scale_rps))
        self._service_url = 'http://rest.ensembl.org'
        self._request_header = {'content-type': 'application/json'}

    def _check_response(self, response):
        """
        :param response:
        :return:
        """
        if response.status_code != 200:
            if response.status_code == 429:
                # this indicates rate limiting
                # header should include waiting time
                # in seconds
                check_result = int(response.headers['Retry-After'])
            else:
                raise ConnectionError('Response header indicates error: {}'.format(response.headers))
        else:
            ratelimit_per_period = int(response.headers['X-RateLimit-Limit'])
            assert ratelimit_per_period == self._total_queries, \
                'Total rate limit changed: {}'.format(response.headers)
            ratelimit_period = int(response.headers['X-RateLimit-Period'])
            assert ratelimit_period == self._querying_period, 'Querying period changed: {}'.format(response.headers)
            check_result = 0
        return check_result

    def plain_rest_query(self, query_type):
        """
        :param query_type:
        :return:
        """
        plain_types = {'species': self._get_species_info}
        return plain_types[query_type]()

    def _get_species_info(self):
        """
        :return:
        """
        query_url = os.path.join(self._service_url, 'info', 'species')
        query_params = {'division': 'EnsemblVertebrates',
                        'hide_strain_info': 0}
        while 1:
            response = requests.get(query_url,
                                    params=query_params,
                                    headers=self._request_header)
            wait = self._check_response(response)
            if wait > 0:
                time.sleep(wait)
                continue
            break
        content = response.json()
        return content


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
    parser.add_argument('--rate-scaling', '-rs', type=float, default=1., dest='rate_scaling',
                        help='Set a scaling factor to further limit the query rates.')
    parser.add_argument('--use-logger', '-ul', type=str,
                        default='default', dest='use_logger',
                        help='Name of logger to use (default).')
    parser.add_argument('--log-config', '-lc', type=str,
                        default=log_config_default_path, dest='log_config',
                        help='Full path to JSON file containing '
                             'configuration parameters for the '
                             'loggers to use. A logger named "debug" '
                             'must be present in the configuration file.')
    parser.add_argument('--output', '-o', type=str, dest='output', required=True,
                        help='Full path to output JSON file to store '
                             'data dump. Non-existing folders will be'
                             ' created.')
    parser.add_argument('--query-type', '-qt', type=str, choices=['species'],
                        required=True, dest='query_type',
                        help='Specify the type of the query to be sent to the '
                             'Ensembl REST service.')

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


def prepare_output_file_path(output_path):
    """
    :param output_path:
    :return:
    """
    if os.path.isfile(output_path):
        logger.warning('Output file already exists: {}'.format(output_path))
    else:
        output_folders = os.path.abspath(os.path.dirname(output_path))
        os.makedirs(output_folders, exist_ok=True)
        logger.debug('Created output folder structure: {}'.format(output_folders))
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    init_logger(args)
    rest_service = EnsemblRestService()
    logger.debug('Connection to Ensembl REST service API established')
    logger.debug('Querying service for data type: {}'.format(args.query_type))
    if args.query_type in ['species']:
        # queries that do not require additional parameters
        response_content = rest_service.plain_rest_query(args.query_type)
    else:
        raise NotImplementedError
    prepare_output_file_path(args.output)
    with open(args.output, 'w') as dump:
        _ = json.dump(response_content, dump, check_circular=True,
                      sort_keys=True, ensure_ascii=True, indent=2)
    logger.debug('JSON output dumped to file')
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
