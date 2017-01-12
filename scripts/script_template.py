#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()

    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    args = parse_command_line()

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
