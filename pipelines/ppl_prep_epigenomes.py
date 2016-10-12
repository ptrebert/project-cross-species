# coding=utf-8

import os as os
import csv as csv

from pipelines.auxmods.auxiliary import collect_full_paths


def build_encode_mddict(mdfile):
    """
    :param mdfile:
    :return:
    """
    metadata = dict()
    with open(mdfile, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            metadata[row['File accession']] = row
    assert metadata, 'No ENCODE metadata records read from file {}'.format(mdfile)
    return metadata


def build_encode_eiddict(idfile):
    """
    :param idfile:
    :return:
    """
    eiddict = dict()
    with open(idfile, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            if row['Project'] != 'ENCODE' or row['type'] != 'epigenome':
                continue
            # TODO continue here


def annotate_encode_files(fpaths, mddict, eids):
    """
    :param fpaths:
    :param mddict:
    :return:
    """
    annfiles = []
    for fp in fpaths:
        fn = os.path.basename(fp).split('.')[0]
        # by construction, no missing keys
        infos = mddict[fn]
        annfiles.append((fp, fn, infos))


def link_encode_epigenomes(folder, mdfile, idfile, outfolder):
    """
    :param folder:
    :param mdfile:
    :param idfile:
    :param outfolder:
    :return:
    """
    md = build_encode_mddict(mdfile)
    bwfiles = collect_full_paths(folder, '*.bigWig', False)
