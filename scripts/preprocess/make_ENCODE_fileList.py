#!/usr/bin/env python3
# coding=utf-8

"""
Isolated script to parse ENCODE metadata file, select relevant information
and create a file listing for batch download of data
"""

import os as os
import sys as sys
import csv as csv
import traceback as trb
import configparser as cfgp
import functools as funct
import subprocess as sp
import tempfile as tmpf
import operator as op
import collections as col


DL_FOLDER = '/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/rawdata/downloads'
LISTFILE = '/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/rawdata/listing_encode.txt'

MDFILE = '/home/pebert/work/code/mpggit/crossspecies/annotation/datasrc/encode/20160713_ENCODE_metadata.tsv'
SELECTFILE = '/home/pebert/work/code/mpggit/crossspecies/scripts/preprocess/ENCODE_md_selectors.ini'
# end of config


ALL_MD_FIELDS = ['File accession', 'File format', 'Output type', 'Experiment accession', 'Assay', 'Biosample term id', 'Biosample term name', 'Biosample type', 'Biosample life stage', 'Biosample sex', 'Biosample organism', 'Biosample treatments', 'Biosample subcellular fraction term name', 'Biosample phase', 'Biosample synchronization stage', 'Experiment target', 'Antibody accession', 'Library made from', 'Library depleted in', 'Library extraction method', 'Library lysis method', 'Library crosslinking method', 'Experiment date released', 'Project', 'RBNS protein concentration', 'Library fragmentation method', 'Library size range', 'Biosample Age', 'Biological replicate(s)', 'Technical replicate', 'Read length', 'Run type', 'Paired end', 'Paired with', 'Derived from', 'Size', 'Lab', 'md5sum', 'File download URL', 'Assembly', 'Platform']
USE_MD_FIELDS = ['File accession', 'File format', 'Output type', 'Experiment accession', 'Assay', 'Biosample term id', 'Biosample term name', 'Biosample type', 'Biosample life stage', 'Biosample organism', 'Experiment target', 'Antibody accession', 'Library made from', 'Library depleted in', 'Library size range', 'Biosample Age', 'Biological replicate(s)', 'Technical replicate', 'Read length', 'Run type', 'Paired end', 'Paired with', 'Derived from', 'Lab', 'File download URL', 'Assembly', 'Platform']


NORMALIZE = {'Biosample term name': {'CH12.LX': 'CH12', 'ES-Bruce4': 'ESB4', 'ES-E14': 'ESE14',
                                     'GM12878': 'GM12878', 'H1-hESC': 'H1hESC', 'HepG2': 'HepG2',
                                     'K562': 'K562', 'MEL cell line': 'MEL', 'liver': 'liver', 'kidney': 'kidney'},
             'Biosample life stage': {'adult': 'ad', 'embryonic': 'em', 'unknown': 'un',
                                      'default': 'un', 'child': 'ch', 'fetal': 'fe'},
             'Biosample organism': {'Homo sapiens': 'hsa', 'Mus musculus': 'mmu'},
             'Lab': {'Barbara Wold, Caltech': 'BWCALT', 'Ross Hardison, PennState': 'RHPSU',
                     'Bradley Bernstein, Broad': 'BBBRD', 'Bing Ren, UCSD': 'BRUCSD',
                     'John Stamatoyannopoulos, UW': 'JSUW'},
             'Run type': {'paired-ended': 'pe', 'single-ended': 'se'},
             'Experiment target': 'foo'}


def normalize_value(field, value):
    """
    :param field:
    :param value:
    :return:
    """
    if field == 'Experiment target':
        normval = value.split('-')[0]
    else:
        normval = NORMALIZE[field][value]
    return normval


def check_entry(selectors, rowdata):
    """
    :param selectors:
    :param rowdata:
    :return:
    """
    req_vals = dict()
    for k, v in selectors.items():
        checklist = v.replace('\n', '').split('@')
        field = k.replace('_',  ' ')
        req_vals[field] = checklist
    for k, v in req_vals.items():
        if rowdata[k] not in v:
            return None
    return rowdata


def get_md5(fp):
    """
    :param fp:
    :return:
    """
    out = sp.check_output('md5sum {}'.format(fp), shell=True, executable='/bin/bash')
    md5, fn = out.decode('ascii').split()
    return md5


def get_updated_md5(metadata):
    """
    :param metadata:
    :return:
    """
    mdtmp = tmpf.NamedTemporaryFile('w', encoding='ascii', delete=True)
    with open(mdtmp.name, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=USE_MD_FIELDS, extrasaction='ignore', delimiter='\t')
        writer.writeheader()
        writer.writerows(metadata)
        out.flush()
        tmp_md5 = get_md5(mdtmp.name)
    return tmp_md5


if __name__ == '__main__':
    try:
        existing_files = os.listdir(DL_FOLDER)
        existing_acc = set([f.split('.')[0] for f in existing_files])

        selectors = cfgp.ConfigParser(allow_no_value=False, interpolation=cfgp.ExtendedInterpolation())
        selectors.optionxform = lambda opt: opt
        selectors.read(SELECTFILE)
        dlfiles = []
        md_out = []
        for category in list(selectors.sections()):
            if category == 'Shared':
                continue
            rowcheck = funct.partial(check_entry, *(selectors[category],))
            with open(MDFILE, 'r') as md:
                rows = csv.DictReader(md, delimiter='\t')
                for row in rows:
                    res = rowcheck(row)
                    if res is not None:
                        if res['File accession'] not in existing_acc:
                            dlfiles.append(res['File download URL'])
                        norm = {}
                        for k, v in res.items():
                            if not v.strip():
                                norm[k.strip()] = 'n/a'
                            elif k in NORMALIZE:
                                try:
                                    norm[k.strip()] = normalize_value(k, v)
                                except KeyError:
                                    print(row)
                                    raise
                            else:
                                norm[k.strip()] = v.strip()
                        md_out.append(norm)
        md_out = sorted(md_out, key=lambda d: d['File accession'])
        if dlfiles:
            with open(LISTFILE, 'w') as listing:
                _ = listing.write('\n'.join(dlfiles))
            print('Listing written to file')
        else:
            try:
                os.unlink(LISTFILE)
            except (OSError, IOError):
                pass

        auto_md = os.path.join(os.path.dirname(MDFILE), 'encode_metadata_ro.tsv')
        auto_md_md5 = get_md5(auto_md)
        tmp_md5 = get_updated_md5(md_out)
        if auto_md_md5 != tmp_md5:
            with open(auto_md, 'w') as mdout:
                writer = csv.DictWriter(mdout, fieldnames=USE_MD_FIELDS, extrasaction='ignore', delimiter='\t')
                writer.writeheader()
                writer.writerows(md_out)
                print('Metadata file updated')

    except Exception as err:
        trb.print_exc()
        sys.exit(1)
    else:
        sys.exit(0)
