# coding=utf-8

import os as os
import csv as csv
import fnmatch as fnm
import datetime as dt


def collect_full_paths(folder, pattern, topdown=True, allow_none=False):
    """
    :param folder:
    :param pattern:
    :param topdown:
    :return:
    """
    if not topdown:
        found = os.listdir(folder)
        found = [os.path.join(folder, fp) for fp in found]
    else:
        found = []
        for root, dirs, files in os.walk(folder):
            if files:
                for f in files:
                    found.append(os.path.join(root, f))
    found = fnm.filter(found, pattern)
    if not allow_none:
        assert found, 'No files found in folder {} with pattern {}'.format(folder, pattern)
    return found


def touch_checkfile(inputfiles, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    timestr = dt.datetime.now().strftime('%A_%Y-%m-%d_%H:%M:%S')
    with open(outputfile, 'w') as outf:
        _ = outf.write(timestr + '\n')
    return outputfile


def dbg_param_list(pl):
    """
    :param pl:
    :return:
    """
    if not isinstance(pl, (list, tuple)):
        raise ValueError('Passed parameter is not list or tuple: {}'.format(type(pl)))
    print('======== First')
    print(pl[0])
    print('======== Last')
    print(pl[-1])
    print('==============')
    raise RuntimeError('Explicit stop (debug)')


def prep_metadata(mdfile, key):
    """
    :param mdfile:
    :return:
    """
    metadata = dict()
    with open(mdfile, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            metadata[row[key]] = row
    assert metadata, 'No metadata records read from file {}'.format(mdfile)
    return metadata


def prep_dset_ids(idfile, project, datatype):
    """
    :param idfile:
    :return:
    """
    dsetids = dict()
    with open(idfile, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            if row['project'] != project or row['type'] != datatype:
                continue
            if datatype == 'transcriptome':
                k = row['biosample'], row['lifestage'], row['lab'], row['experiment']
            else:
                k = row['assembly'], row['biosample'], row['lifestage'], row['lab']
            assert k not in dsetids, 'Duplicate key for IDs: {}'.format(k)
            if row['id'] == 'na':
                assert row['multid'] != 'na', 'Invalid ID for entry: {}'.format(row)
                dsetids[k] = row['multid'].split(',')
            else:
                if datatype == 'epigenome':
                    dsetids[k] = row['id']
                else:
                    dsetids[k] = row['id'], row['assembly']
    assert dsetids, 'No metadata read from ID file {}'.format(idfile)
    return dsetids
