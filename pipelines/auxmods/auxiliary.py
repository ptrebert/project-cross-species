# coding=utf-8

import os as os
import fnmatch as fnm
import datetime as dt


def collect_full_paths(folder, pattern, topdown=True):
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
