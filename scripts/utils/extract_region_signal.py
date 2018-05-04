#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import multiprocessing as mp

import numpy as np
import numpy.ma as msk
import pandas as pd


__fhgfs_base__ = '/TL/deep/fhgfs/projects/pebert/thesis'
__histone_marks__ = ['H3K4me3', 'H3K27ac', 'H3K36me3']


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--date', '-d', type=str, dest='datestamp', default='20180430')
    parser.add_argument('--force', '-f', action='store_true', default=False, dest='force')
    parser.add_argument('--workers', '-w', type=int, default=4, dest='workers')
    parser.add_argument('--cache-folder', '-cache', type=str, dest='cachefolder',
                        default=os.path.join(__fhgfs_base__,
                                             'projects/cross_species/processing/norm/caching/pipe_scripts'))
    parser.add_argument('--assemblies', '-assm', type=str, dest='assemblies',
                        default='/home/pebert/work/code/mpggit/refdata/annotation/assemblies.tsv')
    parser.add_argument('--chrom-folder', '-chrom', type=str, dest='chromfolder',
                        default=os.path.join(__fhgfs_base__,
                                             'refdata/chromsizes/chrom_auto'))

    parser.add_argument('--idx-folder', '-idx', type=str, dest='indexfolder',
                        default=os.path.join(__fhgfs_base__,
                                             'refdata/chainfiles/hdf_map'))
    parser.add_argument('--gene-folder', '-gene', type=str, dest='genefolder',
                        default=os.path.join(__fhgfs_base__,
                                             'refdata/genemodel/subsets/protein_coding/roi_hdf'))
    parser.add_argument('--source-folder', '-source', type=str, dest='sourcefolder',
                        default=os.path.join(__fhgfs_base__,
                                             'projects/cross_species/rawdata/conv/hdf'))
    parser.add_argument('--transfer-folder', '-transfer', type=str, dest='transferfolder',
                        default=os.path.join(__fhgfs_base__,
                                             'projects/cross_species/processing/norm/task_signal_mapping/mapsig'))

    parser.add_argument('--train-data', '-trd', type=str, dest='traindata',
                        default=os.path.join(__fhgfs_base__,
                                             'projects/cross_species/processing/norm/task_traindata_exp/train_datasets'))
    parser.add_argument('--test-data', '-ttd', type=str, dest='testdata',
                        default=os.path.join(__fhgfs_base__,
                                             'projects/cross_species/processing/norm/task_testdata_exp/test_datasets'))

    args = parser.parse_args()
    return args


def read_assembly_file(fpath):
    """
    :param fpath:
    :return:
    """
    use_assemblies = []
    use_species = []
    mapping = {'assembly': dict(), 'species': {}}
    with open(fpath, 'r') as tab:
        _ = tab.readline()
        for line in tab:
            if line.strip():
                a, s = line.strip().split()
                if s in ['frog', 'lizard']:
                    continue
                use_assemblies.append(a)
                use_species.append(s)
                mapping['assembly'][a] = s
                mapping['species'][s] = a
    return use_assemblies, use_species, mapping


def get_chromosome_names(assemblies, chromfolder):

    chroms = set()
    tables = [f for f in os.listdir(chromfolder) if f.endswith('.tsv')]
    for t in tables:
        assm, _ = t.split('_', 1)
        if assm in assemblies:
            fpath = os.path.join(chromfolder, t)
            with open(fpath, 'r') as file:
                for line in file:
                    if line.strip():
                        chroms.add(line.split()[0])
    return chroms


def build_region_mask(params):

    if params[1] == 'target':
        trg, direction, qry, chrom, indexfolder, genefolder = params
        qry = 'noassm'
        idx_filter = '{}_to_'.format(trg)
        filter_idx_keys = '/target/cons/mask/'
    else:
        qry, direction, trg, chrom, indexfolder, genefolder = params
        idx_filter = '{}_to_{}'.format(trg, qry)
        filter_idx_keys = '/query/cons/mask/'

    idx_files = [f for f in os.listdir(indexfolder) if f.startswith(idx_filter) and f.endswith('.idx.h5')]

    chrom_mask = None
    for idx in idx_files:
        idx_path = os.path.join(indexfolder, idx)
        with pd.HDFStore(idx_path, 'r') as hdf:
            load_keys = [k for k in hdf.keys() if k.startswith(filter_idx_keys)]
            for k in load_keys:
                k_chrom = k.split('/')[-1]
                if k_chrom != chrom:
                    continue
                chrom_data = hdf[k]
                if chrom_mask is None:
                    chrom_mask = chrom_data.values
                else:
                    try:
                        chrom_mask = np.logical_and(chrom_data.values, chrom_mask)
                    except ValueError:
                        raise ValueError('Cannot join: {} / {} / {} / {} / {} / {}'.format(idx, chrom, trg, qry, load_keys))
    # not all chromosomes exist for all species
    if chrom_mask is None:
        return '', None, None

    for genefile in os.listdir(genefolder):
        ext = genefile.split('.', 1)[-1]
        if ext in ['body.h5', 'reg5p.h5']:
            assm = genefile.split('_')[1]
            if direction == 'target':
                if assm != trg:
                    continue
            else:
                if assm != qry:
                    continue
            fpath = os.path.join(genefolder, genefile)
            with pd.HDFStore(fpath, 'r') as hdf:
                load_keys = [k for k in hdf.keys() if k.endswith(chrom)]
                assert len(load_keys) == 1, 'No region to mask for {} / {} / {} / {}'.format(ext, chrom, trg, qry)
                for k in load_keys:
                    data = hdf[k]
                    for row in data.itertuples():
                        # the signal from promoter/body
                        # regions will be extracted from
                        # the feature files
                        chrom_mask[row.start:row.end] = True
    chrom_bases = msk.array(np.arange(chrom_mask.size, dtype=np.int32), mask=chrom_mask)
    reg_sizes = np.array([s.stop - s.start for s in msk.clump_unmasked(chrom_bases)], dtype=np.int32)
    out_path = os.path.join(trg, direction, qry, chrom)
    return out_path, chrom_mask, reg_sizes


def prepare_region_masks(args, assemblies, chroms):
    """
    :param args:
    :return:
    """
    cache_name = '{}_cache_region_masks.h5'.format(args.datestamp)
    cache_file = os.path.join(args.cachefolder, cache_name)
    if os.path.isfile(cache_file) and not args.force:
        with pd.HDFStore(cache_file, 'r') as hdf:
            stored_paths = sorted(hdf.keys())
            stored_paths = [s for s in stored_paths if not s.endswith('/sizes')]
        return cache_file, stored_paths

    refs = ['hg19', 'mm9']
    queries = set(assemblies) - set(refs)
    params = []
    for r in refs:
        for c in chroms:
            params.append((r, 'target', None, c, args.indexfolder, args.genefolder))
        for q in queries:
            for c in chroms:
                params.append((q, 'query', r, c, args.indexfolder, args.genefolder))
    stored_paths = []
    try:
        with pd.HDFStore(cache_file, 'w', complevel=9, complib='blosc') as hdf:
            with mp.Pool(args.workers) as pool:
                resit = pool.imap_unordered(build_region_mask, params)
                for path, mask, sizes in resit:
                    if mask is None:
                        continue
                    hdf.put(path, pd.Series(mask), format='fixed')
                    stored_paths.append(path)
                    hdf.put(os.path.join(path, 'sizes'), pd.Series(sizes), format='fixed')
    except Exception as e:
        if os.path.isfile(cache_file):
            os.remove(cache_file)
            raise e
    return cache_file, sorted(stored_paths)


def extract_source_signal(parameters):
    """
    :param parameters:
    :return:
    """
    src_folder, mark, mask, mask_cache = parameters
    _, trg, _, _, chrom = mask.split('/')

    with pd.HDFStore(mask_cache, 'r') as hdf:
        chrom_mask = hdf[mask]

    sigfiles = [f for f in os.listdir(src_folder) if mark in f and trg in f]

    signal_dist = []
    for sf in sigfiles:
        eid, assm, cell, lib = sf.split('.')[0].split('_')
        if eid == 'EE04' and mark == 'H3K7ac':
            continue
        assert assm == trg, 'Species mismatch: {} / {}'.format(mask, sf)
        fpath = os.path.join(src_folder, sf)
        with pd.HDFStore(fpath, 'r') as hdf:
            load_keys = [k for k in hdf.keys() if k.endswith(chrom)]
            assert len(load_keys) == 1, 'Too many data records: {} / {}'.format(load_keys, sf)
            k = load_keys[0]

            chrom_data = msk.array(hdf[k].values, mask=chrom_mask)
            unmasked_slices = msk.clump_unmasked(chrom_data)
            sig_avg = np.array([np.mean(chrom_data[s.start:s.stop]) for s in unmasked_slices], dtype=np.float32)
            signal_dist.append(sig_avg)

    signal_dist = np.concatenate(signal_dist)

    out_path = os.path.join(trg, 'other', mark, chrom)
    return out_path, signal_dist


def extract_transfer_signal(params):
    """
    :param params:
    :return:
    """
    trf_folder, mark, mask, mask_cache = params
    _, trg, _, qry, chrom = mask.split('/')
    subfolder = os.path.join(trf_folder, '{}_from_{}'.format(qry, trg))
    assert os.path.isdir(subfolder), 'Expected folder {} does not exist'.format(subfolder)

    with pd.HDFStore(mask_cache, 'r') as hdf:
        chrom_mask = hdf[mask]

    sigfiles = [f for f in os.listdir(subfolder) if mark in f]

    signal_dist = []
    for sf in sigfiles:
        eid, file_qry, cell, lib = sf.split('.')[0].split('_')
        assert file_qry == qry and trg in sf, 'File mismatch: {} / {}'.format(sf, mask)
        if eid == 'EE04' and mark == 'H3K27ac':
            continue
        fpath = os.path.join(subfolder, sf)
        with pd.HDFStore(fpath, 'r') as hdf:
            load_keys = [k for k in hdf.keys() if k.endswith(chrom)]
            assert len(load_keys) == 1, 'Too many data records: {} / {}'.format(load_keys, sf)
            k = load_keys[0]

            chrom_data = msk.array(hdf[k].values, mask=chrom_mask)
            unmasked_slices = msk.clump_unmasked(chrom_data)
            sig_avg = np.array([np.mean(chrom_data[s.start:s.stop]) for s in unmasked_slices], dtype=np.float32)
            signal_dist.append(sig_avg)

    signal_dist = np.concatenate(signal_dist)
    out_path = os.path.join(trg, qry, 'other', mark, chrom)
    return out_path, signal_dist


def collect_source_signal(args, mask_cache_file, cached_masks):
    """
    :return:
    """
    cache_name = '{}_cache_source_signal.h5'.format(args.datestamp)
    cache_file = os.path.join(args.cachefolder, cache_name)
    if os.path.isfile(cache_file) and not args.force:
        with pd.HDFStore(cache_file, 'r') as hdf:
            stored_paths = sorted([k for k in hdf.keys() if 'other' in k])
        return cache_file, stored_paths

    target_masks = [m for m in cached_masks if 'target' in m]
    params = []
    for t in target_masks:
        for m in __histone_marks__:
            params.append((args.sourcefolder, m, t, mask_cache_file))

    stored_paths = []
    try:
        with pd.HDFStore(cache_file, 'w', complevel=9, complib='blosc') as hdf:
            with mp.Pool(args.workers) as pool:
                resit = pool.imap_unordered(extract_source_signal, params)
                for path, signals in resit:
                    hdf.put(path, pd.Series(signals), format='fixed')
                    stored_paths.append(path)
    except Exception as e:
        if os.path.isfile(cache_file):
            os.remove(cache_file)
            raise e
    return cache_file, stored_paths


def collect_transfer_signal(args, mask_cache_file, cached_masks):
    """
    :param args:
    :param mask_cache_file:
    :param cached_masks:
    :return:
    """
    cache_name = '{}_cache_transfer_signal.h5'.format(args.datestamp)
    cache_file = os.path.join(args.cachefolder, cache_name)
    if os.path.isfile(cache_file) and not args.force:
        with pd.HDFStore(cache_file, 'r') as hdf:
            stored_paths = sorted([k for k in hdf.keys() if 'other' in k])
        return cache_file, stored_paths

    query_masks = [m for m in cached_masks if 'query' in m]
    params = []
    for q in query_masks:
        for m in __histone_marks__:
            params.append((args.transferfolder, m, q, mask_cache_file))
    stored_paths = []
    try:
        with pd.HDFStore(cache_file, 'w', complevel=9, complib='blosc') as hdf:
            with mp.Pool(args.workers) as pool:
                resit = pool.imap_unordered(extract_transfer_signal, params)
                for path, signals in resit:
                    hdf.put(path, pd.Series(signals), format='fixed')
                    stored_paths.append(path)
    except Exception as e:
        if os.path.isfile(cache_file):
            os.remove(cache_file)
            raise e
    return cache_file, stored_paths


def load_signal_averages(parameters):
    """
    :param parameters:
    :return:
    """
    cache_path, root_folder, regtype = parameters
    parts = cache_path.split('/')
    if parts[2] == 'other':
        # target centric - load train datasets
        mark, chrom = parts[3], parts[4]
        filter_folders = '{}_to_'.format(parts[1])
        subfolders = [d for d in os.listdir(root_folder) if d.startswith(filter_folders)]
    else:
        # target/query combination - load test datasets
        mark, chrom = parts[4], parts[5]
        filter_folders = '{}_from_{}'.format(parts[2], parts[1])
        subfolders = [d for d in os.listdir(root_folder) if d.startswith(filter_folders)]

    ftext = 'reg5p' if regtype == 'promoter' else 'body'
    ftname = 'ftmsig_{}_abs_mean_{}'.format(mark, ftext)

    assert len(subfolders) > 0, 'No subfolders found for {} / {}'.format(cache_path, regtype)

    seen_epigenomes = set()
    signal_dist = []
    for sub in subfolders:
        featfiles = os.listdir(os.path.join(root_folder, sub))
        for ff in featfiles:
            if ff.startswith('G99'):
                continue
            gid, eid, tid, species, _ = ff.split('_')
            if eid in seen_epigenomes:
                continue
            if eid == 'EE04' and mark == 'H3K27ac':
                seen_epigenomes.add(eid)
                continue
            fpath = os.path.join(root_folder, sub, ff)
            with pd.HDFStore(fpath, 'r') as hdf:
                load_keys = [k for k in hdf.keys() if k.endswith(chrom)]
                assert len(load_keys) == 1, 'Too many data records: {}'.format(load_keys)
                k = load_keys[0]
                chrom_data = hdf[k]
                signal_dist.append(chrom_data[ftname].values)

    signal_dist = np.concatenate(signal_dist)
    out_path = cache_path.replace('other', regtype)
    return out_path, signal_dist


def augment_signal_cache(args, cached_paths, cache_file, root_folder):
    """
    :param args:
    :param cached_paths:
    :param cache_file:
    :param root_folder:
    :return:
    """
    check_name = '{}_augment-complete.{}.chk'.format(args.datestamp, os.path.split(root_folder)[-1])
    check_file = os.path.join(os.path.dirname(cache_file), check_name)
    if os.path.isfile(check_file):
        return 0
    params = []
    for p in cached_paths:
        params.append((p, root_folder, 'body'))
        params.append((p, root_folder, 'promoter'))

    with pd.HDFStore(cache_file, 'a') as hdf:
        with mp.Pool(args.workers) as pool:
            resit = pool.imap_unordered(load_signal_averages, params)
            for path, signals in resit:
                hdf.put(path, pd.Series(signals), format='fixed')

    with open(check_file, 'w') as check:
        pass
    return 0


def main():
    """
    :return:
    """
    args = parse_command_line()
    assemblies, species, lut = read_assembly_file(args.assemblies)
    chroms = get_chromosome_names(assemblies, args.chromfolder)
    region_mask_cache, stored_masks = prepare_region_masks(args, assemblies, chroms)
    print('Region mask cache built')
    src_signal_cache, stored_source = collect_source_signal(args, region_mask_cache, stored_masks)
    print('Source signal cached')
    trf_signal_cache, stored_transfer = collect_transfer_signal(args, region_mask_cache, stored_masks)
    print('Transfer signal cached')
    _ = augment_signal_cache(args, stored_source, src_signal_cache, args.traindata)
    print('Source signal data augmented')
    _ = augment_signal_cache(args, stored_transfer, trf_signal_cache, args.testdata)
    print('Transfer signal data augmented')
    return 0


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
