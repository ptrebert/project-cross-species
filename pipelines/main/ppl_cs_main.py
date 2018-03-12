# coding=utf-8

import os as os
import sys as sys
import re as re
import json as js
import csv as csv
import collections as col
import fnmatch as fnm
import datetime as dt

from ruffus import *

from pipelines.auxmods.auxiliary import collect_full_paths, touch_checkfile,\
    dbg_param_list, prep_metadata

# all features used at some point for testing
# REGION_TYPE_FEATURES = {'uprr': ['prm', 'kmf', 'gc', 'cpg', 'oecpg', 'rep', 'msig', 'dnm', 'len'],
#                         'reg5p': ['prm', 'kmf', 'gc', 'cpg', 'oecpg', 'rep', 'msig', 'dnm', 'len'],
#                         'body': ['prm', 'kmf', 'gc', 'cpg', 'oecpg', 'rep', 'msig', 'dnm', 'len']}

REGION_TYPE_FEATURES = {'uprr': ['gc', 'cpg', 'oecpg', 'msig', 'len'],
                        'reg5p': ['gc', 'cpg', 'oecpg', 'msig', 'len'],
                        'body': ['gc', 'cpg', 'oecpg', 'msig', 'len']}


def collect_mapfiles(folder, targets, queries):
    """
    :param folder:
    :param targets:
    :param queries:
    :return:
    """
    fpaths = collect_full_paths(folder, '*.idx.h5', False)
    mapfiles = list()
    for fp in fpaths:
        target, query = os.path.basename(fp).split('.')[0].split('_to_')
        if target in targets and query in queries:
            mapfiles.append({'target': target, 'query': query, 'path': fp})
    assert mapfiles, 'No map files collected from path: {}'.format(folder)
    return mapfiles


def collect_roi_files(folder, ext='*.h5'):
    """
    :param folder:
    :return:
    """
    fpaths = collect_full_paths(folder, ext, False)
    roifiles = []
    for fp in fpaths:
        annotid, regtype, _ = os.path.basename(fp).split('.', 2)
        assembly = annotid.split('_')[1]
        roifiles.append({'assembly': assembly, 'regtype': regtype, 'path': fp})
    assert roifiles, 'No ROI files annotated from path: {}'.format(folder)
    return roifiles


def collect_asc_files(folder, ext='*.tsv', key=0):
    """
    :param folder:
    :param ext:
    :return:
    """
    fpaths = collect_full_paths(folder, ext, False)
    ascfiles = dict()
    for fp in fpaths:
        assembly = os.path.basename(fp).split('_')[key]
        ascfiles[assembly] = fp
    assert ascfiles, 'No ASC files annotated from path: {}'.format(folder)
    return ascfiles


def select_compatible_libs(partner1, partner2):
    """
    :param partner1:
    :param partner2:
    :return:
    """
    p1_libs = set([p['lib'] for p in partner1])
    p2_libs = set([p['lib'] for p in partner2])
    isect = p1_libs.intersection(p2_libs)
    dnase = 1 if 'DNase' in isect else 0
    hist = len(isect) - dnase
    p1 = list(p for p in partner1 if p['lib'] in isect)
    p2 = list(p for p in partner2 if p['lib'] in isect)
    assert p1 and p2, 'Empty partner: {} and {}'.format(p1, p2)
    return p1, p2, hist, dnase, sorted(isect)


def make_groups_compatible(groupings, epibundles):
    """
    :param groupings:
    :param epibundles:
    :return:
    """
    groups = dict()
    cmatches = col.defaultdict(set)
    dupgroups = col.defaultdict(list)
    match_types = dict()
    uniq_extids = set()
    with open(groupings, 'r', newline='') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        for r in rows:
            match_types[r['comment']] = r['type']
            foo, bar = r['comment'].split('-')
            match_types[bar + '-' + foo] = r['type']
            # self-match can happen for liver, kidney
            match_types[foo + '-' + foo] = 'pos'
            match_types[bar + '-' + bar] = 'pos'
            try:
                p1, p2, h, d, clibs = select_compatible_libs(epibundles[r['partner1']], epibundles[r['partner2']])
            except AssertionError:
                print(r)
                raise
            clibs = '-'.join(clibs)
            p1_params = ' '.join([x['param'] for x in p1])
            p2_params = ' '.join([x['param'] for x in p2])
            gid = r['gid']
            extgid = gid + str(h) + str(d)
            assert extgid not in uniq_extids, 'Duplicate extended group ID: {}'.format(extgid)
            uniq_extids.add(extgid)
            entry = {'partner1': p1, 'partner2': p2, 'params1': p1_params, 'params2': p2_params, 'extgid': extgid}
            groups[gid] = entry
            c1, c2 = r['comment'].split('-')
            dupgroups[c1 + '-' + clibs].append([gid, extgid])
            dupgroups[c2 + '-' + clibs].append([gid, extgid])
            cmatches[c1].add(c2)
            cmatches[c2].add(c1)
    # workaround for monDom5 / whole blood sample
    # take as compatible to T-naive CD4
    cmatches['ncd4'].add('blood')
    match_types['ncd4-blood'] = 'pos'
    match_types['blood-ncd4'] = 'pos'
    assert groups, 'No group dictionary constructed'
    return groups, cmatches, dupgroups, match_types
            

def bundle_epigenomes(folder, mapped=False):
    """
    :param folder:
    :return:
    """
    fpaths = collect_full_paths(folder, '*/E*.h5', allow_none=True)
    collector = col.defaultdict(list)
    for fp in fpaths:
        try:
            if mapped:
                sample, foo, target, bar, ext = os.path.basename(fp).split('.')
                assert foo == 'from', 'Unexpected file name: {}'.format(os.path.basename(fp))
                eid, query, cell, lib = sample.split('_')
                key = eid, target, query
                param = '{}::{}'.format(lib, fp)
                infos = {'EID': eid, 'target': target, 'query': query, 'lib': lib,
                         'cell': cell, 'param': param, 'path': fp}
            else:
                eid, assm, cell, lib = os.path.basename(fp).split('.')[0].split('_')
                key = eid
                param = '{}::{}'.format(lib, fp)
                infos = {'EID': eid, 'assembly': assm, 'lib': lib,
                         'cell': cell, 'param': param, 'path': fp}

        except ValueError as ve:
            raise ValueError('Cannot process file {} - {}'.format(fp, str(ve)))
        collector[key].append(infos)
    return collector


def annotate_training_groups(mapfiles, roifiles, ascfiles, epidir, groupfile, outraw, cmd, jobcall):
    """
    :param mapfiles:
    :param roifiles:
    :param epidir:
    :param groupfile:
    :param outraw:
    :param cmd:
    :param jobcall:
    :return:
    """
    epibundles = bundle_epigenomes(epidir)
    groupings, matchings, dups, mtypes = make_groups_compatible(groupfile, epibundles)
    arglist = []
    uniq = set()
    for gid, partners in groupings.items():
        assm1, assm2 = partners['partner1'][0]['assembly'], partners['partner2'][0]['assembly']
        cell1, cell2 = partners['partner1'][0]['cell'], partners['partner2'][0]['cell']
        assert assm1 != assm2, \
            'Assembly self-match for training data: {} and {}'.format(partners['partner1'], partners['partner2'])
        compat_roi = [roi for roi in roifiles if roi['assembly'] == assm1]
        assm1_maps = [mapf for mapf in mapfiles if mapf['target'] == assm1]
        assert assm1_maps, 'No map files selected for target {}'.format(assm1)
        for mapf in assm1_maps:
            target, query = mapf['target'], mapf['query']
            if query not in ['mm9', 'hg19']:
                # avoid using cell lines for non-human/murine datasets
                # would be pretty pointless
                if cell1 not in ['liver', 'kidney', 'hepa', 'heart', 'ncd4'] and \
                                cell2 not in ['liver', 'kidney', 'hepa', 'heart', 'ncd4']:
                    continue
            for roi in compat_roi:
                outfolder = outraw.format(**{'target': target, 'query': query})
                outfile = '_'.join([partners['extgid'], partners['partner1'][0]['EID'],
                                    assm1, partners['partner1'][0]['cell']])
                outfile += '.to.{}.{}.h5'.format(query, roi['regtype'])
                outpath = os.path.join(outfolder, outfile)
                assert outpath not in uniq, 'Created duplicate: {}'.format(outpath)
                uniq.add(outpath)
                calc_feat = ' '.join(REGION_TYPE_FEATURES[roi['regtype']])
                if roi['regtype'] != 'reg5p':
                    params = {'target': target, 'query': query, 'genome': assm1,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params1'], 'info': partners['extgid'],
                              'ascregfile': '', 'ascregfeat': '',
                              'features': calc_feat}
                else:
                    params = {'target': target, 'query': query, 'genome': assm1,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params1'], 'info': partners['extgid'],
                              'ascregfile': '--asc-regions enh:nogroup:' + ascfiles[target],
                              'ascregfeat': 'asc', 'features': calc_feat}
                tmp = cmd.format(**params)
                arglist.append([roi['path'], outpath, tmp, jobcall])
        # same for other group partner
        compat_roi = [roi for roi in roifiles if roi['assembly'] == assm2]
        assm2_maps = [mapf for mapf in mapfiles if mapf['target'] == assm2]
        for mapf in assm2_maps:
            target, query = mapf['target'], mapf['query']
            if query not in ['mm9', 'hg19']:
                # same as above
                if cell1 not in ['liver', 'kidney', 'hepa', 'heart', 'ncd4'] and \
                                cell2 not in ['liver', 'kidney', 'hepa', 'heart', 'ncd4']:
                    continue
            for roi in compat_roi:
                outfolder = outraw.format(**{'target': target, 'query': query})
                outfile = '_'.join([partners['extgid'], partners['partner2'][0]['EID'],
                                    assm2, partners['partner2'][0]['cell']])
                outfile += '.to.{}.{}.h5'.format(query, roi['regtype'])
                outpath = os.path.join(outfolder, outfile)
                assert outpath not in uniq, 'Created duplicate: {}'.format(outpath)
                uniq.add(outpath)
                calc_feat = ' '.join(REGION_TYPE_FEATURES[roi['regtype']])
                if roi['regtype'] != 'reg5p':
                    params = {'target': target, 'query': query, 'genome': assm2,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params2'], 'info': partners['extgid'],
                              'ascregfile': '', 'ascregfeat': '',
                              'features': calc_feat}
                else:
                    params = {'target': target, 'query': query, 'genome': assm2,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params2'], 'info': partners['extgid'],
                              'ascregfile': '--asc-regions enh:nogroup:' + ascfiles[target],
                              'ascregfeat': 'asc', 'features': calc_feat}
                tmp = cmd.format(**params)
                arglist.append([roi['path'], outpath, tmp, jobcall])
    # store group info
    store_group_info(groupfile, groupings, dups)
    # store cell type matches
    store_cell_matches(groupfile, matchings, mtypes)
    if mapfiles:
        assert arglist, 'No calls created: annotate_training_groups'
    return arglist


def store_group_info(groupfile, grouping, dupgroups):
    """
    :param groupfile:
    :param grouping:
    :param dupgroups:
    :return:
    """
    timestr = dt.datetime.now().strftime('%A_%Y-%m-%d_%H:%M:%S')
    filepath = os.path.join(os.path.dirname(groupfile), 'groupinfo_ro.json')
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'groupinfo': {}, 'duplicates': {}}, dump, indent=1, sort_keys=True)
    rewrite = False
    with open(filepath, 'r') as dumped:
        gi = js.load(dumped)['groupinfo']
        if gi != grouping:
            rewrite = True
    if rewrite:
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr,
                     'groupinfo': grouping,
                     'duplicates': dupgroups}, dump, indent=1, sort_keys=True)
    return


def store_cell_matches(groupfile, matchings, matchtypes):
    """
    :param groupfile:
    :return:
    """
    primary_mouse = ['ESE14', 'liver', 'heart', 'kidney', 'ncd4', 'brain']
    # extend human list with blood - single whole blood sample for monDom5
    primary_human = ['H1hESC', 'hepa', 'ncd4', 'blood']
    primary_tissues = primary_human + primary_mouse
    serialize = dict()
    for key, vals in matchings.items():
        if key in primary_tissues:
            serialize[key] = sorted(set(primary_tissues))
        else:
            serialize[key] = sorted(list(vals))
    timestr = dt.datetime.now().strftime('%A_%Y-%m-%d_%H:%M:%S')
    filepath = os.path.join(os.path.dirname(groupfile), 'cellmatches_ro.json')
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'matchings': {}, 'matchtypes': {}}, dump, indent=1, sort_keys=True)
    rewrite = False
    with open(filepath, 'r') as dumped:
        full_dump = js.load(dumped)
        cmi = full_dump['matchings']
        if cmi != serialize:
            rewrite = True
    if rewrite:
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'matchings': serialize, 'matchtypes': matchtypes}, dump, indent=1, sort_keys=True)
    return


def make_sigmap_input(inputfiles, mapfiles, baseout, targets, queries, cmd, jobcall):
    """
    :param inputfiles:
    :param mapfiles:
    :param baseout:
    :param targets:
    :param queries:
    :param cmd:
    :param jobcall:
    :return:
    """
    arglist = []
    out_done = set()
    for mpf in sorted(mapfiles):
        fn = os.path.basename(mpf)
        trg, to, qry = fn.split('.')[0].split('_')
        if trg in targets and qry in queries:
            sigfiles = fnm.filter(inputfiles, '*/E*_{}_*.h5'.format(trg))
            for sgf in sigfiles:
                sgfn = os.path.basename(sgf)
                eid, assembly, cell, lib = sgfn.split('.')[0].split('_')
                assert trg == assembly, 'Filtering for target assembly failed: {} - {}'.format(fn, sgfn)
                mapfile = '_'.join([eid, qry, cell, lib])
                mapfile += '.'.join(['', 'from', trg, 'mapsig', 'h5'])
                mapout = os.path.join(baseout, '{}_from_{}'.format(qry, trg), mapfile)
                if mapout in out_done:
                    # for assemblies w/o cell type matches (pig, cow etc), duplicates
                    # can arise since the original cell type is used (of the source signal)
                    # simply ignore these
                    continue
                out_done.add(mapout)
                tmp = cmd.format(**{'mapfile': mpf, 'query': qry})
                arglist.append([sgf, mapout, tmp, jobcall])
    if inputfiles:
        assert arglist, 'No argument list created for signal mapping'
    return arglist


def merge_augment_featfiles(featdir, expdir, outraw, cmd, jobcall, mapped=False, matchings=None):
    """
    :param featdir:
    :param expdir:
    :param outraw:
    :param cmd:
    :param jobcall:
    :return:
    """
    if mapped:
        assert matchings is not None, 'Need to specify cell type matchings for mapped data'
        matchings = js.load(open(matchings, 'r'))['matchings']
    featfiles = collect_full_paths(featdir, '*/G[0-9]*.h5', True, True)
    expfiles = collect_full_paths(expdir, '*/T*genes.h5', False, False)
    collector = col.defaultdict(list)
    for fp in featfiles:
        parts = os.path.basename(fp).split('.')
        # k is the sample ID
        # i is the region type
        k, i = '.'.join(parts[:3]), parts[3]
        sub = os.path.split(os.path.dirname(fp))[1]
        collector[(k, sub)].append((fp, '{}::{}'.format(i, fp)))
    arglist = []
    uniq = set()
    for (group, sub), files in collector.items():
        assert len(files) == 3, 'Unexpected number of feature files {}: {}'.format(group, files)
        if mapped:
            # the pattern for mapped files consists of QRY-ASSM_TRG-CELL
            # for that combo, there is obv. no expression data, so fix here
            # to get QRY-ASSM_QRY-CELL and add TRG-CELL for ambig. cases like liver and kidney
            pat = group.split('.')[0].split('_', 2)[-1]
            query, cell_to_match = pat.split('_')
            matched_cell = matchings[cell_to_match]
            matched_cell.append(cell_to_match)
            match_exp = []
            for mc in set(matched_cell):
                foo = fnm.filter(expfiles, '*/T*_{}_{}_mRNA*'.format(query, mc))
                match_exp.extend(foo)
        else:
            pat = group.split('.')[0].split('_', 2)[-1]
            match_exp = fnm.filter(expfiles, '*/T*_{}_mRNA*'.format(pat))
        try:
            assert match_exp, 'No matched expression files for group {} / pattern {} / {}'.format(group, pat, files)
        except AssertionError as ae:
            if not mapped:
                continue
            else:
                raise ae
        mergefiles = ' '.join([t[1] for t in files])
        gid, eid, assm, cell = group.split('.')[0].split('_')
        tmp = os.path.split(os.path.dirname(files[0][0]))[1]
        if mapped:
            qry, foo, trg = tmp.split('_')
            assert foo == 'from', 'Unexpected folder structure: {}'.format(tmp)
        else:
            trg, foo, qry = tmp.split('_')
            assert foo == 'to', 'Unexpected folder structure: {}'.format(tmp)
        for fp in match_exp:
            fn = os.path.basename(fp)
            tid = fn.split('_')[0]
            if mapped:
                og = os.path.join(qry, 'from', trg, gid, 'features')
            else:
                og = os.path.join(trg, 'to', qry, gid, 'features')
            params = {'outgroup': og, 'mergefiles': mergefiles}
            outbase = outraw.format(**{'target': trg, 'query': qry})
            if mapped:
                outfile = '_'.join([gid, eid, tid, assm, cell]) + '.from.{}.feat.h5'.format(trg)
            else:
                outfile = '_'.join([gid, eid, tid, assm, cell]) + '.to.{}.feat.h5'.format(qry)
            outpath = os.path.join(outbase, outfile)
            assert outpath not in uniq, 'Duplicate created: {}'.format(outpath)
            uniq.add(outpath)
            tmp = cmd.format(**params)
            arglist.append([fp, outpath, tmp, jobcall])
    if featfiles:
        assert arglist, 'No merge arguments created'
    return arglist


def build_lola_parameter_set(inputfolder, outputfolder, cmd, jobcall):
    """
    :param inputfolder:
    :param outputfolder:
    :param cmd:
    :param jobcall:
    :return:
    """
    tissue_map = {'H1hESC': 'esc', 'ESE14': 'esc'}
    bed_files = [fn for fn in os.listdir(inputfolder) if fn.endswith('.bed')]
    if not bed_files:
        return []
    args = []
    for bf in bed_files:
        if any([x in bf for x in ['human', 'hsa']]):
            assm = 'hg19'
            dbs = ['core', 'custom']
        elif any([x in bf for x in ['mouse', 'mmu']]):
            assm = 'mm9'
            dbs = ['custom']
        else:
            raise ValueError('Species not recognized: {}'.format(bf))
        infile = os.path.join(inputfolder, bf)
        if inputfolder.endswith('unaln_genes'):
            out_infix = 'unaln-genes'
        else:
            tissue = bf.split('_')[1]
            tissue = tissue_map.get(tissue, tissue)
            out_infix = '{}-uniq-tp-genes'.format(tissue)
        for db in dbs:
            subfolder = '{}_promoter_{}_{}'.format(assm, db, out_infix)
            outfile = os.path.join(outputfolder, subfolder, 'allEnrichments.tsv')
            tmp = cmd.format(**{'infix': out_infix, 'database': db,
                                'assembly': assm})
            args.append([infile, outfile, tmp, jobcall])
    return args


def annotate_test_datasets(mapfiles, roifiles, ascfiles, mapepidir, expfiles, groupfile, outraw, cmd, jobcall):
    """
    :param mapfiles:
    :param roifiles:
    :param mapepidir:
    :param expfiles:
    :param groupfile:
    :param outraw:
    :param cmd:
    :param jobcall:
    :return:
    """
    mapepibundles = bundle_epigenomes(mapepidir, mapped=True)
    groupings = js.load(open(groupfile, 'r'))['groupinfo']
    arglist = []
    uniq = set()
    debug_record = dict()
    for mapf in mapfiles:
        target, query = mapf['target'], mapf['query']
        for (eid, trg, qry), bundle in mapepibundles.items():
            if trg == target and qry == query:
                cell = bundle[0]['cell']
                compat_exp = fnm.filter(expfiles, '*_{}_*'.format(query))
                if not compat_exp:
                    continue
                # there are ENCODE expression data for mm9 kidney/liver, but no
                # associated histone data - so make this unique here, otherwise
                # this would create duplicates
                exp_cells = set([os.path.basename(fp).split('_')[2] for fp in compat_exp])
                compat_roi = [roi for roi in roifiles if roi['assembly'] == query]
                for exp_cell in exp_cells:
                    for gid, ginfo in groupings.items():
                        if (ginfo['partner1'][0]['cell'] == exp_cell or ginfo['partner2'][0]['cell'] == exp_cell) and \
                                (ginfo['partner1'][0]['EID'] == eid or ginfo['partner2'][0]['EID'] == eid):
                            extgid = ginfo['extgid']
                            libs = [entry['lib'] for entry in ginfo['partner1']]
                            # the filter for actually available expression data exists to avoid creating
                            # test datasets on which the learner performance cannot be evaluated
                            for roi in compat_roi:
                                use_mapped = [item for item in bundle if item['lib'] in libs]
                                datafiles = ' '.join([x['param'] for x in use_mapped])

                                outfolder = outraw.format(**{'target': target, 'query': query})
                                outfile = '_'.join([extgid, eid, query, cell])
                                outfile += '.'.join(['', 'from', target, roi['regtype'], 'h5'])
                                outpath = os.path.join(outfolder, outfile)
                                try:
                                    assert outpath not in uniq, 'Created duplicate: {} / {}'.format(outpath, mapf)
                                except AssertionError:
                                    if extgid in ['G1130', 'G1330', 'G1930', 'G1730',
                                                  'G1830', 'G2330', 'G2430']:
                                        sys.stderr.write('\nWarning: skipping group duplicate {}\n'.format(outfile))
                                        # Problem here: ncd4 datasets are the only ones that
                                        # exist with the same name for mouse and human, so the
                                        # hepa-ncd4 group is selected twice in the above if statement
                                        # Manual fix to avoid too much rewriting
                                        continue
                                    else:
                                        raise AssertionError('Duplicate created: {}\n\n{}\n\n{}\n'.format(outpath, debug_record[outpath], exp_cells))
                                uniq.add(outpath)
                                debug_record[outpath] = (exp_cells, extgid, roi, ginfo)
                                calc_feat = ' '.join(REGION_TYPE_FEATURES[roi['regtype']])
                                if roi['regtype'] != 'reg5p' or query not in ascfiles:
                                    params = {'target': target, 'query': query, 'genome': query,
                                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                                              'datafiles': datafiles, 'info': extgid,
                                              'ascregfile': '', 'ascregfeat': '',
                                              'features': calc_feat}
                                else:
                                    params = {'target': target, 'query': query, 'genome': query,
                                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                                              'datafiles': datafiles, 'info': extgid,
                                              'ascregfile': '--asc-regions enh:nogroup:' + ascfiles[query],
                                              'ascregfeat': 'asc', 'features': calc_feat}
                                tmp = cmd.format(**params)
                                arglist.append([roi['path'], outpath, tmp, jobcall])
                        else:
                            #print('Skipped {} with {} and {}'.format(exp_cell, gid, eid))
                            pass
    if mapfiles and mapepibundles:
        assert arglist, 'No argument list for building test datasets created'
    return arglist


def params_predict_testdata(models, datasets, outroot, cmd, jobcall, addmd=False, mdpath=''):
    """
    :param models:
    :param datasets:
    :param outroot:
    :param cmd:
    :param jobcall:
    :param addmd:
    :param mdpath:
    :return:
    """
    if isinstance(models, str) and os.path.isdir(models):
        all_models = collect_full_paths(models, '*.pck', True, True)
    else:
        all_models = models
    if isinstance(datasets, str) and os.path.isdir(datasets):
        all_testdata = collect_full_paths(datasets, '*feat.h5', True, True)
    else:
        all_testdata = datasets

    arglist = []
    done = set()
    test_data_missing = False
    for model in all_models:
        fp, fn = os.path.split(model)
        subdir = os.path.split(fp)[-1]
        try:
            mtrg, to, mqry = subdir.split('_')
        except ValueError:
            # subdirs do not yet exist
            test_data_missing = True
            break
        assert to == 'to', 'Unexpected folder structure: {}'.format(fp)
        modelparts = fn.split('.')[0].split('_')
        groupid, eid, tid, massm, mcell = modelparts[:5]
        modeldesc = (fn.split('.', 1)[-1]).rsplit('.', 1)[0]
        out_prefix = groupid + '_trg-{}_qry-{}_'.format(mtrg, mqry)
        out_model = '_model-{}-{}-{}'.format(eid, tid, mcell) + '.' + modeldesc
        datasets = fnm.filter(all_testdata, '*/{}_*'.format(groupid))
        if not datasets:
            test_data_missing = True
            continue
        # assert datasets, 'No test datasets for group: {}'.format(groupid)
        for dat in datasets:
            fp2, fn2 = os.path.split(dat)
            assert fn2.startswith(groupid), 'Group mismatch: {} vs {}'.format(groupid, fn2)
            dqry, fr, dtrg = (os.path.split(fp2)[-1]).split('_')
            assert fr == 'from', 'Unexpected folder structure: {}'.format(fp2)
            # pairing all models and data files seems not to provide additional/useful
            # information, so restrict to matching pairs
            if dqry != mqry or dtrg != mtrg:
                continue
            dataparts = fn2.split('.')[0].split('_')
            out_data = 'data-{}-{}-{}-{}-from-{}'.format(dataparts[1], dataparts[2], dataparts[4], dqry, dtrg)
            outname = out_prefix + out_data + out_model + '.json'
            if addmd:
                mdname = outname.replace('rfreg', 'rfcls')
                mdfull = os.path.join(mdpath, groupid, mdname)
            outbase = outroot.format(**{'groupid': groupid})
            outpath = os.path.join(outbase, outname)
            assert outpath not in done, 'Duplicate created: {}'.format(outname)
            done.add(outpath)
            if addmd:
                tmp = cmd.format(**{'mdfile': mdfull})
                arglist.append([[dat, model], outpath, tmp, jobcall])
            else:
                arglist.append([[dat, model], outpath, cmd, jobcall])
    if all_models and not test_data_missing:
        assert arglist, 'No apply parameters created'
    return arglist


def make_srcsig_pairs(inputfiles, assemblies, roifiles, outbase, cmd, jobcall):
    """
    :param inputfiles:
    :param assemblies:
    :param roifiles:
    :param outbase:
    :param cmd:
    :param jobcall:
    :return:
    """
    done = set()
    arglist = []
    for trg in assemblies:
        qry = assemblies[1] if trg == assemblies[0] else assemblies[0]
        epigenomes = fnm.filter(inputfiles, '*/E*{}*.h5'.format(trg))
        for e1 in epigenomes:
            en1 = os.path.basename(e1).split('.')[0]
            _, _, cell, _ = en1.split('_')
            for e2 in epigenomes:
                en2 = os.path.basename(e2).split('.')[0]
                _, _, cell, _ = en2.split('_')
                if (e1, e2) in done:
                    continue
                done.add((e1, e2))
                done.add((e2, e1))
                assm_roi = [(d['path'], d['regtype']) for d in roifiles if d['assembly'] == trg]
                assert len(assm_roi) == 3, 'ROI file missing'
                for roip, reg in assm_roi:
                    tmp = cmd.format(**{'roifile': roip, 'target': trg, 'query': qry})
                    outfile = '_'.join(['roicorr', en1, 'vs', en2, reg]) + '.json'
                    outpath = os.path.join(outbase, outfile)
                    arglist.append([[e1, e2], outpath, tmp, jobcall])
    if inputfiles:
        assert arglist, 'No arguments created for signal correlation'
    return arglist


def make_corr_pairs(rawfiles, mapfiles, roifiles, assemblies, fp_cellmatches, outbase, cmd, jobcall):
    """
    :param rawfiles:
    :param mapfiles:
    :param roifiles:
    :param assemblies:
    :param fp_cellmatches:
    :param cmd:
    :param jobcall:
    :return:
    """
    with open(fp_cellmatches, 'r') as cm:
        cellmatches = js.load(cm)['matchings']
    done = set()
    arglist = []
    for rf in rawfiles:
        fp, fn = os.path.split(rf)
        _, qry, cell, _ = fn.split('.')[0].split('_')
        trg = assemblies[0] if qry == assemblies[1] else assemblies[1]
        if cell in ['K562', 'GM12878', 'CH12', 'MEL']:
            continue
        compatible = cellmatches[cell]
        pair_files = fnm.filter(mapfiles, '*_{}_*.from.{}.mapsig.h5'.format(qry, trg))
        for pf in pair_files:
            fp2, fn2 = os.path.split(pf)
            _, _, cell2, _ = fn2.split('_')
            if cell2 in compatible:
                if (rf, pf) in done:
                    continue
                done.add((rf, pf))
                assm_roi = [(d['path'], d['regtype']) for d in roifiles if d['assembly'] == qry]
                assert len(assm_roi) == 3, 'ROI file missing'
                for roip, reg in assm_roi:
                    tmp = cmd.format(**{'roifile': roip, 'target': trg, 'query': qry})
                    outfile = '_'.join(['roicorr', fn.split('.')[0], 'vs', fn2.rsplit('.', 1)[0], reg]) + '.json'
                    outpath = os.path.join(outbase, outfile)
                    arglist.append([[rf, pf], outpath, tmp, jobcall])
    if rawfiles and mapfiles and roifiles:
        assert arglist, 'No arguments created for signal map correlation'
    return arglist


def build_pipeline(args, config, sci_obj, pipe):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """
    if pipe is None:
        pipe = Pipeline(name=config.get('Pipeline', 'name'))
    else:
        # remove all previous tasks from pipeline
        pipe.clear()
        # turns out clear() seems NOT to have the effect
        # of really clearing the pipeline object, do it manually...
        pipe.task_names = set()
        pipe.tasks = set()

    # folder containing input data
    dir_indata = config.get('Pipeline', 'indata')
    dir_maps = config.get('Pipeline', 'refmaps')

    # base work path for processing
    # workbase = config.get('Pipeline', 'workdir')
    workbase = os.path.join(config.get('Pipeline', 'workdir'), 'norm')

    # collect raw input data in HDF format
    inputfiles = []
    tmpf = collect_full_paths(dir_indata, '*.h5', False)
    inputfiles.extend(tmpf)

    mapfiles = []
    tmpf = collect_full_paths(dir_maps, '*.idx.h5', False)
    mapfiles.extend(tmpf)

    # targets and queries to process
    use_targets = config.get('Pipeline', 'targets').split()
    use_queries = config.get('Pipeline', 'queries').split()

    # step 0: initiate pipeline with input signal files
    init = pipe.originate(task_func=lambda x: x, output=inputfiles, name='init')

    # ==========================
    # Major task: signal correlation
    # to get a feeling for the "upper" bound for the correlation across species,
    # compute pairwise correlation within species

    # compute correlations of signal in regions of interest

    dir_task_sigcorr = os.path.join(workbase, 'task_signal_correlation')
    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'corrsigroi').replace('\n', ' ')
    dir_sub_sigcorr_roi = os.path.join(dir_task_sigcorr, 'sub_roi')
    dir_refroi_bed = os.path.join(os.path.split(config.get('Pipeline', 'refroiexp'))[0], 'roi_bed')
    corrsigroi = pipe.files(sci_obj.get_jobf('inpair_out'),
                            make_srcsig_pairs(inputfiles, use_targets,
                                              collect_roi_files(dir_refroi_bed, ext='*.bed.gz'),
                                              dir_sub_sigcorr_roi, cmd, jobcall),
                            name='corrsigroi')
    corrsigroi = corrsigroi.mkdir(dir_sub_sigcorr_roi)
    corrsigroi = corrsigroi.follows(init)

    run_task_sigcorr = pipe.merge(task_func=touch_checkfile,
                                  name='task_sigcorr',
                                  input=output_from(corrsigroi),
                                  output=os.path.join(dir_task_sigcorr, 'run_task_sigcorr.chk'))

    #
    # END: major task signal correlation
    # ============================

    # ==========================
    # Major task: signal mapping
    #
    dir_task_sigmap = os.path.join(workbase, 'task_signal_mapping')

    sci_obj.set_config_env(dict(config.items('BigMemJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # map all signal tracks
    cmd = config.get('Pipeline', 'mapsig')
    dir_sub_signal = os.path.join(dir_task_sigmap, 'mapsig')
    mapsig = pipe.files(sci_obj.get_jobf('in_out'),
                        make_sigmap_input(inputfiles, mapfiles,
                                          dir_sub_signal, use_targets, use_queries,
                                          cmd, jobcall),
                        name='mapsig').mkdir(dir_sub_signal)

    # Subtask
    # compute correlations of mapped to true signal across species in regions of interest
    cmd = config.get('Pipeline', 'corrmaproi')
    dir_sub_mapcorr = os.path.join(dir_task_sigmap, 'mapcorr', 'corr_roi')
    cellmatches_file = os.path.join(os.path.split(config.get('Annotations', 'groupfile'))[0], 'cellmatches_ro.json')
    corrmaproi = pipe.files(sci_obj.get_jobf('inpair_out'),
                            make_corr_pairs(collect_full_paths(dir_indata, '*/E*.h5', topdown=False, allow_none=False),
                                            collect_full_paths(dir_sub_signal, '*/E*.h5', topdown=True, allow_none=True),
                                            collect_roi_files(dir_refroi_bed, ext='*.bed.gz'),
                                            use_targets, cellmatches_file, dir_sub_mapcorr,
                                            cmd, jobcall),
                            name='corrmaproi')
    corrmaproi = corrmaproi.mkdir(dir_sub_mapcorr)
    corrmaproi = corrmaproi.follows(mapsig)
    corrmaproi = corrmaproi.active_if(os.path.isfile(cellmatches_file))

    task_sigmap = pipe.merge(task_func=touch_checkfile,
                             name='task_sigmap',
                             input=output_from(mapsig, corrmaproi),
                             output=os.path.join(dir_task_sigmap, 'task_sigmap.chk'))

    #
    # END: major task signal mapping
    # ==============================

    # # ==================================
    # # Major task: generate training data for predicting gene expression
    # #
    dir_task_traindata_exp = os.path.join(workbase, 'task_traindata_exp')
    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # At first, consider only training groups, i.e., matched cell types between human and mouse
    # This set will give rise to machine learning models that can applied within and across species

    dir_subtask_traindata_exp_groups = os.path.join(dir_task_traindata_exp, 'compfeat_groups', '{target}_to_{query}')
    mapfiles = collect_mapfiles(dir_maps, use_targets, use_queries)
    roifiles = collect_roi_files(config.get('Pipeline', 'refroiexp'))
    ascfiles = collect_asc_files(config.get('Pipeline', 'refascreg'), key=0)
    groupfile = config.get('Annotations', 'groupfile')
    cmd = config.get('Pipeline', 'traindataexp').replace('\n', ' ')
    traindataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                     annotate_training_groups(mapfiles, roifiles, ascfiles, dir_indata,
                                                              groupfile, dir_subtask_traindata_exp_groups,
                                                              cmd, jobcall),
                                     name='traindataexp_groups')
    traindataexp_groups = traindataexp_groups.mkdir(os.path.split(dir_subtask_traindata_exp_groups)[0])

    dir_mrg_train_datasets = os.path.join(dir_task_traindata_exp, 'train_datasets', '{target}_to_{query}')
    cmd = config.get('Pipeline', 'mrgtraindataexp').replace('\n', ' ')
    mrgtraindataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                        merge_augment_featfiles(os.path.split(dir_subtask_traindata_exp_groups)[0],
                                                                dir_indata, dir_mrg_train_datasets, cmd, jobcall),
                                        name='mrgtraindataexp_groups')
    mrgtraindataexp_groups = mrgtraindataexp_groups.mkdir(os.path.split(dir_mrg_train_datasets)[0])
    mrgtraindataexp_groups = mrgtraindataexp_groups.follows(traindataexp_groups)

    task_traindata_exp = pipe.merge(task_func=touch_checkfile,
                                    name='task_traindata_exp',
                                    input=output_from(traindataexp_groups, mrgtraindataexp_groups),
                                    output=os.path.join(dir_task_traindata_exp, 'run_task_traindata_exp.chk'))

    # END: major task generate training data
    # ==============================

    # ==================================
    # Major task: generate test data for gene expression prediction (based on mapped signal)
    #
    dir_task_testdataexp_groups = os.path.join(workbase, 'task_testdata_exp')
    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # compute features for all regions of interest separately and
    # then merge region data and add expression values
    dir_sub_cmpf_testdataexp = os.path.join(dir_task_testdataexp_groups, 'compfeat_groups', '{query}_from_{target}')
    mapfiles = collect_mapfiles(dir_maps, use_targets, use_queries)
    roifiles = collect_roi_files(config.get('Pipeline', 'refroiexp'))
    expfiles = collect_full_paths(dir_indata, '*/T*genes.h5')
    groupinfos = os.path.join(os.path.dirname(config.get('Annotations', 'groupfile')), 'groupinfo_ro.json')
    matchings = os.path.join(os.path.dirname(config.get('Annotations', 'groupfile')), 'cellmatches_ro.json')
    cmd = config.get('Pipeline', 'testdataexp')
    testdataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                    annotate_test_datasets(mapfiles, roifiles, ascfiles, dir_sub_signal, expfiles,
                                                           groupinfos, dir_sub_cmpf_testdataexp,
                                                           cmd, jobcall),
                                    name='testdataexp_groups')
    testdataexp_groups = testdataexp_groups.follows(task_sigmap)
    testdataexp_groups = testdataexp_groups.mkdir(os.path.split(dir_sub_cmpf_testdataexp)[0])

    dir_mrg_test_datasets = os.path.join(dir_task_testdataexp_groups, 'test_datasets', '{query}_from_{target}')
    cmd = config.get('Pipeline', 'mrgtestdataexp')
    params_mrgtestdataexp_groups = merge_augment_featfiles(os.path.split(dir_sub_cmpf_testdataexp)[0],
                                                           dir_indata, dir_mrg_test_datasets, cmd, jobcall,
                                                           True, matchings)
    mrgtestdataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                       params_mrgtestdataexp_groups,
                                       name='mrgtestdataexp_groups')
    mrgtestdataexp_groups = mrgtestdataexp_groups.mkdir(os.path.split(dir_mrg_test_datasets)[0])
    mrgtestdataexp_groups = mrgtestdataexp_groups.follows(testdataexp_groups)

    task_testdata_exp = pipe.merge(task_func=touch_checkfile,
                                   name='task_testdata_exp',
                                   input=output_from(testdataexp_groups, mrgtestdataexp_groups),
                                   output=os.path.join(dir_task_testdataexp_groups, 'run_task_testdata_exp.chk'))

    # END: major task generate test data
    # ==============================

    # ==================================
    # Major task: train and apply models for gene expression prediction
    #
    dir_task_trainmodel_exp = os.path.join(workbase, 'task_trainmodel_exp')
    dir_task_applymodel_exp = os.path.join(workbase, 'task_applymodel_exp')

    # Set subdirectories for train/apply for gene expression status prediction
    dir_train_expstat = os.path.join(dir_task_trainmodel_exp, 'sub_status')
    dir_apply_expstat = os.path.join(dir_task_applymodel_exp, 'sub_status')

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'trainmodel_expstat_gcf')
    trainmodel_expstat_gcf = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_expstat_gcf',
                                            input=output_from(mrgtraindataexp_groups),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                            output=os.path.join(dir_train_expstat,
                                                                'gcf',
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.gbcls.gcf.all.pck'),
                                            extras=[cmd, jobcall])
    trainmodel_expstat_gcf = trainmodel_expstat_gcf.mkdir(dir_train_expstat, 'gcf')
    trainmodel_expstat_gcf = trainmodel_expstat_gcf.active_if(config.getboolean('Pipeline', 'trainmodel_expstat_gcf_run'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'apply_expstat_gcf')
    params_apply_expstat_gcf = None
    params_apply_expstat_gcf = params_predict_testdata(os.path.join(dir_train_expstat, 'gcf'),
                                                       os.path.split(dir_mrg_test_datasets)[0],
                                                       os.path.join(dir_apply_expstat, 'gcf', '{groupid}'),
                                                       cmd, jobcall)
    apply_expstat_gcf = pipe.files(sci_obj.get_jobf('inpair_out'),
                                   params_apply_expstat_gcf,
                                   name='apply_expstat_gcf')
    apply_expstat_gcf = apply_expstat_gcf.mkdir(os.path.join(dir_apply_expstat, 'gcf'))
    apply_expstat_gcf = apply_expstat_gcf.follows(task_testdata_exp)
    apply_expstat_gcf = apply_expstat_gcf.follows(trainmodel_expstat_gcf)

    # =========================================================================
    # train and apply chromatin model (using only canonical chromatin features)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'trainmodel_expstat_can').replace('\n', ' ')
    trainmodel_expstat_can = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_expstat_can',
                                            input=output_from(mrgtraindataexp_groups),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                            output=os.path.join(dir_train_expstat,
                                                                'can',
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.gbcls.can.all.pck'),
                                            extras=[cmd, jobcall])
    trainmodel_expstat_can = trainmodel_expstat_can.mkdir(os.path.join(dir_train_expstat, 'can'))
    trainmodel_expstat_can = trainmodel_expstat_can.active_if(config.getboolean('Pipeline', 'trainmodel_expstat_can_run'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # apply canonical chromatin model to test data
    cmd = config.get('Pipeline', 'apply_expstat_can')
    params_apply_expstat_can = None
    params_apply_expstat_can = params_predict_testdata(os.path.join(dir_train_expstat, 'can'),
                                                       os.path.split(dir_mrg_test_datasets)[0],
                                                       os.path.join(dir_apply_expstat, 'can', '{groupid}'),
                                                       cmd, jobcall)
    apply_expstat_can = pipe.files(sci_obj.get_jobf('inpair_out'),
                                   params_apply_expstat_can,
                                   name='apply_expstat_can')
    apply_expstat_can = apply_expstat_can.mkdir(os.path.join(dir_apply_expstat, 'can'))
    apply_expstat_can = apply_expstat_can.follows(task_testdata_exp)
    apply_expstat_can = apply_expstat_can.follows(trainmodel_expstat_can)

    # # =============================================================================
    # # =============================================================================

    cmd = config.get('Pipeline', 'summ_perf_status')
    dir_summary = os.path.join(workbase, 'task_summarize')
    summ_perf_status = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                                  name='summ_perf_status',
                                  input=output_from(apply_expstat_can, apply_expstat_gcf),
                                  output=os.path.join(dir_summary, 'agg_expstat_est.h5'),
                                  extras=[cmd, jobcall])
    summ_perf_status = summ_perf_status.mkdir(dir_summary)
    summ_perf_status = summ_perf_status.follows(task_testdata_exp)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'lola_enrich')
    dir_task_lola = os.path.join(workbase, 'task_lola')

    # before this step, need to run Notebook
    # /notebooks/utils/extract_uniq-tp-genes.ipynb
    params_lola_enrich_uniq = build_lola_parameter_set(os.path.join(dir_task_lola, 'uniq_tp_genes'),
                                                       dir_task_lola, cmd, jobcall)
    lola_enrich_uniq = pipe.files(sci_obj.get_jobf('in_out'),
                                  params_lola_enrich_uniq,
                                  name='lola_enrich_uniq')
    lola_enrich_uniq = lola_enrich_uniq.follows(summ_perf_status)
    lola_enrich_uniq = lola_enrich_uniq.active_if(os.path.isdir(dir_task_lola))

    # before this step, need to run Notebook
    # /notebooks/utils/extract_unaln_genes.ipynb
    # and perform the manual intersection steps
    # described in there
    params_lola_enrich_unaln = build_lola_parameter_set(os.path.join(dir_task_lola, 'unaln_genes'),
                                                        dir_task_lola, cmd, jobcall)
    lola_enrich_unaln = pipe.files(sci_obj.get_jobf('in_out'),
                                   params_lola_enrich_unaln,
                                   name='lola_enrich_unaln')
    lola_enrich_unaln = lola_enrich_unaln.follows(summ_perf_status)
    lola_enrich_unaln = lola_enrich_unaln.active_if(os.path.isdir(dir_task_lola))

    task_sig_cls_model = pipe.merge(task_func=touch_checkfile,
                                    name='task_sig_cls_model',
                                    input=output_from(trainmodel_expstat_can,
                                                      trainmodel_expstat_gcf,
                                                      apply_expstat_can,
                                                      apply_expstat_gcf,
                                                      summ_perf_status,
                                                      lola_enrich_uniq,
                                                      lola_enrich_unaln),
                                    output=os.path.join(workbase, 'task_sig_cls_model.chk'))
    return pipe
