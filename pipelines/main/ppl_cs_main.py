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
                if roi['regtype'] != 'reg5p':
                    params = {'target': target, 'query': query, 'genome': assm1,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params1'], 'info': partners['extgid'],
                              'ascregfile': '', 'ascregfeat': ''}
                else:
                    params = {'target': target, 'query': query, 'genome': assm1,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params1'], 'info': partners['extgid'],
                              'ascregfile': '--asc-regions enh:nogroup:' + ascfiles[target],
                              'ascregfeat': 'asc'}
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
                if roi['regtype'] != 'reg5p':
                    params = {'target': target, 'query': query, 'genome': assm2,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params2'], 'info': partners['extgid'],
                              'ascregfile': '', 'ascregfeat': ''}
                else:
                    params = {'target': target, 'query': query, 'genome': assm2,
                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                              'datafiles': partners['params2'], 'info': partners['extgid'],
                              'ascregfile': '--asc-regions enh:nogroup:' + ascfiles[target],
                              'ascregfeat': 'asc'}
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
    primary_human = ['H1hESC', 'hepa', 'ncd4']
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


def make_meta_traindata(infiles, outdir, cmd, jobcall):
    """
    :param infiles:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    sort_files = col.defaultdict(list)
    for inpf in infiles:
        base, to, query, _, _ = os.path.basename(inpf).split('.')
        assert to == 'to', 'Unexpected filename: {}'.format(inpf)
        _, eid, tid, target, cell = base.split('_')
        sort_files[(target, query)].append(('_'.join([eid, tid, cell]), inpf))
    arglist = []
    done = set()
    for (target, query), values in sort_files.items():
        call_inputs = []
        label_inputs = ''
        for label, fp in values:
            call_inputs.append(fp)
            label_inputs += '{}::{} '.format(label, fp)
        tmp = cmd.format(**{'mergefiles': label_inputs, 'target': target, 'query': query})
        outfile = '{}_to_{}_meta.feat.h5'.format(target, query)
        outpath = os.path.join(outdir, outfile)
        assert outpath not in done, 'Created duplicate: {}'.format(outpath)
        arglist.append([call_inputs, outpath, tmp, jobcall])
    if infiles:
        assert arglist, 'No call parameters create for meta train datasets'
    return arglist


def match_prior_traindata(traindata, priordata, suffix, cmd, jobcall):
    """
    :param traindata:
    :param priordata:
    :param suffix:
    :param cmd:
    :param jobcall:
    :return:
    """
    if traindata and priordata:
        #assert len(traindata) == len(priordata),\
        #    'Length mismatch: {} train vs {} priors'.format(len(traindata), len(priordata))
        if len(traindata) != len(priordata):
            sys.stderr.write('\nWarning: only {} derived feature datasets for {} training datasets. '
                             'Return empty.\n'.format(len(priordata), len(traindata)))
            return []
    arglist = []
    for tf, pf in zip(sorted(traindata), sorted(priordata)):
        tfp, tfn = os.path.split(tf)
        pfp, pfn = os.path.split(pf)
        gid = pfn.split('_')[0]
        trg, _, qry = os.path.split(tfp)[1].split('_')
        assert os.path.split(tfp)[1] == os.path.split(pfp)[1], 'Sub path does not match: {} vs {}'.format(tfp, pfp)
        assert tfn.split('.')[0] == pfn.split('.')[0], 'Dataset mismatch: {} vs {}'.format(tfn, pfn)
        outname = tfn.replace('.h5', '{}.h5'.format(suffix))
        outpath = os.path.join(tfp, outname)
        tmp = cmd.format(**{'featdata': tf, 'priordata': pf, 'outputfile': outpath,
                            'outgroup': '{}/to/{}/{}/feat'.format(trg, qry, gid)})
        arglist.append([[tf, pf], outpath, tmp, jobcall])
    if priordata and traindata:
        assert arglist, 'No argument calls created'
    return arglist


def match_prior_testdata(testdata, priordata, suffix, cmd, jobcall):
    """
    :param testdata:
    :param priordata:
    :param suffix:
    :param cmd:
    :param jobcall:
    :return:
    """
    if not (testdata and priordata):
        return []
    arglist = []
    # G0951_data-ED13-TE11-hepa-mm9-from-hg19.seq-cprob.h5
    for tf in testdata:
        tfp, tfn = os.path.split(tf)
        _, species_match = os.path.split(tfp)
        qry, fr, trg = species_match.split('_')
        assert fr == 'from', 'Unexpected folder structure: {}'.format(tfp)
        try:
            dset, _, dtrg, _ = tfn.split('.', 3)
        except ValueError:
            raise ValueError('Unexpected number of components: {}'.format(tfn))
        assert dtrg == trg, 'Target species mismatch: {} vs {}'.format(species_match, tfn)
        gid, eid, tid, dqry, cell = dset.split('_')
        assert dqry == qry, 'Query species mismatch: {} vs {}'.format(species_match, tfn)
        use_priors = fnm.filter(priordata, '*{}*{}-{}-{}-{}-from-{}.*.h5'.format(gid, eid, tid, cell, qry, trg))
        try:
            assert len(use_priors) == 1, 'No/ambig. metadata with run information selected: {}'.format(use_priors)
        except AssertionError:
            if len(use_priors) == 0:
                sys.stderr.write('\n{}: no priors for {} - {} - {} - {}'.format(suffix, gid, trg, qry, tfn))
                # this can happen during incomplete pipeline runs, i.e. some prior information
                # is available and match_prior_testdata is executed; return empty here instead of
                # raising AssertionError
                # raise AssertionError('No prior information available for group: {} - TRG {} - QRY {}'.format(gid, trg, qry))
                arglist = []
                break
            else:
                raise AssertionError('Ambig. prior information selected: {}'.format(use_priors))
        pdata_file = use_priors[0]
        outname = tfn.replace('.h5', '{}.h5'.format(suffix))
        outpath = os.path.join(tfp, outname)
        inputfiles = use_priors + [tf]
        tmp = cmd.format(**{'featdata': tf, 'outputfile': outpath,
                            'priordata': pdata_file, 'outgroup': '{}/from/{}/{}/feat'.format(qry, trg, gid)})
        arglist.append([inputfiles, outpath, tmp, jobcall])
    return arglist


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
                                    if extgid in ['G1140', 'G1340', 'G1940', 'G1730',
                                                  'G1830', 'G2340', 'G2440']:
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
                                if roi['regtype'] != 'reg5p':
                                    params = {'target': target, 'query': query, 'genome': query,
                                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                                              'datafiles': datafiles, 'info': extgid,
                                              'ascregfile': '', 'ascregfeat': ''}
                                else:
                                    params = {'target': target, 'query': query, 'genome': query,
                                              'regtype': roi['regtype'], 'mapfile': mapf['path'],
                                              'datafiles': datafiles, 'info': extgid,
                                              'ascregfile': '--asc-regions enh:nogroup:' + ascfiles[query],
                                              'ascregfeat': 'asc'}
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


def params_predict_testdata_prior(modelroot, dataroot, outroot, cmd, jobcall, addmd=False, mdpath=''):
    """
    :param modelroot:
    :param dataroot:
    :param cmd:
    :param jobcall:
    :return:
    """
    all_models = collect_full_paths(modelroot, '*psig*.pck', True, True)
    all_testdata = collect_full_paths(dataroot, '*feat.cprob.h5', True, True)

    arglist = []
    done = set()
    test_data_missing = False
    for model in all_models:
        if '.psig.' not in model:
            continue
        fp, fn = os.path.split(model)
        subdir = os.path.split(fp)[-1]
        mtrg, to, mqry = subdir.split('_')
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


def params_predict_testdata_meta(modelroot, dataroot, outroot, cmd, jobcall, addmd=False, mdpath=''):
    """
    :param modelroot:
    :param dataroot:
    :param cmd:
    :param jobcall:
    :return:
    """
    all_models = collect_full_paths(modelroot, '*.pck', True, True)
    all_testdata = collect_full_paths(dataroot, '*feat.h5', True, True)

    arglist = []
    done = set()
    test_data_missing = False
    for model in all_models:
        fp, fn = os.path.split(model)
        subdir = os.path.split(fp)[-1]
        assert subdir == 'meta', 'Unexpected path for meta models: {}'.format(model)
        mtrg, to, mqry = fn.split('_')[:3]
        assert to == 'to', 'Unexpected filename: {}'.format(fn)
        modeldesc = (fn.split('.', 1)[-1]).rsplit('.', 1)[0]
        out_model = '_model-meta.to.{}'.format(mqry) + '.' + modeldesc
        datasets = fnm.filter(all_testdata, '*/{}_from_{}/G*'.format(mqry, mtrg))
        if not datasets:
            test_data_missing = True
            continue
        # assert datasets, 'No test datasets for group: {}'.format(groupid)
        for dat in datasets:
            fp2, fn2 = os.path.split(dat)
            assert '.from.{}'.format(mtrg) in fn2 and mqry in fn2, 'Target/query mismatch: {} vs {}'.format(fn, fn2)
            dqry, fr, dtrg = (os.path.split(fp2)[-1]).split('_')
            assert fr == 'from', 'Unexpected folder structure: {}'.format(fp2)
            dataparts = fn2.split('.')[0].split('_')
            groupid = dataparts[0]
            out_data = 'data-{}-{}-{}-{}-from-{}'.format(dataparts[1], dataparts[2], dataparts[4], dqry, dtrg)
            out_prefix = '{}_trg-{}_qry-{}_'.format(groupid, mtrg, mqry)
            outname = out_prefix + out_data + out_model + '.json'
            if addmd:
                mdname = outname.replace('rfreg', 'rfcls')
                mdfull = os.path.join(mdpath, mdname)
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


def annotate_modelfile(fpath, datasets, groupings):
    """
    :param fpath:
    :return:
    """
    path, name = os.path.split(fpath)
    sample_desc, model_desc = name.split('.', 1)
    groupid, model_epi, model_trans, target, model_cell = sample_desc.split('_')
    group, hist, dnase = groupid[:3], groupid[3], groupid[4]
    _, query, model_type = model_desc.split('.', 2)
    model_type = model_type.rsplit('.', 1)[0]
    test_type = groupings[group]['type']
    epi_cell = datasets[model_epi]['biosample']
    trans_cell = datasets[model_trans]['biosample']
    epi_proj = datasets[model_epi]['project']
    trans_proj = datasets[model_trans]['project']
    # since model == training dataset, the cell types
    # have to match
    assert epi_cell == trans_cell, 'Biosample mismatch for model file: {}'.format(fpath)
    md = {'path': fpath, 'group': group, 'num_hist': hist, 'num_dnase': dnase,
          'group_test_type': test_type, 'model_spec': model_type, 'model_target': target,
          'model_query': query, 'model_cell': epi_cell, 'model_epigenome': model_epi,
          'model_transcriptome': model_trans, 'model_epi_project': epi_proj,
          'model_trans_project': trans_proj}
    return groupid, md


def annotate_test_dataset(fpath, datasets, groupings, cmatchtypes):
    """
    :param fpath:
    :param datasets:
    :param groupings:
    :return:
    """
    path, name = os.path.split(fpath)
    sample_desc, _, target, _, _ = name.split('.')
    groupid, data_epi, data_trans, query, data_cell = sample_desc.split('_')
    group, hist, dnase = groupid[:3], groupid[3], groupid[4]
    epi_cell = datasets[data_epi]['biosample']
    trans_cell = datasets[data_trans]['biosample']
    epi_proj = datasets[data_epi]['project']
    trans_proj = datasets[data_trans]['project']
    try:
        test_type = cmatchtypes[epi_cell + '-' + trans_cell]
    except KeyError:
        if epi_cell == trans_cell:
            test_type = 'pos'
        else:
            raise
    group_test = groupings[group]['type']
    epi_assm = datasets[data_epi]['assembly']
    trans_assm = datasets[data_trans]['assembly']
    assert trans_assm == query, 'Assembly mismatch for transcriptome / query in test_dataset: {}'.format(fpath)
    assert epi_assm == target, 'Assembly mismatch for epigenome / target in test dataset: {}'.format(fpath)
    assert epi_assm != trans_assm, 'Assembly match in test dataset: {}'.format(fpath)
    md = {'path': fpath, 'group': group, 'num_hist': hist, 'num_dnase': dnase,
          'group_test_type': group_test, 'data_test_type': test_type,
          'data_epigenome': data_epi, 'data_transcriptome': data_trans,
          'data_target': target, 'data_query': query, 'data_epi_cell': epi_cell,
          'data_trans_cell': trans_cell, 'data_epi_project': epi_proj,
          'data_trans_project': trans_proj}
    return groupid, md


def annotate_test_output(paramset, datasets, groupings, cellmatches, task, outfile):
    """
    :param paramset:
    :param datasets:
    :param groupings:
    :param cellmatches:
    :param outfile:
    :return:
    """
    if not os.path.isdir(os.path.dirname(outfile)):
        # pipeline potentially restarted, folder does not yet exist
        return 1
    cmatches = cellmatches['matchtypes']
    collect_runs = col.defaultdict(list)
    for p in paramset:
        assert isinstance(p, list), 'Unexpected parameter set: {}'.format(p)
        assert isinstance(p[0], list) and len(p[0]) == 2, 'Expected pair of files: {}'.format(p[0])
        datafile = p[0][0]
        modelfile = p[0][1]
        gid1, data_md = annotate_test_dataset(datafile, datasets, groupings, cmatches)
        gid2, model_md = annotate_modelfile(modelfile, datasets, groupings)
        assert gid1 == gid2, 'Group ID mismatch for files {} / {}'.format(datafile, modelfile)
        mtype, rtype = categorize_test_run(data_md, model_md, cmatches)
        collect_runs[gid1].append({'run_test_type': rtype, 'run_spec_match': mtype, 'setting': task,
                                   'model_metadata': model_md, 'data_metadata': data_md,
                                   'run_file': p[1]})
    with open(outfile, 'w') as outf:
        js.dump(collect_runs, outf, indent=1, sort_keys=True)
    return 0


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
            if cell in ['K562', 'MEL', 'CH12', 'GM12878']:
                continue
            for e2 in epigenomes:
                en2 = os.path.basename(e2).split('.')[0]
                _, _, cell, _ = en2.split('_')
                if cell in ['K562', 'MEL', 'CH12', 'GM12878']:
                    continue
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


def categorize_test_run(datamd, modelmd, cmatches):
    """
    :param datamd:
    :param modelmd:
    :param cmatches:
    :return:
    """
    try:
        model_dataepi = cmatches[modelmd['model_cell'] + '-' + datamd['data_epi_cell']]
    except KeyError:
        model_dataepi = 'neg'
    try:
        model_datatrans = cmatches[modelmd['model_cell'] + '-' + datamd['data_trans_cell']]
    except KeyError:
        model_datatrans = 'neg'

    if (model_dataepi, model_datatrans) == ('pos', 'pos'):
        run_test_type = 'pos'
    elif (model_dataepi, model_datatrans) == ('pos', 'neg'):
        run_test_type = 'expneg'
    elif (model_dataepi, model_datatrans) == ('neg', 'pos'):
        run_test_type = 'epineg'
    elif (model_dataepi, model_datatrans) == ('neg', 'neg'):
        run_test_type = 'rand'
    else:
        run_test_type = 'undefined'
    if datamd['data_target'] == modelmd['model_target'] and \
        datamd['data_query'] == modelmd['model_query']:
        run_spec_match = 'cons'
    elif datamd['data_query'] == modelmd['model_target']:
        run_spec_match = 'self'
    elif datamd['data_target'] != modelmd['model_target'] and \
        datamd['data_query'] == modelmd['model_query']:
        run_spec_match = 'trgmix'
    elif datamd['data_target'] == modelmd['model_target'] and \
        datamd['data_query'] != modelmd['model_query']:
        run_spec_match = 'qrymix'
    elif datamd['data_target'] != modelmd['model_target'] and \
        datamd['data_query'] != modelmd['model_query']:
        run_spec_match = 'dblmix'
    else:
        run_spec_match = 'undefined'
    return run_spec_match, run_test_type


def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

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

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # map all signal tracks
    cmd = config.get('Pipeline', 'mapsig').replace('\n', ' ')
    dir_sub_signal = os.path.join(dir_task_sigmap, 'mapsig')
    mapsig = pipe.files(sci_obj.get_jobf('in_out'),
                        make_sigmap_input(inputfiles, mapfiles,
                                          dir_sub_signal, use_targets, use_queries,
                                          cmd, jobcall),
                        name='mapsig').mkdir(dir_sub_signal)

    # Subtask
    # compute correlations of mapped to true signal across species in regions of interest
    cmd = config.get('Pipeline', 'corrmaproi').replace('\n', ' ')
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

    #
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
    cmd = config.get('Pipeline', 'testdataexp').replace('\n', ' ')
    testdataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                    annotate_test_datasets(mapfiles, roifiles, ascfiles, dir_sub_signal, expfiles,
                                                           groupinfos, dir_sub_cmpf_testdataexp,
                                                           cmd, jobcall),
                                    name='testdataexp_groups')
    testdataexp_groups = testdataexp_groups.follows(task_sigmap)
    testdataexp_groups = testdataexp_groups.mkdir(os.path.split(dir_sub_cmpf_testdataexp)[0])

    dir_mrg_test_datasets = os.path.join(dir_task_testdataexp_groups, 'test_datasets', '{query}_from_{target}')
    cmd = config.get('Pipeline', 'mrgtestdataexp').replace('\n', ' ')
    mrgtestdataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                       merge_augment_featfiles(os.path.split(dir_sub_cmpf_testdataexp)[0],
                                                               dir_indata, dir_mrg_test_datasets, cmd, jobcall, True, matchings),
                                       name='mrgtestdataexp_groups')\
                                        .mkdir(os.path.split(dir_mrg_test_datasets)[0])\
                                        .follows(testdataexp_groups)

    task_testdata_exp = pipe.merge(task_func=touch_checkfile,
                                   name='task_testdata_exp',
                                   input=output_from(testdataexp_groups, mrgtestdataexp_groups),
                                   output=os.path.join(dir_task_testdataexp_groups, 'run_task_testdata_exp.chk'))

    #
    # END: major task generate test data
    # ==============================

    # ==================================
    # Major task: train and apply models for gene expression prediction
    #
    dir_task_trainmodel_exp = os.path.join(workbase, 'task_trainmodel_exp')
    dir_task_applymodel_exp = os.path.join(workbase, 'task_applymodel_exp')

    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # train sequence model to classify genes as active/inactive (TPM >= 1)
    dir_train_expstat = os.path.join(dir_task_trainmodel_exp, 'sub_status')

    cmd = config.get('Pipeline', 'trainmodel_expstat_seq').replace('\n', ' ')
    trainmodel_expstat_seq = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_expstat_seq',
                                            input=output_from(mrgtraindataexp_groups),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                            output=os.path.join(dir_train_expstat,
                                                                'seq',
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.rfcls.seq.all.pck'),
                                            extras=[cmd, jobcall])
    trainmodel_expstat_seq = trainmodel_expstat_seq.mkdir(dir_train_expstat, 'seq')

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_seq_out').replace('\n', ' ')
    extrain_seq_out = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                     name='extrain_seq_out',
                                     input=output_from(trainmodel_expstat_seq),
                                     filter=formatter('(?P<SAMPLE>[\w\.]+).rfcls.seq.all.pck'),
                                     output='{path[0]}/{SAMPLE[0]}.seq-cprob.h5',
                                     extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_seq_out').replace('\n', ' ')
    params_addtrain_seq_out = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                       pattern='*.feat.h5', allow_none=True),
                                                    collect_full_paths(os.path.join(dir_train_expstat, 'seq'),
                                                                       pattern='*seq-cprob.h5',
                                                                       allow_none=True),
                                                    '.seqout', cmd, jobcall)

    addtrain_seq_out = pipe.files(sci_obj.get_jobf('ins_out'),
                                  params_addtrain_seq_out,
                                  name='addtrain_seq_out')
    addtrain_seq_out = addtrain_seq_out.follows(extrain_seq_out)
    addtrain_seq_out = addtrain_seq_out.follows(task_traindata_exp)
    params_addtrain_seq_out = list()

    # apply sequence model to test data
    dir_apply_expstat = os.path.join(dir_task_applymodel_exp, 'sub_status')

    cmd = config.get('Pipeline', 'apply_expstat_seq').replace('\n', ' ')
    params_apply_expstat_seq = params_predict_testdata(os.path.join(dir_train_expstat, 'seq'),
                                                       os.path.split(dir_mrg_test_datasets)[0],
                                                       os.path.join(dir_apply_expstat, 'seq', '{groupid}'),
                                                       cmd, jobcall)
    apply_expstat_seq = pipe.files(sci_obj.get_jobf('inpair_out'),
                                   params_apply_expstat_seq,
                                   name='apply_expstat_seq')
    apply_expstat_seq = apply_expstat_seq.mkdir(os.path.join(dir_apply_expstat, 'seq'))
    apply_expstat_seq = apply_expstat_seq.follows(mrgtestdataexp_groups)
    apply_expstat_seq = apply_expstat_seq.follows(trainmodel_expstat_seq)
    params_apply_expstat_seq = list()

    # extract sequence model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_seq_out').replace('\n', ' ')
    extest_seq_out = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                  name='extest_seq_out',
                                  input=output_from(apply_expstat_seq),
                                  filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+seq\.all\.json'),
                                  output='{path[0]}/{GID[0]}_{DSET[0]}.seq-cprob.h5',
                                  extras=[cmd, jobcall])

    # add sequence model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_seq_out').replace('\n', ' ')
    params_addtest_seq_out = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                     '*.feat.h5', allow_none=True),
                                                  collect_full_paths(os.path.join(dir_apply_expstat, 'seq'),
                                                                     '*.seq-cprob.h5', allow_none=True),
                                                  '.seqout', cmd, jobcall)
    addtest_seq_out = pipe.files(sci_obj.get_jobf('ins_out'),
                                 params_addtest_seq_out,
                                 name='addtest_seq_out')
    addtest_seq_out = addtest_seq_out.follows(apply_expstat_seq)
    addtest_seq_out = addtest_seq_out.follows(extest_seq_out)
    addtest_seq_out = addtest_seq_out.follows(task_testdata_exp)
    params_addtest_seq_out = list()

    task_seq_model = pipe.merge(task_func=touch_checkfile,
                                name='task_seq_model',
                                input=output_from(trainmodel_expstat_seq,
                                                  extrain_seq_out,
                                                  addtrain_seq_out,
                                                  apply_expstat_seq,
                                                  extest_seq_out,
                                                  addtest_seq_out),
                                output=os.path.join(workbase, 'task_seq_model.chk'))

    # ========================================
    # subtask: train model on enhancers only

    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'trainmodel_expstat_enh').replace('\n', ' ')
    trainmodel_expstat_enh = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_expstat_enh',
                                            input=output_from(mrgtraindataexp_groups),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                            output=os.path.join(dir_train_expstat,
                                                                'enh',
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.rfcls.enh.all.pck'),
                                            extras=[cmd, jobcall])
    trainmodel_expstat_enh = trainmodel_expstat_enh.mkdir(dir_train_expstat, 'enh')
    trainmodel_expstat_enh = trainmodel_expstat_enh.follows(task_seq_model)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_enh_out').replace('\n', ' ')
    extrain_enh_out = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                     name='extrain_enh_out',
                                     input=output_from(trainmodel_expstat_enh),
                                     filter=formatter('(?P<SAMPLE>[\w\.]+).rfcls.enh.all.pck'),
                                     output='{path[0]}/{SAMPLE[0]}.enh-cprob.h5',
                                     extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_enh_out').replace('\n', ' ')
    params_addtrain_enh_out = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                       pattern='*.feat.seqout.h5', allow_none=True),
                                                    collect_full_paths(os.path.join(dir_train_expstat, 'enh'),
                                                                       pattern='*enh-cprob.h5',
                                                                       allow_none=True),
                                                    '.enhout', cmd, jobcall)

    addtrain_enh_out = pipe.files(sci_obj.get_jobf('ins_out'),
                                  params_addtrain_enh_out,
                                  name='addtrain_enh_out')
    addtrain_enh_out = addtrain_enh_out.follows(extrain_enh_out)
    addtrain_enh_out = addtrain_enh_out.follows(task_traindata_exp)
    params_addtrain_enh_out = list()

    # apply enhancer model to test data
    
    cmd = config.get('Pipeline', 'apply_expstat_enh').replace('\n', ' ')
    params_apply_expstat_enh = params_predict_testdata(os.path.join(dir_train_expstat, 'enh'),
                                                       os.path.split(dir_mrg_test_datasets)[0],
                                                       os.path.join(dir_apply_expstat, 'enh', '{groupid}'),
                                                       cmd, jobcall)
    apply_expstat_enh = pipe.files(sci_obj.get_jobf('inpair_out'),
                                   params_apply_expstat_enh,
                                   name='apply_expstat_enh')
    apply_expstat_enh = apply_expstat_enh.mkdir(os.path.join(dir_apply_expstat, 'enh'))
    apply_expstat_enh = apply_expstat_enh.follows(mrgtestdataexp_groups)
    apply_expstat_enh = apply_expstat_enh.follows(trainmodel_expstat_enh)
    params_apply_expstat_enh = list()
    
    # extract sequence model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_enh_out').replace('\n', ' ')
    extest_enh_out = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                  name='extest_enh_out',
                                  input=output_from(apply_expstat_enh),
                                  filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+enh\.all\.json'),
                                  output='{path[0]}/{GID[0]}_{DSET[0]}.enh-cprob.h5',
                                  extras=[cmd, jobcall])
    
    # add enhancer model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_enh_out').replace('\n', ' ')
    params_addtest_enh_out = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                     '*.feat.seqout.h5', allow_none=True),
                                                  collect_full_paths(os.path.join(dir_apply_expstat, 'enh'),
                                                                     '*.enh-cprob.h5', allow_none=True),
                                                  '.enhout', cmd, jobcall)
    addtest_enh_out = pipe.files(sci_obj.get_jobf('ins_out'),
                                 params_addtest_enh_out,
                                 name='addtest_enh_out')
    addtest_enh_out = addtest_enh_out.follows(apply_expstat_enh)
    addtest_enh_out = addtest_enh_out.follows(extest_enh_out)
    addtest_enh_out = addtest_enh_out.follows(task_testdata_exp)
    params_addtest_enh_out = list()
    
    task_enh_model = pipe.merge(task_func=touch_checkfile,
                                name='task_enh_model',
                                input=output_from(trainmodel_expstat_enh,
                                                  extrain_enh_out,
                                                  addtrain_enh_out,
                                                  apply_expstat_enh,
                                                  extest_enh_out,
                                                  addtest_enh_out),
                                output=os.path.join(workbase, 'task_enh_model.chk'))

    # Subtask
    # train augmented signal model to classify genes as active/inactive (TPM >= 1)
    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'trainmodel_expstat_asig').replace('\n', ' ')
    trainmodel_expstat_asig = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                             name='trainmodel_expstat_asig',
                                             input=output_from(addtrain_enh_out),
                                             filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.seqout\.enhout\.h5'),
                                             output=os.path.join(dir_train_expstat,
                                                                 'asig',
                                                                 '{subdir[0][0]}',
                                                                 '{SAMPLE[0]}.rfcls.asig.all.pck'),
                                             extras=[cmd, jobcall])
    trainmodel_expstat_asig = trainmodel_expstat_asig.mkdir(os.path.join(dir_train_expstat, 'asig'))
    trainmodel_expstat_asig = trainmodel_expstat_asig.follows(task_enh_model)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_asig_out').replace('\n', ' ')
    extrain_asig_out = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='extrain_asig_out',
                                      input=output_from(trainmodel_expstat_asig),
                                      filter=formatter('(?P<SAMPLE>[\w\.]+).rfcls.asig.all.pck'),
                                      output='{path[0]}/{SAMPLE[0]}.asig-cprob.h5',
                                      extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_asig_out').replace('\n', ' ')
    params_addtrain_asig_out = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                        pattern='*feat.seqout.enhout.h5', allow_none=True),
                                                     collect_full_paths(os.path.join(dir_train_expstat, 'asig'),
                                                                        pattern='*asig-cprob.h5',
                                                                        allow_none=True),
                                                     '.asigout', cmd, jobcall)

    addtrain_asig_out = pipe.files(sci_obj.get_jobf('ins_out'),
                                   params_addtrain_asig_out,
                                   name='addtrain_asig_out')
    addtrain_asig_out = addtrain_asig_out.follows(extrain_asig_out)
    addtrain_asig_out = addtrain_asig_out.follows(task_seq_model)
    params_addtrain_asig_out = list()

    # apply asig model to test data
    cmd = config.get('Pipeline', 'apply_expstat_asig').replace('\n', ' ')
    params_apply_expstat_asig = params_predict_testdata(os.path.join(dir_train_expstat, 'asig'),
                                                        collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                           pattern='*seqout.enhout.h5',
                                                                           allow_none=True),
                                                        os.path.join(dir_apply_expstat, 'asig', '{groupid}'),
                                                        cmd, jobcall)

    apply_expstat_asig = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    params_apply_expstat_asig,
                                    name='apply_expstat_asig')
    apply_expstat_asig = apply_expstat_asig.mkdir(os.path.join(dir_apply_expstat, 'asig'))
    apply_expstat_asig = apply_expstat_asig.follows(addtrain_asig_out)
    params_apply_expstat_asig = list()

    # extract asig model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_asig_cls_out').replace('\n', ' ')
    extest_asig_cls_out = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                       name='extest_asig_cls_out',
                                       input=output_from(apply_expstat_asig),
                                       filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+asig\.all\.json'),
                                       output='{path[0]}/{GID[0]}_{DSET[0]}.asig-cprob.h5',
                                       extras=[cmd, jobcall])

    # add asig model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_asig_cls_out').replace('\n', ' ')
    params_addtest_asig_cls_out = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                          '*seqout.enhout.h5', allow_none=True),
                                                       collect_full_paths(os.path.join(dir_apply_expstat, 'asig'),
                                                                          '*.asig-cprob.h5', allow_none=True),
                                                       '.asigout', cmd, jobcall)

    addtest_asig_cls_out = pipe.files(sci_obj.get_jobf('ins_out'),
                                      params_addtest_asig_cls_out,
                                      name='addtest_asig_cls_out')
    addtest_asig_cls_out = addtest_asig_cls_out.follows(apply_expstat_asig)
    addtest_asig_cls_out = addtest_asig_cls_out.follows(extest_asig_cls_out)
    params_addtest_asig_cls_out = list()

    cmd = config.get('Pipeline', 'summ_perf_status')
    dir_summary = os.path.join(workbase, 'task_summarize')
    summ_perf_status = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                                  name='summ_perf_status',
                                  input=output_from(apply_expstat_seq, apply_expstat_asig),
                                  output=os.path.join(dir_summary, 'agg_expstat_est.h5'),
                                  extras=[cmd, jobcall])
    summ_perf_status = summ_perf_status.mkdir(dir_summary)
    summ_perf_status = summ_perf_status.active_if(False)

    task_asig_cls_model = pipe.merge(task_func=touch_checkfile,
                                     name='task_asig_cls_model',
                                     input=output_from(trainmodel_expstat_asig,
                                                       extrain_asig_out,
                                                       addtrain_asig_out,
                                                       apply_expstat_asig,
                                                       extest_asig_cls_out,
                                                       addtest_asig_cls_out,
                                                       summ_perf_status),
                                     output=os.path.join(workbase, 'task_asig_cls_model.chk'))

    # above: status prediction
    # =======================================
    # below: rank/value prediction

    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # train augmented signal model to regress gene expression rank
    dir_train_exprank = os.path.join(dir_task_trainmodel_exp, 'sub_rank')

    cmd = config.get('Pipeline', 'trainmodel_exprank_all').replace('\n', ' ')
    trainmodel_exprank_all = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_exprank_all',
                                            input=output_from(addtrain_asig_out),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.asigout\.h5'),
                                            output=os.path.join(dir_train_exprank,
                                                                'all',
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.rfreg.rk.all.pck'),
                                            extras=[cmd, jobcall])
    trainmodel_exprank_all = trainmodel_exprank_all.mkdir(os.path.join(dir_train_exprank, 'all'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_rank_all').replace('\n', ' ')
    extrain_rank_all = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='extrain_rank_all',
                                      input=output_from(trainmodel_exprank_all),
                                      filter=formatter('(?P<SAMPLE>[\w\.]+).rfreg.rk.all.pck'),
                                      output='{path[0]}/{SAMPLE[0]}.rank-all.h5',
                                      extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_rank_all').replace('\n', ' ')
    params_addtrain_rank_all = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                        pattern='*feat.h5', allow_none=True),
                                                     collect_full_paths(os.path.join(dir_train_exprank, 'all'),
                                                                        pattern='*rank-all.h5',
                                                                        allow_none=True),
                                                     '.rk-all', cmd, jobcall)

    addtrain_rank_all = pipe.files(sci_obj.get_jobf('ins_out'),
                                   params_addtrain_rank_all,
                                   name='addtrain_rank_all')
    addtrain_rank_all = addtrain_rank_all.follows(extrain_rank_all)
    addtrain_rank_all = addtrain_rank_all.follows(task_asig_cls_model)
    params_addtrain_rank_all = list()

    # apply rank model to test data
    dir_apply_exprank = os.path.join(dir_task_applymodel_exp, 'sub_rank')

    cmd = config.get('Pipeline', 'apply_exprank_all').replace('\n', ' ')
    params_apply_exprank_all = params_predict_testdata(os.path.join(dir_train_exprank, 'all'),
                                                       collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                          pattern='*asigout.h5',
                                                                          allow_none=True),
                                                       os.path.join(dir_apply_exprank, 'all', '{groupid}'),
                                                       cmd, jobcall)

    apply_exprank_all = pipe.files(sci_obj.get_jobf('inpair_out'),
                                   params_apply_exprank_all,
                                   name='apply_exprank_all')
    apply_exprank_all = apply_exprank_all.mkdir(os.path.join(dir_apply_exprank, 'all'))
    apply_exprank_all = apply_exprank_all.follows(addtrain_rank_all)
    params_apply_exprank_all = list()

    # extract asig model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_rank_all').replace('\n', ' ')
    extest_rank_all = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                   name='extest_rank_all',
                                   input=output_from(apply_exprank_all),
                                   filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+rk\.all\.json'),
                                   output='{path[0]}/{GID[0]}_{DSET[0]}.rk-all.h5',
                                   extras=[cmd, jobcall])

    # add extracted asig model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_rank_all').replace('\n', ' ')
    params_addtest_rank_all = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                      '*.feat.h5', allow_none=True),
                                                   collect_full_paths(os.path.join(dir_apply_exprank, 'all'),
                                                                      '*.rk-all.h5', allow_none=True),
                                                   '.rkall', cmd, jobcall)
    addtest_rank_all = pipe.files(sci_obj.get_jobf('ins_out'),
                                  params_addtest_rank_all,
                                  name='addtest_rank_all')
    addtest_rank_all = addtest_rank_all.follows(apply_exprank_all)
    addtest_rank_all = addtest_rank_all.follows(extest_rank_all)
    params_addtest_rank_all = list()

    task_rank_all_model = pipe.merge(task_func=touch_checkfile,
                                     name='task_rank_all_model',
                                     input=output_from(trainmodel_exprank_all,
                                                       extrain_rank_all,
                                                       addtrain_rank_all,
                                                       apply_exprank_all,
                                                       extest_rank_all,
                                                       addtest_rank_all),
                                     output=os.path.join(workbase, 'task_rank_all_model.chk'))

    # Subtask
    # Using sample ranks as additional feature, make final prediction of actual expression level (TPM)

    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_train_explevel = os.path.join(dir_task_trainmodel_exp, 'sub_level')

    cmd = config.get('Pipeline', 'trainmodel_explevel_all').replace('\n', ' ')
    trainmodel_explevel_all = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                             name='trainmodel_explevel_all',
                                             input=output_from(addtrain_rank_all),
                                             filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.rk-all\.h5'),
                                             output=os.path.join(dir_train_explevel,
                                                                 'all',
                                                                 '{subdir[0][0]}',
                                                                 '{SAMPLE[0]}.rfreg.tpm.all.pck'),
                                             extras=[cmd, jobcall])
    trainmodel_explevel_all = trainmodel_explevel_all.mkdir(os.path.join(dir_train_explevel, 'all'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_level_all').replace('\n', ' ')
    extrain_level_all = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='extrain_level_all',
                                       input=output_from(trainmodel_explevel_all),
                                       filter=formatter('(?P<SAMPLE>[\w\.]+).rfreg.tpm.all.pck'),
                                       output='{path[0]}/{SAMPLE[0]}.level-all.h5',
                                       extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_level_all').replace('\n', ' ')
    params_addtrain_level_all = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                         pattern='*feat.rk-all.h5', allow_none=True),
                                                      collect_full_paths(os.path.join(dir_train_explevel, 'all'),
                                                                         pattern='*level-all.h5',
                                                                         allow_none=True),
                                                      '.lvl-all', cmd, jobcall)

    addtrain_level_all = pipe.files(sci_obj.get_jobf('ins_out'),
                                    params_addtrain_level_all,
                                    name='addtrain_level_all')
    addtrain_level_all = addtrain_level_all.follows(extrain_level_all)
    addtrain_level_all = addtrain_level_all.follows(task_rank_all_model)
    params_addtrain_level_all = list()

    # apply exp level model to test data
    dir_apply_explevel = os.path.join(dir_task_applymodel_exp, 'sub_level')

    cmd = config.get('Pipeline', 'apply_explevel_all').replace('\n', ' ')
    params_apply_explevel_all = params_predict_testdata(os.path.join(dir_train_explevel, 'all'),
                                                        collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                           pattern='*rkall.h5',
                                                                           allow_none=True),
                                                        os.path.join(dir_apply_explevel, 'all', '{groupid}'),
                                                        cmd, jobcall)

    apply_explevel_all = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    params_apply_explevel_all,
                                    name='apply_explevel_all')
    apply_explevel_all = apply_explevel_all.mkdir(os.path.join(dir_apply_explevel, 'all'))
    apply_explevel_all = apply_explevel_all.follows(addtrain_level_all)
    params_apply_explevel_all = list()

    # extract asig model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_level_all').replace('\n', ' ')
    extest_level_all = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                    name='extest_level_all',
                                    input=output_from(apply_explevel_all),
                                    filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+tpm\.all\.json'),
                                    output='{path[0]}/{GID[0]}_{DSET[0]}.level-all.h5',
                                    extras=[cmd, jobcall])

    # add extracted asig model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_level_all').replace('\n', ' ')
    params_addtest_level_all = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                       '*.feat.rkall.h5', allow_none=True),
                                                    collect_full_paths(os.path.join(dir_apply_explevel, 'all'),
                                                                       '*.level-all.h5', allow_none=True),
                                                    '.lvl-all', cmd, jobcall)

    addtest_level_all = pipe.files(sci_obj.get_jobf('ins_out'),
                                   params_addtest_level_all,
                                   name='addtest_level_all')
    addtest_level_all = addtest_level_all.follows(apply_explevel_all)
    addtest_level_all = addtest_level_all.follows(extest_level_all)
    params_addtest_level_all = list()

    task_level_all_model = pipe.merge(task_func=touch_checkfile,
                                      name='task_level_all_model',
                                      input=output_from(trainmodel_explevel_all,
                                                        extrain_level_all,
                                                        addtrain_level_all,
                                                        apply_explevel_all,
                                                        extest_level_all,
                                                        addtest_level_all),
                                      output=os.path.join(workbase, 'task_level_all_model.chk'))

    # above: rank/level prediction for all samples
    # ===============================================
    # below: rank/level prediction for active subset

    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # train augmented signal model to regress gene expression rank
    dir_train_exprank = os.path.join(dir_task_trainmodel_exp, 'sub_rank')

    cmd = config.get('Pipeline', 'trainmodel_exprank_act').replace('\n', ' ')
    trainmodel_exprank_act = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_exprank_act',
                                            input=output_from(addtrain_asig_out),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.asigout\.h5'),
                                            output=os.path.join(dir_train_exprank,
                                                                'act',
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.rfreg.rk.act.pck'),
                                            extras=[cmd, jobcall])
    trainmodel_exprank_act = trainmodel_exprank_act.mkdir(os.path.join(dir_train_exprank, 'act'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_rank_act').replace('\n', ' ')
    extrain_rank_act = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                      name='extrain_rank_act',
                                      input=output_from(trainmodel_exprank_act),
                                      filter=formatter('(?P<SAMPLE>[\w\.]+).rfreg.rk.act.pck'),
                                      output='{path[0]}/{SAMPLE[0]}.rank-act.h5',
                                      extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_rank_act').replace('\n', ' ')
    params_addtrain_rank_act = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                        pattern='*feat.h5', allow_none=True),
                                                     collect_full_paths(os.path.join(dir_train_exprank, 'act'),
                                                                        pattern='*rank-act.h5',
                                                                        allow_none=True),
                                                     '.rk-act', cmd, jobcall)

    addtrain_rank_act = pipe.files(sci_obj.get_jobf('ins_out'),
                                   params_addtrain_rank_act,
                                   name='addtrain_rank_act')
    addtrain_rank_act = addtrain_rank_act.follows(extrain_rank_act)
    addtrain_rank_act = addtrain_rank_act.follows(task_asig_cls_model)
    params_addtrain_rank_act = list()

    # apply rank model to test data
    dir_apply_exprank = os.path.join(dir_task_applymodel_exp, 'sub_rank')

    cmd = config.get('Pipeline', 'apply_exprank_act').replace('\n', ' ')
    params_apply_exprank_act = params_predict_testdata(os.path.join(dir_train_exprank, 'act'),
                                                       collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                          pattern='*asigout.h5',
                                                                          allow_none=True),
                                                       os.path.join(dir_apply_exprank, 'act', '{groupid}'),
                                                       cmd, jobcall)

    apply_exprank_act = pipe.files(sci_obj.get_jobf('inpair_out'),
                                   params_apply_exprank_act,
                                   name='apply_exprank_act')
    apply_exprank_act = apply_exprank_act.mkdir(os.path.join(dir_apply_exprank, 'act'))
    apply_exprank_act = apply_exprank_act.follows(addtrain_rank_act)
    params_apply_exprank_act = list()

    # extract asig model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_rank_act').replace('\n', ' ')
    extest_rank_act = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                   name='extest_rank_act',
                                   input=output_from(apply_exprank_act),
                                   filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+rk\.act\.json'),
                                   output='{path[0]}/{GID[0]}_{DSET[0]}.rk-act.h5',
                                   extras=[cmd, jobcall])

    # add extracted asig model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_rank_act').replace('\n', ' ')
    params_addtest_rank_act = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                      '*.feat.h5', allow_none=True),
                                                   collect_full_paths(os.path.join(dir_apply_exprank, 'act'),
                                                                      '*.rk-act.h5', allow_none=True),
                                                   '.rkact', cmd, jobcall)

    addtest_rank_act = pipe.files(sci_obj.get_jobf('ins_out'),
                                  params_addtest_rank_act,
                                  name='addtest_rank_act')
    addtest_rank_act = addtest_rank_act.follows(apply_exprank_act)
    addtest_rank_act = addtest_rank_act.follows(extest_rank_act)
    params_addtest_rank_act = list()

    task_rank_act_model = pipe.merge(task_func=touch_checkfile,
                                     name='task_rank_act_model',
                                     input=output_from(trainmodel_exprank_act,
                                                       extrain_rank_act,
                                                       addtrain_rank_act,
                                                       apply_exprank_act,
                                                       extest_rank_act,
                                                       addtest_rank_act),
                                     output=os.path.join(workbase, 'task_rank_act_model.chk'))

    # Subtask
    # Using sample ranks as additional feature, make final prediction of actual expression level (TPM)

    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    dir_train_explevel = os.path.join(dir_task_trainmodel_exp, 'sub_level')

    cmd = config.get('Pipeline', 'trainmodel_explevel_act').replace('\n', ' ')
    trainmodel_explevel_act = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                             name='trainmodel_explevel_act',
                                             input=output_from(addtrain_rank_act),
                                             filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.rk-act\.h5'),
                                             output=os.path.join(dir_train_explevel,
                                                                 'act',
                                                                 '{subdir[0][0]}',
                                                                 '{SAMPLE[0]}.rfreg.tpm.act.pck'),
                                             extras=[cmd, jobcall])
    trainmodel_explevel_act = trainmodel_explevel_act.mkdir(os.path.join(dir_train_explevel, 'act'))

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'extrain_level_act').replace('\n', ' ')
    extrain_level_act = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='extrain_level_act',
                                       input=output_from(trainmodel_explevel_act),
                                       filter=formatter('(?P<SAMPLE>[\w\.]+).rfreg.tpm.act.pck'),
                                       output='{path[0]}/{SAMPLE[0]}.level-act.h5',
                                       extras=[cmd, jobcall])

    cmd = config.get('Pipeline', 'addtrain_level_act').replace('\n', ' ')
    params_addtrain_level_act = match_prior_traindata(collect_full_paths(os.path.split(dir_mrg_train_datasets)[0],
                                                                         pattern='*feat.rk-act.h5',
                                                                         allow_none=True),
                                                      collect_full_paths(os.path.join(dir_train_explevel, 'act'),
                                                                         pattern='*level-act.h5',
                                                                         allow_none=True),
                                                      '.lvl-act', cmd, jobcall)

    addtrain_level_act = pipe.files(sci_obj.get_jobf('ins_out'),
                                    params_addtrain_level_act,
                                    name='addtrain_level_act')
    addtrain_level_act = addtrain_level_act.follows(extrain_level_act)
    addtrain_level_act = addtrain_level_act.follows(task_rank_act_model)
    params_addtrain_level_act = list()

    # apply exp level model to test data
    dir_apply_explevel = os.path.join(dir_task_applymodel_exp, 'sub_level')

    cmd = config.get('Pipeline', 'apply_explevel_act').replace('\n', ' ')
    params_apply_explevel_act = params_predict_testdata(os.path.join(dir_train_explevel, 'act'),
                                                        collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                           pattern='*rkact.h5',
                                                                           allow_none=True),
                                                        os.path.join(dir_apply_explevel, 'act', '{groupid}'),
                                                        cmd, jobcall)

    apply_explevel_act = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    params_apply_explevel_act,
                                    name='apply_explevel_act')
    apply_explevel_act = apply_explevel_act.mkdir(os.path.join(dir_apply_explevel, 'act'))
    apply_explevel_act = apply_explevel_act.follows(addtrain_level_act)
    params_apply_explevel_act = list()

    # extract asig model predictions for test data and add to test datasets
    cmd = config.get('Pipeline', 'extest_level_act').replace('\n', ' ')
    extest_level_act = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                                    name='extest_level_act',
                                    input=output_from(apply_explevel_act),
                                    filter=formatter('(?P<GID>G[0-9]{4})_.+(?P<DSET>data-[a-zA-Z0-9-]+)_model.+tpm\.act\.json'),
                                    output='{path[0]}/{GID[0]}_{DSET[0]}.level-act.h5',
                                    extras=[cmd, jobcall])

    # add extracted asig model predictions to test datasets
    cmd = config.get('Pipeline', 'addtest_level_act').replace('\n', ' ')
    params_addtest_level_act = match_prior_testdata(collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
                                                                       '*.feat.rkact.h5',
                                                                       allow_none=True),
                                                    collect_full_paths(os.path.join(dir_apply_explevel, 'act'),
                                                                       '*.level-act.h5', allow_none=True),
                                                    '.lvl-act', cmd, jobcall)

    addtest_level_act = pipe.files(sci_obj.get_jobf('ins_out'),
                                   params_addtest_level_act,
                                   name='addtest_level_act')
    addtest_level_act = addtest_level_act.follows(apply_explevel_act)
    addtest_level_act = addtest_level_act.follows(extest_level_act)
    params_addtest_level_act = list()

    task_level_act_model = pipe.merge(task_func=touch_checkfile,
                                      name='task_level_act_model',
                                      input=output_from(trainmodel_explevel_act,
                                                        extrain_level_act,
                                                        addtrain_level_act,
                                                        apply_explevel_act,
                                                        extest_level_act,
                                                        addtest_level_act),
                                      output=os.path.join(workbase, 'task_level_act_model.chk'))

    # training and applying models done

    run_all = pipe.merge(task_func=touch_checkfile,
                         name='run_task_all',
                         input=output_from(task_sigmap, task_traindata_exp, task_testdata_exp,
                                           task_seq_model, task_asig_cls_model,
                                           task_rank_all_model, task_rank_act_model,
                                           task_level_all_model, task_level_act_model),
                         output=os.path.join(workbase, 'run_task_all.chk'))

    # aggregate training and testing performances
    cmd = config.get('Pipeline', 'summ_perf').replace('\n', ' ')
    summarize_perf = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                                name='summarize_perf',
                                input=output_from(trainmodel_expstat_seq, trainmodel_expstat_asig,
                                                  trainmodel_exprank_all, trainmodel_exprank_act,
                                                  trainmodel_explevel_all, trainmodel_explevel_act,
                                                  apply_expstat_seq, apply_expstat_asig,
                                                  apply_exprank_all, apply_exprank_act,
                                                  apply_explevel_all, apply_explevel_act),
                                output=os.path.join(dir_summary, 'train_test_perf_agg.h5'),
                                extras=[cmd, jobcall])
    summarize_perf = summarize_perf.mkdir(dir_summary)
    summarize_perf = summarize_perf.follows(run_all)

    cmd = config.get('Pipeline', 'summ_expest')
    summarize_exp_est = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                       name='summarize_exp_est',
                                       input=output_from(summarize_perf),
                                       filter=formatter(),
                                       output=os.path.join(dir_summary, 'exp_est_agg.h5'),
                                       extras=[cmd, jobcall])

    # # summarize statistics
    # dir_summarize = os.path.join(workbase, 'task_summarize')
    # cmd = config.get('Pipeline', 'summtt_hg19').replace('\n', ' ')
    # summtt_hg19 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                              name='summtt_hg19',
    #                              input=output_from(task_level_act_model),
    #                              filter=formatter(),
    #                              output=os.path.join(dir_summarize, 'hg19_train_test_stat.tsv'),
    #                              extras=[cmd, jobcall]).mkdir(dir_summarize)
    #
    # cmd = config.get('Pipeline', 'summtt_mm9').replace('\n', ' ')
    # summtt_mm9 = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
    #                             name='summtt_mm9',
    #                             input=output_from(task_level_act_model),
    #                             filter=formatter(),
    #                             output=os.path.join(dir_summarize, 'mm9_train_test_stat.tsv'),
    #                             extras=[cmd, jobcall]).mkdir(dir_summarize)

    # cmd = config.get('Pipeline', 'collect_train_metrics').replace('\n', ' ')
    # collect_train_metrics = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
    #                                    name='train_metrics',
    #                                    input=output_from(trainmodel_expone_seq, trainmodel_expone_psig,
    #                                                      trainmodel_expone_full, trainmodel_expvalone_seq,
    #                                                      trainmodel_expvalone_psig_all, trainmodel_expvalone_psig_sub,
    #                                                      trainmodel_expvalone_full),
    #                                    output=os.path.join(dir_task_trainmodel_exp, 'models_train_metrics.h5'),
    #                                    extras=[cmd, jobcall])
    #

    #
    # # ==============================
    # # Major task: apply models to predict gene status (on/off) and expression level
    # #
    #
    # # NB: this does not perform permutation testing at the moment (takes a loooong time)
    # # iff that is enabled, change job config to CVParallelJobConfig!
    # dir_task_applymodels_exp = os.path.join(workbase, 'task_applymodels_exp')
    # sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # # ==================================================
    # # annotations needed to categorize test output
    # # matchings json from above
    # # groupfile from above
    # datasetfile = config.get('Annotations', 'dsetfile')
    # groupings = prep_metadata(groupfile, 'gid')
    # datasets = prep_metadata(datasetfile, 'id')
    # cellmatches = js.load(open(matchings, 'r'))
    # # ===================================================
    #
    # # subtask: predict gene status (on/off) for subset TPM >= 1
    # dir_apply_statusone = os.path.join(dir_task_applymodels_exp, 'sub_status', 'active_geq1', '{groupid}')
    # cmd = config.get('Pipeline', 'applyexpstatone').replace('\n', ' ')
    # params_applyexpstatone = params_predict_testdata(dir_exp_statone, os.path.split(dir_mrg_test_datasets)[0],
    #                                                  dir_apply_statusone, cmd, jobcall)
    #
    # statone_out = os.path.join(dir_task_applymodels_exp, 'eval_models_stat-one.json')
    # _ = annotate_test_output(params_applyexpstatone, datasets, groupings, cellmatches, 'Status >= 1', statone_out)
    #
    # applyexpstatone = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                              params_applyexpstatone,
    #                              name='applyexpstatone')\
    #                                 .mkdir(os.path.split(dir_apply_statusone)[0])\
    #                                 .follows(mrgtestdataexp_groups).follows(task_trainmodel_exp)
    # params_applyexpstatone = list()
    #

    #
    # cmd = config.get('Pipeline', 'applyexpstatone').replace('\n', ' ')
    # params_applyexpstatone_prior = params_predict_testdata_prior(dir_exp_statone, os.path.split(dir_mrg_test_datasets)[0],
    #                                                              dir_apply_statusone, cmd, jobcall)
    #
    # applyexpstatone_prior = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                                    params_applyexpstatone_prior,
    #                                    name='applyexpstatone_prior')\
    #                                     .mkdir(os.path.split(dir_apply_statusone)[0])\
    #                                     .follows(mrgtestdataexp_groups)\
    #                                     .follows(task_trainmodel_exp)\
    #                                     .follows(addcprob_testdata)
    #
    # # ==========
    # # subtask: predict gene expression rank on subset: predicted status 1
    # dir_apply_exprank_sub = os.path.join(dir_task_applymodels_exp, 'sub_rank', 'active_geq1', '{groupid}')
    # cmd = config.get('Pipeline', 'applyexprank_sub').replace('\n', ' ')
    # params_applyexprank_sub = params_predict_testdata(collect_full_paths(dir_exp_valone_sub,
    #                                                                      pattern='*psig*.pck',
    #                                                                      allow_none=True),
    #                                                   collect_full_paths(os.path.split(dir_mrg_test_datasets)[0],
    #                                                                      pattern='*feat.cprob.h5',
    #                                                                      allow_none=True),
    #                                                   dir_apply_exprank_sub,
    #                                                   cmd, jobcall)
    #
    # #valonetrue_out = os.path.join(dir_task_applymodels_exp, 'eval_models_val-one-true.json')
    # #_ = annotate_test_output(params_applyexpvalone_true, datasets, groupings, cellmatches, 'Value TPM >= 1', valonetrue_out)
    #
    # applyexprank_sub = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                               params_applyexprank_sub,
    #                               name='applyexpranksub')\
    #                               .mkdir(os.path.split(dir_apply_exprank_sub)[0])\
    #                               .follows(mrgtestdataexp_groups).follows(task_trainmodel_exp)
    # params_applyexpvalone_true = list()
    #
    # # subtask: predict gene expression level (using predicted TPM-based status)
    # dir_apply_valueonepred = os.path.join(dir_task_applymodels_exp, 'sub_value', 'pred_status', 'active_geq1', '{groupid}')
    # cmd = config.get('Pipeline', 'applyexpvalonepred').replace('\n', ' ')
    # load_metadata = os.path.split(dir_apply_statusone)[0]
    # params_applyexpvalone_est = params_predict_testdata(dir_exp_valone_sub, os.path.split(dir_mrg_test_datasets)[0],
    #                                                     dir_apply_valueonepred, cmd, jobcall, True, load_metadata)
    #
    # valoneest_out = os.path.join(dir_task_applymodels_exp, 'eval_models_val-one-est.json')
    # _ = annotate_test_output(params_applyexpvalone_est, datasets, groupings, cellmatches, 'Value Est. >= 1', valoneest_out)
    #
    # applyexpvalonepred = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                                 params_applyexpvalone_est,
    #                                 name='applyexpvalonepred')\
    #                                     .mkdir(os.path.split(dir_apply_valueonepred)[0])\
    #                                     .follows(mrgtestdataexp_groups).follows(task_trainmodel_exp)
    # params_applyexpvalone_est = list()
    #
    # cmd = config.get('Pipeline', 'mrgtestmd')
    # mrgtestmd = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
    #                        name='mrgtestmd',
    #                        input=[statone_out, valonetrue_out, valoneest_out],
    #                        output=os.path.join(dir_task_applymodels_exp, 'eval_models_test_merged.json'),
    #                        extras=[cmd, jobcall]).follows(applyexpvalonepred)
    #
    # run_task_applymodels_exp = pipe.merge(task_func=touch_checkfile,
    #                                       name='task_applymodels_exp',
    #                                       input=output_from(applyexpstatone,
    #                                                         applyexpvalone,
    #                                                         applyexpvalonepred,
    #                                                         mrgtestmd),
    #                                       output=os.path.join(dir_task_applymodels_exp, 'task_applymodels_exp.chk'))
    #
    # #
    # # END: major task apply gene expression models
    # # ==============================
    #
    #
    #
    #
    #

    return pipe
