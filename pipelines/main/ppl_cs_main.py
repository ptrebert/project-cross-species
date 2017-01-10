# coding=utf-8

import os as os
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
    fpaths = collect_full_paths(folder, '*.map.h5', False)
    mapfiles = list()
    for fp in fpaths:
        target, query = os.path.basename(fp).split('.')[0].split('_to_')
        if target in targets and query in queries:
            mapfiles.append({'target': target, 'query': query, 'path': fp})
    assert mapfiles, 'No map files collected from path: {}'.format(folder)
    return mapfiles


def collect_roi_files(folder):
    """
    :param folder:
    :return:
    """
    fpaths = collect_full_paths(folder, '*.h5', False)
    roifiles = []
    for fp in fpaths:
        annotid, regtype, _ = os.path.basename(fp).split('.')
        assembly = annotid.split('_')[1]
        roifiles.append({'assembly': assembly, 'regtype': regtype, 'path': fp})
    assert roifiles, 'No ROI files annotated from path: {}'.format(folder)
    return roifiles


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
    with open(groupings, 'r', newline='') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        for r in rows:
            match_types[r['comment']] = r['type']
            foo, bar = r['comment'].split('-')
            match_types[bar + '-' + foo] = r['type']
            # self-match can happen for liver, kidney
            match_types[foo + '-' + foo] = 'pos'
            match_types[bar + '-' + bar] = 'pos'
            p1, p2, h, d, clibs = select_compatible_libs(epibundles[r['partner1']], epibundles[r['partner2']])
            clibs = '-'.join(clibs)
            p1_params = ' '.join([x['param'] for x in p1])
            p2_params = ' '.join([x['param'] for x in p2])
            gid = r['gid']
            extgid = gid + str(h) + str(d)
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
    fpaths = collect_full_paths(folder, '*/E0*.h5')
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


def annotate_training_groups(mapfiles, roifiles, epidir, groupfile, outraw, cmd, jobcall):
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
                if cell1 not in ['liver', 'kidney'] and cell2 not in ['liver', 'kidney']:
                    continue
            for roi in compat_roi:
                outfolder = outraw.format(**{'target': target, 'query': query})
                outfile = '_'.join([partners['extgid'], partners['partner1'][0]['EID'],
                                    assm1, partners['partner1'][0]['cell']])
                outfile += '.to.{}.{}.h5'.format(query, roi['regtype'])
                outpath = os.path.join(outfolder, outfile)
                assert outpath not in uniq, 'Created duplicate: {}'.format(outpath)
                uniq.add(outpath)
                params = {'target': target, 'query': query, 'genome': assm1,
                          'regtype': roi['regtype'], 'mapfile': mapf['path'],
                          'datafiles': partners['params1'], 'info': partners['extgid']}
                tmp = cmd.format(**params)
                arglist.append([roi['path'], outpath, tmp, jobcall])
        # same for other group partner
        compat_roi = [roi for roi in roifiles if roi['assembly'] == assm2]
        assm2_maps = [mapf for mapf in mapfiles if mapf['target'] == assm2]
        for mapf in assm2_maps:
            target, query = mapf['target'], mapf['query']
            if query not in ['mm9', 'hg19']:
                if cell1 not in ['liver', 'kidney'] and cell2 not in ['liver', 'kidney']:
                    continue
            for roi in compat_roi:
                outfolder = outraw.format(**{'target': target, 'query': query})
                outfile = '_'.join([partners['extgid'], partners['partner2'][0]['EID'],
                                    assm2, partners['partner2'][0]['cell']])
                outfile += '.to.{}.{}.h5'.format(query, roi['regtype'])
                outpath = os.path.join(outfolder, outfile)
                assert outpath not in uniq, 'Created duplicate: {}'.format(outpath)
                uniq.add(outpath)
                params = {'target': target, 'query': query, 'genome': assm2,
                          'regtype': roi['regtype'], 'mapfile': mapf['path'],
                          'datafiles': partners['params2'], 'info': partners['extgid']}
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
    serialize = dict()
    for key, vals in matchings.items():
        serialize[key] = list(vals)
    timestr = dt.datetime.now().strftime('%A_%Y-%m-%d_%H:%M:%S')
    filepath = os.path.join(os.path.dirname(groupfile), 'cellmatches_ro.json')
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'matchings': {}, 'matchtypes': {}}, dump, indent=1, sort_keys=True)
    rewrite = False
    with open(filepath, 'r') as dumped:
        cmi = js.load(dumped)['matchings']
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
            sigfiles = fnm.filter(inputfiles, '*/E0*_{}_*.h5'.format(trg))
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
    expfiles = collect_full_paths(expdir, '*/T0*.h5', False, False)
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
        assert len(files) == 4, 'Unexpected number of feature files {}: {}'.format(group, files)
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
                foo = fnm.filter(expfiles, '*/T0*_{}_{}_mRNA*'.format(query, mc))
                match_exp.extend(foo)
        else:
            pat = group.split('.')[0].split('_', 2)[-1]
            match_exp = fnm.filter(expfiles, '*/T0*_{}_mRNA*'.format(pat))
        try:
            assert match_exp, 'No matched expression files for group {} / pattern {} / {}'.format(group, pat, files)
        except AssertionError as ae:
            if '_ESE14' in pat:
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


def annotate_test_datasets(mapfiles, roifiles, mapepidir, expfiles, groupfile, outraw, cmd, jobcall):
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
                                assert outpath not in uniq, 'Created duplicate: {} / {}'.format(outpath, mapf)
                                uniq.add(outpath)
                                params = {'target': target, 'query': query, 'genome': query,
                                          'regtype': roi['regtype'], 'mapfile': mapf['path'],
                                          'datafiles': datafiles, 'info': extgid}
                                tmp = cmd.format(**params)
                                arglist.append([roi['path'], outpath, tmp, jobcall])
    if mapfiles:
        assert arglist, 'No argument list for building test datasets created'
    return arglist


def params_predict_testdata(modelroot, dataroot, outroot, cmd, jobcall, addmd=False, mdpath=''):
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
    for model in all_models:
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
        assert datasets, 'No test datasets for group: {}'.format(groupid)
        for dat in datasets:
            fp2, fn2 = os.path.split(dat)
            assert fn2.startswith(groupid), 'Group mismatch: {} vs {}'.format(groupid, fn2)
            dqry, fr, dtrg = (os.path.split(fp2)[-1]).split('_')
            assert fr == 'from', 'Unexpected folder structure: {}'.format(fp2)
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
    if all_models:
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
    workbase = config.get('Pipeline', 'workdir')

    # collect raw input data in HDF format
    inputfiles = []
    tmpf = collect_full_paths(dir_indata, '*.h5', False)
    inputfiles.extend(tmpf)

    mapfiles = []
    tmpf = collect_full_paths(dir_maps, '*.map.h5', False)
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
    # dir_task_sigcorr = os.path.join(workbase, 'task_signal_correlation')
    #
    # # Subtask
    # # compute correlations of signal in regions of interest
    # sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # cmd = config.get('Pipeline', 'corrsigroi').replace('\n', ' ')
    # dir_sub_sigcorr_roi = os.path.join(dir_task_sigcorr, 'sub_roi')
    # corrsigroi = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                         make_srcsig_pairs(inputfiles, use_targets,
    #                                           get_generic_roifiles(config.get('Pipeline', 'refbase')),
    #                                           dir_sub_sigcorr_roi, cmd, jobcall),
    #                         name='corrsigroi').mkdir(dir_sub_sigcorr_roi).active_if(False)
    #
    # run_task_sigcorr = pipe.merge(task_func=touch_checkfile,
    #                               name='run_task_sigcorr',
    #                               input=output_from(corrsigroi),
    #                               output=os.path.join(dir_task_sigcorr, 'run_task_sigcorr.chk')).active_if(False)

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

    # # Subtask
    # # compute correlations of mapped to true signal across species in regions of interest
    # cmd = config.get('Pipeline', 'corrmaproi').replace('\n', ' ')
    # dir_sub_mapext_corr = os.path.join(dir_task_sigmap, 'sub_mapext', 'corr_roi')
    # corrmaproi = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                         make_corr_pairs(inputfiles, use_targets, dir_sub_signal,
    #                                         get_generic_roifiles(config.get('Pipeline', 'refbase')),
    #                                         dir_sub_mapext_corr, cmd, jobcall),
    #                         name='corrmaproi').mkdir(dir_sub_mapext_corr).follows(mapsig).active_if(False)

    task_sigmap = pipe.merge(task_func=touch_checkfile,
                             name='task_sigmap',
                             input=output_from(mapsig),
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
    groupfile = config.get('Annotations', 'groupfile')
    cmd = config.get('Pipeline', 'traindataexp').replace('\n', ' ')

    traindataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                     annotate_training_groups(mapfiles, roifiles, dir_indata,
                                                              groupfile, dir_subtask_traindata_exp_groups,
                                                              cmd, jobcall),
                                     name='traindataexp_groups').mkdir(os.path.split(dir_subtask_traindata_exp_groups)[0])

    dir_mrg_train_datasets = os.path.join(dir_task_traindata_exp, 'train_datasets', '{target}_to_{query}')
    cmd = config.get('Pipeline', 'mrgtraindataexp').replace('\n', ' ')
    mrgtraindataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                        merge_augment_featfiles(os.path.split(dir_subtask_traindata_exp_groups)[0],
                                                                dir_indata, dir_mrg_train_datasets, cmd, jobcall),
                                        name='mrgtraindataexp_groups').mkdir(os.path.split(dir_mrg_train_datasets)[0]).follows(traindataexp_groups)

    task_traindata_exp = pipe.merge(task_func=touch_checkfile,
                                    name='task_traindata_exp',
                                    input=output_from(traindataexp_groups, mrgtraindataexp_groups),
                                    output=os.path.join(dir_task_traindata_exp, 'run_task_traindata_exp.chk'))

    #
    # END: major task generate training data
    # ==============================

    # ==================================
    # Major task: train model for gene expression prediction
    #
    dir_task_trainmodel_exp = os.path.join(workbase, 'task_trainmodel_exp')
    sci_obj.set_config_env(dict(config.items('CVParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # Subtask
    # train models to classify genes as active/inactive (TPM >= 1)
    dir_exp_statone = os.path.join(dir_task_trainmodel_exp, 'sub_status', 'active_geq1')

    cmd = config.get('Pipeline', 'trainmodel_expone_seq').replace('\n', ' ')
    trainmodel_expone_seq = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                           name='trainmodel_expone_seq',
                                           input=output_from(mrgtraindataexp_groups),
                                           filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                           output=os.path.join(dir_exp_statone,
                                                               '{subdir[0][0]}',
                                                               '{SAMPLE[0]}.rfcls.seq.uw.geq1.pck'),
                                           extras=[cmd, jobcall]).mkdir(dir_exp_statone)

    cmd = config.get('Pipeline', 'trainmodel_expone_sig').replace('\n', ' ')
    trainmodel_expone_sig = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                           name='trainmodel_expone_sig',
                                           input=output_from(mrgtraindataexp_groups),
                                           filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                           output=os.path.join(dir_exp_statone,
                                                               '{subdir[0][0]}',
                                                               '{SAMPLE[0]}.rfcls.sig.uw.geq1.pck'),
                                           extras=[cmd, jobcall]).mkdir(dir_exp_statone)

    cmd = config.get('Pipeline', 'trainmodel_expone_full').replace('\n', ' ')
    trainmodel_expone_full = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_expone_full',
                                            input=output_from(mrgtraindataexp_groups),
                                            filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                            output=os.path.join(dir_exp_statone,
                                                                '{subdir[0][0]}',
                                                                '{SAMPLE[0]}.rfcls.full.uw.geq1.pck'),
                                            extras=[cmd, jobcall]).mkdir(dir_exp_statone)

    # Subtask
    # train models to regress expression value for subset TPM >= 1
    dir_exp_valone = os.path.join(dir_task_trainmodel_exp, 'sub_value', 'subset_geq1')

    cmd = config.get('Pipeline', 'trainmodel_expvalone_seq').replace('\n', ' ')
    trainmodel_expvalone_seq = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='trainmodel_expvalone_seq',
                                              input=output_from(mrgtraindataexp_groups),
                                              filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                              output=os.path.join(dir_exp_valone,
                                                                  '{subdir[0][0]}',
                                                                  '{SAMPLE[0]}.rfreg.seq.uw.geq1.pck'),
                                              extras=[cmd, jobcall]).mkdir(dir_exp_valone)

    cmd = config.get('Pipeline', 'trainmodel_expvalone_sig').replace('\n', ' ')
    trainmodel_expvalone_sig = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='trainmodel_expvalone_sig',
                                              input=output_from(mrgtraindataexp_groups),
                                              filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                              output=os.path.join(dir_exp_valone,
                                                                  '{subdir[0][0]}',
                                                                  '{SAMPLE[0]}.rfreg.sig.uw.geq1.pck'),
                                              extras=[cmd, jobcall]).mkdir(dir_exp_valone)

    cmd = config.get('Pipeline', 'trainmodel_expvalone_full').replace('\n', ' ')
    trainmodel_expvalone_full = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                               name='trainmodel_expvalone_full',
                                               input=output_from(mrgtraindataexp_groups),
                                               filter=formatter('(?P<SAMPLE>[\w\.]+)\.feat\.h5'),
                                               output=os.path.join(dir_exp_valone,
                                                                   '{subdir[0][0]}',
                                                                   '{SAMPLE[0]}.rfreg.full.uw.geq1.pck'),
                                               extras=[cmd, jobcall]).mkdir(dir_exp_valone)

    task_trainmodel_exp = pipe.merge(task_func=touch_checkfile,
                                     name='task_trainmodel_exp',
                                     input=output_from(trainmodel_expone_seq, trainmodel_expone_sig, trainmodel_expone_full,
                                                       trainmodel_expvalone_seq, trainmodel_expvalone_sig, trainmodel_expvalone_full),
                                     output=os.path.join(dir_task_trainmodel_exp, 'task_trainmodel_exp.chk'))

    #
    # END: major task train model for gene expression prediction
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
    expfiles = collect_full_paths(dir_indata, '*/T0*.h5')
    groupinfos = os.path.join(os.path.dirname(config.get('Annotations', 'groupfile')), 'groupinfo_ro.json')
    matchings = os.path.join(os.path.dirname(config.get('Annotations', 'groupfile')), 'cellmatches_ro.json')
    cmd = config.get('Pipeline', 'testdataexp').replace('\n', ' ')
    testdataexp_groups = pipe.files(sci_obj.get_jobf('in_out'),
                                    annotate_test_datasets(mapfiles, roifiles, dir_sub_signal, expfiles,
                                                           groupinfos, dir_sub_cmpf_testdataexp,
                                                           cmd, jobcall),
                                    name='testdataexp_groups').mkdir(os.path.split(dir_sub_cmpf_testdataexp)[0])

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

    # ==============================
    # Major task: apply models to predict gene status (on/off) and expression level
    #

    # NB: this does not perform permutation testing at the moment (takes a loooong time)
    # iff that is enabled, change job config to CVParallelJobConfig!
    dir_task_applymodels_exp = os.path.join(workbase, 'task_applymodels_exp')
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    # ==================================================
    # annotations needed to categorize test output
    # matchings json from above
    # groupfile from above
    datasetfile = config.get('Annotations', 'dsetfile')
    groupings = prep_metadata(groupfile, 'gid')
    datasets = prep_metadata(datasetfile, 'id')
    cellmatches = js.load(open(matchings, 'r'))
    # ===================================================

    # subtask: predict gene status (on/off) for subset TPM >= 1
    dir_apply_statusone = os.path.join(dir_task_applymodels_exp, 'sub_status', 'active_geq1', '{groupid}')
    cmd = config.get('Pipeline', 'applyexpstatone').replace('\n', ' ')
    params_applyexpstatone = params_predict_testdata(dir_exp_statone, os.path.split(dir_mrg_test_datasets)[0],
                                                     dir_apply_statusone, cmd, jobcall)

    statone_out = os.path.join(dir_task_applymodels_exp, 'eval_models_stat-one.json')
    _ = annotate_test_output(params_applyexpstatone, datasets, groupings, cellmatches, 'Status >= 1', statone_out)

    applyexpstatone = pipe.files(sci_obj.get_jobf('inpair_out'),
                                 params_applyexpstatone,
                                 name='applyexpstatone')\
                                    .mkdir(os.path.split(dir_apply_statusone)[0])\
                                    .follows(mrgtestdataexp_groups).follows(task_trainmodel_exp)
    params_applyexpstatone = list()

    # subtask: predict gene expression level (using true TPM-based status)
    dir_apply_valueone = os.path.join(dir_task_applymodels_exp, 'sub_value', 'true_status', 'active_geq1', '{groupid}')
    cmd = config.get('Pipeline', 'applyexpvalone').replace('\n', ' ')
    params_applyexpvalone_true = params_predict_testdata(dir_exp_valone, os.path.split(dir_mrg_test_datasets)[0],
                                                         dir_apply_valueone, cmd, jobcall)

    valonetrue_out = os.path.join(dir_task_applymodels_exp, 'eval_models_val-one-true.json')
    _ = annotate_test_output(params_applyexpvalone_true, datasets, groupings, cellmatches, 'Value TPM >= 1', valonetrue_out)

    applyexpvalone = pipe.files(sci_obj.get_jobf('inpair_out'),
                                params_applyexpvalone_true,
                                name='applyexpvalone')\
                                    .mkdir(os.path.split(dir_apply_valueone)[0])\
                                    .follows(mrgtestdataexp_groups).follows(task_trainmodel_exp)
    params_applyexpvalone_true = list()

    # subtask: predict gene expression level (using predicted TPM-based status)
    dir_apply_valueonepred = os.path.join(dir_task_applymodels_exp, 'sub_value', 'pred_status', 'active_geq1', '{groupid}')
    cmd = config.get('Pipeline', 'applyexpvalonepred').replace('\n', ' ')
    load_metadata = os.path.split(dir_apply_statusone)[0]
    params_applyexpvalone_est = params_predict_testdata(dir_exp_valone, os.path.split(dir_mrg_test_datasets)[0],
                                                        dir_apply_valueonepred, cmd, jobcall, True, load_metadata)

    valoneest_out = os.path.join(dir_task_applymodels_exp, 'eval_models_val-one-est.json')
    _ = annotate_test_output(params_applyexpvalone_est, datasets, groupings, cellmatches, 'Value Est. >= 1', valoneest_out)

    applyexpvalonepred = pipe.files(sci_obj.get_jobf('inpair_out'),
                                    params_applyexpvalone_est,
                                    name='applyexpvalonepred')\
                                        .mkdir(os.path.split(dir_apply_valueonepred)[0])\
                                        .follows(mrgtestdataexp_groups).follows(task_trainmodel_exp)
    params_applyexpvalone_est = list()

    run_task_applymodels_exp = pipe.merge(task_func=touch_checkfile,
                                          name='task_applymodels_exp',
                                          input=output_from(applyexpstatone,
                                                            applyexpvalone,
                                                            applyexpvalonepred),
                                          output=os.path.join(dir_task_applymodels_exp, 'task_applymodels_exp.chk'))

    #
    # END: major task apply gene expression models
    # ==============================
    #
    # run_all = pipe.merge(task_func=touch_checkfile,
    #                      name='run_task_all',
    #                      input=output_from(run_task_sigcorr, run_task_sigmap,
    #                                        run_task_traindata_exp, run_task_testdata_exp,
    #                                        run_task_trainmodel_exp, run_task_applymodels_exp),
    #                      output=os.path.join(workbase, 'run_task_all.chk'))

    return pipe
