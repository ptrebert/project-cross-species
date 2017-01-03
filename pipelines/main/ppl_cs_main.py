# coding=utf-8

import os as os
import re as re
import json as js
import csv as csv
import collections as col
import fnmatch as fnm
import datetime as dt

from ruffus import *

from pipelines.auxmods.auxiliary import collect_full_paths, touch_checkfile, dbg_param_list


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
    return p1, p2, hist, dnase


def make_groups_compatible(groupings, epibundles):
    """
    :param groupings:
    :param epibundles:
    :return:
    """
    groups = dict()
    cmatches = col.defaultdict(set)
    with open(groupings, 'r', newline='') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        for r in rows:
            p1, p2, h, d = select_compatible_libs(epibundles[r['partner1']], epibundles[r['partner2']])
            p1_params = ' '.join([x['param'] for x in p1])
            p2_params = ' '.join([x['param'] for x in p2])
            gid = r['gid']
            extgid = gid + str(h) + str(d)
            entry = {'partner1': p1, 'partner2': p2, 'params1': p1_params, 'params2': p2_params, 'extgid': extgid}
            groups[gid] = entry
            c1, c2 = r['comment'].split('-')
            cmatches[c1].add(c2)
            cmatches[c2].add(c1)
    assert groups, 'No group dictionary constructed'
    return groups, cmatches
            

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
    groupings, matchings = make_groups_compatible(groupfile, epibundles)
    arglist = []
    uniq = set()
    for mapf in mapfiles:
        target, query = mapf['target'], mapf['query']
        for gid, partners in groupings.items():
            assm1, assm2 = partners['partner1'][0]['assembly'], partners['partner2'][0]['assembly']
            compat_roi = [roi for roi in roifiles if roi['assembly'] == assm1]
            mapref = 'target' if assm1 == target else 'query'
            for roi in compat_roi:
                outfolder = outraw.format(**{'target': target, 'query': query})
                outfile = '_'.join([partners['extgid'], partners['partner1'][0]['EID'],
                                    assm1, partners['partner1'][0]['cell']])
                outfile += '.' + roi['regtype'] + '.h5'
                outpath = os.path.join(outfolder, outfile)
                assert outpath not in uniq, 'Created duplicate: {}'.format(outpath)
                uniq.add(outpath)
                params = {'target': target, 'query': query, 'mapref': mapref, 'genome': assm1,
                          'regtype': roi['regtype'], 'mapfile': mapf['path'],
                          'datafiles': partners['params1'], 'info': partners['extgid']}
                tmp = cmd.format(**params)
                arglist.append([roi['path'], outpath, tmp, jobcall])
            if assm2 == assm1:
                continue
            compat_roi = [roi for roi in roifiles if roi['assembly'] == assm2]
            mapref = 'target' if assm2 == target else 'query'
            for roi in compat_roi:
                outfolder = outraw.format(**{'target': target, 'query': query})
                outfile = '_'.join([partners['extgid'], partners['partner2'][0]['EID'],
                                    assm2, partners['partner2'][0]['cell']])
                outfile += '.' + roi['regtype'] + '.h5'
                outpath = os.path.join(outfolder, outfile)
                assert outpath not in uniq, 'Created duplicate: {}'.format(outpath)
                uniq.add(outpath)
                params = {'target': target, 'query': query, 'mapref': mapref, 'genome': assm2,
                          'regtype': roi['regtype'], 'mapfile': mapf['path'],
                          'datafiles': partners['params2'], 'info': partners['extgid']}
                tmp = cmd.format(**params)
                arglist.append([roi['path'], outpath, tmp, jobcall])
    # store group info
    store_group_info(groupfile, groupings)
    # store cell type matches
    store_cell_matches(groupfile, matchings)
    if mapfiles:
        assert arglist, 'No calls created: annotate_training_groups'
    return arglist


def store_group_info(groupfile, grouping):
    """
    :param groupfile:
    :param grouping:
    :return:
    """
    timestr = dt.datetime.now().strftime('%A_%Y-%m-%d_%H:%M:%S')
    filepath = os.path.join(os.path.dirname(groupfile), 'groupinfo_ro.json')
    if not os.path.isfile(filepath):
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'groupinfo': {}}, dump, indent=1, sort_keys=True)
    rewrite = False
    with open(filepath, 'r') as dumped:
        gi = js.load(dumped)['groupinfo']
        if gi != grouping:
            rewrite = True
    if rewrite:
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'groupinfo': grouping}, dump, indent=1, sort_keys=True)
    return


def store_cell_matches(groupfile, matchings):
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
            js.dump({'timestamp': timestr, 'matchings': {}}, dump, indent=1, sort_keys=True)
    rewrite = False
    with open(filepath, 'r') as dumped:
        cmi = js.load(dumped)['matchings']
        if cmi != serialize:
            rewrite = True
    if rewrite:
        with open(filepath, 'w') as dump:
            js.dump({'timestamp': timestr, 'matchings': serialize}, dump, indent=1, sort_keys=True)
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
    featfiles = collect_full_paths(featdir, '*/G[0-9]*.h5', True)
    expfiles = collect_full_paths(expdir, '*/T0*.h5')
    collector = col.defaultdict(list)
    for fp in featfiles:
        if mapped:
            parts = os.path.basename(fp).split('.')
            k, i = '.'.join(parts[:3]), parts[3]
        else:
            k, i = os.path.basename(fp).split('.')[:2]
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
            pat = group.split('_', 2)[-1]
            match_exp = fnm.filter(expfiles, '*/T0*_{}_mRNA*'.format(pat))
        try:
            assert match_exp, 'No matched expression files for group {} / pattern {} / {}'.format(group, pat, files)
        except AssertionError as ae:
            if '_ESE14' in pat:
                continue
            else:
                raise ae
        mergefiles = ' '.join([t[1] for t in files])
        if mapped:
            gid, eid, assm, cell = group.split('.')[0].split('_')
        else:
            gid, eid, assm, cell = group.split('_')
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
                og = os.path.join(qry, trg, gid, 'features')
            else:
                og = os.path.join(trg, qry, gid, 'features')
            params = {'outgroup': og, 'mergefiles': mergefiles}
            outbase = outraw.format(**{'target': trg, 'query': qry})
            if mapped:
                outfile = '_'.join([gid, eid, tid, assm, cell]) + '.from.{}.feat.h5'.format(trg)
            else:
                outfile = '_'.join([gid, eid, tid, assm, cell]) + '.feat.h5'
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
                                params = {'target': target, 'query': query, 'mapref': 'query', 'genome': query,
                                          'regtype': roi['regtype'], 'mapfile': mapf['path'],
                                          'datafiles': datafiles, 'info': extgid}
                                tmp = cmd.format(**params)
                                arglist.append([roi['path'], outpath, tmp, jobcall])
    if mapfiles:
        assert arglist, 'No argument list for building test datasets created'
    return arglist


# ===============
# everything below needs to be revised
# ===============


def make_corr_pairs(inputfiles, targets, mapdata, roifiles, outbase, cmd, jobcall):
    """
    :param inputfiles:
    :param targets:
    :param mapdata:
    :param roifiles:
    :param outbase:
    :param cmd:
    :param jobcall:
    :return:
    """
    mates = []
    for root, dirs, files in os.walk(mapdata):
        if files:
            basedir, subdir = os.path.split(root)
            query, from_assm, target = subdir.split('_')
            if query in targets and target in targets:
                # there is assayed signal data just for some species/assemblies
                for f in files:
                    mates.append(os.path.join(root, f))
    params = []
    # assayed: ENCSR000CGJ_mm9_CH12_H3K27ac_BRUCSD.srcsig.h5
    # mapped: ENCSR000AMO_mm9_CH12_H3K27ac_BBBRD.from.hg19.HepG2.mapsig.h5
    done = set()
    for inf in inputfiles:
        if not inf.endswith('.srcsig.h5'):
            continue
        _, qry, qcell, lib, qlab = inf.split('.')[0].split('_')
        matched_files = fnm.filter(mates, '*_' + qry + '_*' + '_' + lib + '_*')
        assert matched_files, 'No matches for input file {}'.format(inf)
        for mf in matched_files:
            mfname = os.path.basename(mf)
            _, _, trg, tcell, _, _ = mfname.split('.')
            if (os.path.basename(inf), trg, tcell) in done:
                continue
            for roi in roifiles:
                roilabel, roifile = roi.split(':')
                this_roi = roifile.format(**{'assembly': qry})
                tmp = cmd.format(**{'query': qry, 'target': trg, 'roifile': this_roi})
                outname = 'corr_{}_{}_to_{}.json'.format(roilabel,
                                                         os.path.basename(inf).strip('.h5'),
                                                         os.path.basename(mfname).strip('.h5'))
                outpath = os.path.join(outbase, '{}_from_{}'.format(qry, trg), outname)
                params.append([[inf, mf], outpath, tmp, jobcall])
            done.add((os.path.basename(inf), trg, tcell))
    # assert params, 'No parameter pairs created for signal/mapping correlation'
    return params


def make_srcsig_pairs(inputfiles, targets, roifiles, outbase, cmd, jobcall):
    """
    :param inputfiles:
    :param targets:
    :param roifiles:
    :param outbase:
    :param cmd:
    :param jobcall:
    :return:
    """
    params = []
    # assayed: ENCSR000CGJ_mm9_CH12_H3K27ac_BRUCSD.srcsig.h5
    done = set()
    for inf in inputfiles:
        if not inf.endswith('.srcsig.h5'):
            continue
        _, trg, cell, lib, lab = inf.split('.')[0].split('_')
        # the following assumes two targets: hg19, mm9
        assert len(targets) == 2, 'More than two targets: {}'.format(targets)
        if targets[0] == trg:
            qry = targets[1]
        else:
            qry = targets[0]
        matched_files = fnm.filter(inputfiles, '*_' + trg + '_*' + '_' + lib + '_*.srcsig.h5')
        assert matched_files, 'No matches for input file {}'.format(inf)
        for mf in matched_files:
            mfname = os.path.basename(mf)
            if mf == inf or (inf, mf) in done or (mf, inf) in done:
                continue
            for roi in roifiles:
                roilabel, roifile = roi.split(':')
                this_roi = roifile.format(**{'assembly': trg})
                tmp = cmd.format(**{'query': qry, 'target': trg, 'roifile': this_roi})
                outname = 'corr_{}_{}_to_{}.json'.format(roilabel,
                                                         os.path.basename(inf).strip('.h5'),
                                                         os.path.basename(mfname).strip('.h5'))
                outpath = os.path.join(outbase, '{}_from_{}'.format(trg, qry), outname)
                params.append([[inf, mf], outpath, tmp, jobcall])
            done.add((inf, mf))
    # assert params, 'No parameter pairs created for signal/mapping correlation'
    return params


def annotate_histone_files(fpaths):
    """
    :param fpaths:
    :return:
    """
    ret = []
    for filep in fpaths:
        fp, fn = os.path.split(filep)
        expid_parts = fn.split('.')[0].split('_')
        run_parts = fn.split('.')
        assm, cell, mark, lab = expid_parts[1:]
        if run_parts[-2] == 'mapsig':
            source = run_parts[-3]
            if source in CTYPE_MATCH:
                ret.append((cell, source, lab, mark, filep))
            else:
                ret.append((cell, '', lab, mark, filep))
        else:
            ret.append((cell, '', lab, mark, filep))
    assert ret, 'Annotating data files failed for {}'.format(fpaths)
    return ret


def make_predtestdata_params(modelroot, dataroot, outroot, cmd, jobcall, addmd=False):
    """
    :param modelroot:
    :param dataroot:
    :param cmd:
    :param jobcall:
    :return:
    """
    all_models = collect_full_paths(modelroot, '*.pck')
    all_testdata = collect_full_paths(dataroot, '*genes*.h5')

    arglist = []
    done = set()
    for testdata in all_testdata:
        testname = os.path.basename(testdata)
        parts = testname.split('.')
        from_cell = parts[-2]
        from_assm = parts[-3]
        to_assm = parts[0].split('_')[1]
        models = fnm.filter(all_models, '*_{}_{}_*.to.{}*pck'.format(from_assm, from_cell, to_assm))
        if not models:
            continue
        for mdf in models:
            temp_out = outroot.format(**{'from': from_assm, 'to': to_assm})
            modelname = os.path.basename(mdf)
            if (testname, modelname) in done:
                continue
            outfile = testname.replace('.h5', '') + '_with_' + modelname.replace('.pck', '') + '.json'
            outpath = os.path.join(temp_out, outfile)
            if addmd:
                mdfile = outpath.replace('rfreg', 'rfcls').replace('sub_value', 'sub_status').replace('pred_status/', '')
                tmp = cmd.format(**{'mdfile': mdfile})
                arglist.append([[testdata, mdf], outpath, tmp, jobcall])
            else:
                arglist.append([[testdata, mdf], outpath, cmd, jobcall])
            done.add((testname, modelname))
    return arglist


def make_peak_testdata(peakdir, peakext, signalroot, signalext, outbasedir, cmd, jobcall):

    peakfiles = os.listdir(peakdir)
    peakfiles = fnm.filter(peakfiles, peakext)

    all_sigfiles = []
    for root, dirs, files in os.walk(signalroot):
        if files:
            sigfiles = fnm.filter(files, signalext)
            sigfiles = [os.path.join(root, sf) for sf in sigfiles]
            all_sigfiles.extend(sigfiles)

    arglist = []
    done_pairs = set()
    for pf in peakfiles:
        parts = pf.split('.')[0].split('_')
        base_filter = '*' + parts[1] + '_'
        ctype_match = CTYPE_MATCH[parts[2]]
        for ctm in ctype_match:
            match_filter = base_filter + ctm + '_' + parts[3] + '*'
            matched_sigfiles = fnm.filter(all_sigfiles, match_filter)
            for msig in matched_sigfiles:
                if (pf, msig) in done_pairs:
                    continue
                sigparts = os.path.basename(msig).split('.')
                assert sigparts[1] == 'trg', 'Not a mapped/estimated signal file'
                trgassembly = sigparts[2]
                sigparts = sigparts[0].split('_')
                subdir = outbasedir.format(**{'assembly': parts[1], 'query': trgassembly})
                outname = '_'.join([parts[0], parts[1], parts[2] + '-' + sigparts[2], parts[3], parts[4] + '-' + sigparts[4]])
                outname += '.trg.' + trgassembly + '.testdat.pk.h5'
                outpath = os.path.join(subdir, outname)
                mycmd = cmd.format(**{'signalfile': msig, 'mark': parts[3], 'assembly': parts[1], 'query': trgassembly})
                arglist.append([os.path.join(peakdir, pf), outpath, mycmd, jobcall])
    return arglist


def make_compfeat_genes(expdir, regdir, sigdir, outdir, cmd, jobcall):
    """
    :return:
    """
    expression_files = collect_full_paths(expdir, '*.h5')
    gene_regions = collect_full_paths(regdir, '*.h5')
    signalfiles = collect_full_paths(sigdir, '*.h5')

    arglist = []
    for expf in expression_files:
        # ENCSR000CPE_hg19_HepG2_mRNA_L13.genc-v19.h5
        _, assembly, cell, _, _ = expf.split('_')
        gencode_version = expf.split('.')[1].split('-')[1]
        matched_cells = CTYPE_MATCH[cell]
        gene_files = fnm.filter(gene_regions, '*' + gencode_version + '_' + assembly + '*.h5')
        assert len(gene_files) == 2, 'Missing gene regions for expression file: {}'.format(expf)
        for genf in gene_files:
            _, _, regtype, _ = genf.split('.')
            msigfiles = filter(lambda fp: '_' + assembly + '_' in fp and any(['_' + ct + '_' in fp for ct in matched_cells]), signalfiles)
            for msig in msigfiles:
                _, _, mcell, mark, labext = os.path.basename(msig).split('_')
                if mark in ['DNaseI', 'InpControl']:
                    continue
                lab, _, target, _, _ = labext.split('.')
                outgroup = os.path.join('gencode', gencode_version, regtype, mark)
                if regtype == 'body':
                    win = '1000'
                    step = '500'
                elif regtype == 'promoter':
                    win = '500'
                    step = '250'
                else:
                    raise ValueError('Unsupported region type for windowing: {}'.format(regtype))
                tmp = cmd.format(**{'assembly': assembly, 'query': target,
                                    'outputgroup': outgroup, 'mark': mark, 'window': win, 'stepsize': step})
                outfile = os.path.join(outdir, '_'.join(['gencode', gencode_version, regtype, assembly, mark, cell + '-' + mcell, lab]))
                outfile += '.'.join(['', 'trg', target, 'scnfeat', 'h5'])
                args = [[genf, msig], outfile, tmp, jobcall]
                arglist.append(args)
    return arglist


def make_model_featdata_pairs(modeldir, datadir, outdir, cmd, jobcall):
    """
    :param modeldir:
    :param datadir:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    arglist = []
    all_models = []
    for root, dirs, files in os.walk(modeldir):
        if files:
            files = fnm.filter(files, '*.full.pck')
            for f in files:
                all_models.append(os.path.join(root, f))
    all_featdata = os.listdir(datadir)
    for fp in all_featdata:
        fn, _, target, _, _ = fp.split('.')
        _, gencver, regtype, assembly, mark, cells, lab = fn.split('_')
        if 'Control' in mark:
            continue
        cell, matched_cell = cells.split('-')
        model_select = '*' + '_'.join([target, matched_cell, mark, lab]) + '*' + assembly + '*'
        use_models = fnm.filter(all_models, model_select)
        assert len(use_models) > 0, 'No model files selected for file / pattern: {} / {}'.format(fp, model_select)
        # for debug purposes only
        assert len(use_models) == 1, 'Too many models for the same prediction task: {}'.format(use_models)
        # end debug
        outfilename = fp.replace('scnfeat', 'pkpred')
        outfilepath = os.path.join(outdir, outfilename)
        ingroup = os.path.join('gencode', gencver, regtype, mark)
        outgroup = ingroup
        infilepath = os.path.join(datadir, fp)
        tmp = cmd.format(**{'inputfile': infilepath, 'outputfile': outfilepath, 'inputgroup': ingroup,
                            'outputgroup': outgroup, 'modelfile': use_models[0]})
        arglist.append([infilepath, outfilepath, tmp, jobcall])
    return arglist


def build_histone_feature_map(histfiles, assayed=True, defaultext='.5hm.'):
    """
    :param histfiles:
    :return:
    """
    use_marks = ['H3K4me1', 'H3K4me3', 'H3K27me3', 'H3K36me3', 'H3K9me3']
    get_mark = re.compile('_(?P<MARK>H3K[0-9a-z]+)_')
    get_lab = re.compile('_(?P<LAB>L[0-9]{2})')
    collect = col.defaultdict(list)
    for hf in histfiles:
        mobj = get_mark.search(hf)
        if mobj is None or mobj.group('MARK') not in use_marks:
            continue
        mark = mobj.group('MARK')
        hfbase = os.path.basename(hf)
        if assayed:
            assembly = hfbase.split('_')[1]
            cell = hfbase.split('_')[2]
        else:
            assembly = hfbase.split('_')[3]
            cell = hfbase.split('_')[5]
        lab = get_lab.search(hf).group('LAB')
        collect[(assembly, cell)].append((lab, mark, hf))
    file_map = col.defaultdict(list)
    for key, files in collect.items():
        for lab in sorted(set([info[0] for info in files])):
            # default ext: 5hm = 5 histone marks / peaks
            # 5ht = 5 histone tracks / signal
            ext = defaultext
            lcount = col.Counter()
            filelist = ''
            found_marks = set()
            for l, m, fp in files:
                if l != lab:
                    continue
                if m in found_marks:
                    continue
                filelist += ' ' + ':'.join([m, '', fp])
                found_marks.add(m)
                lcount[l] += 1
            if len(found_marks) < len(use_marks):
                for l, m, fp in files:
                    if l == lab:
                        continue
                    if m in found_marks:
                        continue
                    filelist += ' ' + ':'.join([m, '', fp])
                    found_marks.add(m)
                    lcount[l] += 1
                    if len(found_marks) == len(use_marks):
                        break
            for l in sorted(lcount.keys()):
                ext += str(lcount[l]) + 'x' + l + '.'
            file_map[key].append((ext, filelist.strip()))
    return file_map


def make_traindata_geneparts_params(genedir, featmap, outdir, cmd, jobcall, targets=''):
    """
    :param genedir:
    :param featmap:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    genefiles = os.listdir(genedir)
    genefiles = fnm.filter(genefiles, '*.h5')
    params = []
    for gf in genefiles:
        _, gencver, assembly = gf.split('.')[0].split('_')
        regtype = gf.split('.')[2]
        for key, all_featfiles in featmap.items():
            if key[0] != assembly:
                continue
            cell = key[1]
            for featfile in all_featfiles:
                outfile = '_'.join(['gencode', gencver, assembly, cell])
                outfile += '.' + regtype + featfile[0] + 'h5'
                outpath = os.path.join(outdir, outfile)
                inpath = os.path.join(genedir, gf)
                outgroup = '/'.join(['gencode', gencver, regtype])
                if targets:
                    trg_assemblies = targets.split()
                    query = [t for t in trg_assemblies if t != assembly]
                    assert len(query) == 1, 'Multiple queries for target: {} - {}'.format(query, assembly)
                    query = query[0]
                    tmp = cmd.format(**{'featfiles': featfile[1], 'outputgroup': outgroup,
                                        'assembly': assembly, 'query': query})
                else:
                    tmp = cmd.format(**{'featfiles': featfile[1], 'outputgroup': outgroup, 'assembly': assembly})
                params.append([inpath, outpath, tmp, jobcall])
    return params


def make_merge_geneparts_params(geneparts, expression, outdir, cmd, jobcall, assayed=True):
    """
    :param geneparts:
    :param expression:
    :param outdir:
    :param cmd:
    :param jobcall:
    :return:
    """
    arglist = []
    done = set()
    for fp in geneparts:
        if fp in done:
            continue
        if '.body.' in fp:
            infile1 = fp.replace('.body.', '.promoter.')
            infile2 = fp
        else:
            infile1 = fp
            infile2 = fp.replace('.promoter.', '.body.')
        done.add(infile1)
        done.add(infile2)
        cell = fp.split('.')[0].split('_')[3]
        if not assayed:
            matched_cell = cell
            if len(cell.split('-')) == 1:
                # TODO for mapped sig tracks, should use the same naming: cell-mapped_cell
                cell = CTYPE_MATCH[matched_cell][0]
                matched_cell = cell + '-' + matched_cell
            else:
                cell = cell.split('-')[0]
        # for ESE14, no expression data available
        if cell == 'ESE14':
            continue
        expfiles = fnm.filter(expression, '*' + cell + '*.h5')
        assert expfiles, 'No expression files for cell {}'.format(cell)
        for ef in expfiles:
            sample, gencver, _ = os.path.basename(ef).split('.')
            if not assayed:
                sample = sample.replace(cell, matched_cell)
            labext = os.path.basename(infile1).split('.')[2:]
            outfile = '.'.join([sample, gencver] + labext)
            outpath = os.path.join(outdir, outfile)
            tmp = cmd.format(**{'valuefile': ef, 'outgroup': os.path.join('gencode', gencver.split('-')[1], 'gene')})
            arglist.append([[infile1, infile2], outpath, tmp, jobcall])
    return arglist


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
    group_partners = use_targets
    mapfiles = collect_mapfiles(dir_maps, group_partners, group_partners)
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
                                        name='mrgtraindataexp_groups').mkdir(os.path.split(dir_mrg_train_datasets)[0])

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
                                           filter=formatter('(?P<SAMPLE>\w+)\.feat\.h5'),
                                           output=os.path.join(dir_exp_statone,
                                                               '{subdir[0][0]}',
                                                               '{SAMPLE[0]}.rfcls.seq.uw.geq1.pck'),
                                           extras=[cmd, jobcall]).mkdir(dir_exp_statone)

    cmd = config.get('Pipeline', 'trainmodel_expone_sig').replace('\n', ' ')
    trainmodel_expone_sig = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                           name='trainmodel_expone_sig',
                                           input=output_from(mrgtraindataexp_groups),
                                           filter=formatter('(?P<SAMPLE>\w+)\.feat\.h5'),
                                           output=os.path.join(dir_exp_statone,
                                                               '{subdir[0][0]}',
                                                               '{SAMPLE[0]}.rfcls.sig.uw.geq1.pck'),
                                           extras=[cmd, jobcall]).mkdir(dir_exp_statone)

    cmd = config.get('Pipeline', 'trainmodel_expone_full').replace('\n', ' ')
    trainmodel_expone_full = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                            name='trainmodel_expone_full',
                                            input=output_from(mrgtraindataexp_groups),
                                            filter=formatter('(?P<SAMPLE>\w+)\.feat\.h5'),
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
                                              filter=formatter('(?P<SAMPLE>\w+)\.feat\.h5'),
                                              output=os.path.join(dir_exp_valone,
                                                                  '{subdir[0][0]}',
                                                                  '{SAMPLE[0]}.rfreg.seq.uw.geq1.pck'),
                                              extras=[cmd, jobcall]).mkdir(dir_exp_valone)

    cmd = config.get('Pipeline', 'trainmodel_expvalone_sig').replace('\n', ' ')
    trainmodel_expvalone_sig = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                              name='trainmodel_expvalone_sig',
                                              input=output_from(mrgtraindataexp_groups),
                                              filter=formatter('(?P<SAMPLE>\w+)\.feat\.h5'),
                                              output=os.path.join(dir_exp_valone,
                                                                  '{subdir[0][0]}',
                                                                  '{SAMPLE[0]}.rfreg.sig.uw.geq1.pck'),
                                              extras=[cmd, jobcall]).mkdir(dir_exp_valone)

    cmd = config.get('Pipeline', 'trainmodel_expvalone_full').replace('\n', ' ')
    trainmodel_expvalone_full = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                               name='trainmodel_expvalone_full',
                                               input=output_from(mrgtraindataexp_groups),
                                               filter=formatter('(?P<SAMPLE>\w+)\.feat\.h5'),
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
                                       name='mrgtestdataexp_groups').mkdir(os.path.split(dir_mrg_test_datasets)[0])

    task_testdata_exp = pipe.merge(task_func=touch_checkfile,
                                   name='task_testdata_exp',
                                   input=output_from(testdataexp_groups, mrgtestdataexp_groups),
                                   output=os.path.join(dir_task_traindata_exp, 'run_task_testdata_exp.chk'))

    #
    # END: major task generate test data
    # ==============================
    #
    # # ==============================
    # # Major task: apply models to predict gene status (on/off) and expression level
    # #
    #
    # # NB: this does not perform permutation testing at the moment (takes a loooong time)
    # # iff that is enabled, change job config to CVParallelJobConfig!
    # dir_task_applymodels_exp = os.path.join(workbase, 'task_applymodels_exp')
    # sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('EnvConfig')))
    # if args.gridmode:
    #     jobcall = sci_obj.ruffus_gridjob()
    # else:
    #     jobcall = sci_obj.ruffus_localjob()
    #
    # # subtask: predict gene status (on/off) for subsets TPM > 0 OR TPM >= 1
    # dir_apply_statuszero = os.path.join(dir_task_applymodels_exp, 'sub_status', 'active_grt0')
    # cmd = config.get('Pipeline', 'applyexpstatzero').replace('\n', ' ')
    # applyexpstatzero = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                               make_predtestdata_params(dir_exp_statzero,
    #                                                        os.path.split(dir_base_testdataexp)[0],
    #                                                        dir_apply_statuszero, cmd, jobcall),
    #                               name='applyexpstatzero').mkdir(dir_apply_statuszero).follows(run_task_testdata_exp)
    #
    # dir_apply_statusone = os.path.join(dir_task_applymodels_exp, 'sub_status', 'active_geq1')
    # cmd = config.get('Pipeline', 'applyexpstatone').replace('\n', ' ')
    # applyexpstatone = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                              make_predtestdata_params(dir_exp_statone,
    #                                                       os.path.split(dir_base_testdataexp)[0],
    #                                                       dir_apply_statusone, cmd, jobcall),
    #                              name='applyexpstatone').mkdir(dir_apply_statusone).follows(run_task_testdata_exp)
    #
    # # subtask: predict gene expression level (using true TPM-based status)
    # dir_apply_valuezero = os.path.join(dir_task_applymodels_exp, 'sub_value', 'true_status', 'active_grt0')
    # cmd = config.get('Pipeline', 'applyexpvalzero').replace('\n', ' ')
    # applyexpvalzero = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                              make_predtestdata_params(dir_exp_valzero,
    #                                                       os.path.split(dir_base_testdataexp)[0],
    #                                                       dir_apply_valuezero, cmd, jobcall),
    #                              name='applyexpvalzero').mkdir(dir_apply_valuezero).follows(run_task_testdata_exp)
    #
    # dir_apply_valueone = os.path.join(dir_task_applymodels_exp, 'sub_value', 'true_status', 'active_geq1')
    # cmd = config.get('Pipeline', 'applyexpvalone').replace('\n', ' ')
    # applyexpvalone = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                             make_predtestdata_params(dir_exp_valone,
    #                                                      os.path.split(dir_base_testdataexp)[0],
    #                                                      dir_apply_valueone, cmd, jobcall),
    #                             name='applyexpvalone').mkdir(dir_apply_valueone).follows(run_task_testdata_exp)
    #
    # # subtask: predict gene expression level (using predicted TPM-based status)
    # dir_apply_valuezeropred = os.path.join(dir_task_applymodels_exp, 'sub_value', 'pred_status', 'active_grt0')
    # cmd = config.get('Pipeline', 'applyexpvalzeropred').replace('\n', ' ')
    # applyexpvalzeropred = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                                  make_predtestdata_params(dir_exp_valzero,
    #                                                           os.path.split(dir_base_testdataexp)[0],
    #                                                           dir_apply_valuezeropred, cmd, jobcall, True),
    #                                  name='applyexpvalzeropred').mkdir(dir_apply_valuezeropred).follows(run_task_testdata_exp)
    #
    # dir_apply_valueonepred = os.path.join(dir_task_applymodels_exp, 'sub_value', 'pred_status', 'active_geq1')
    # cmd = config.get('Pipeline', 'applyexpvalonepred').replace('\n', ' ')
    # applyexpvalonepred = pipe.files(sci_obj.get_jobf('inpair_out'),
    #                                 make_predtestdata_params(dir_exp_valone,
    #                                                          os.path.split(dir_base_testdataexp)[0],
    #                                                          dir_apply_valueonepred, cmd, jobcall, True),
    #                                 name='applyexpvalonepred').mkdir(dir_apply_valueonepred).follows(run_task_testdata_exp)
    #
    # run_task_applymodels_exp = pipe.merge(task_func=touch_checkfile,
    #                                       name='run_task_applymodels_exp',
    #                                       input=output_from(applyexpstatzero, applyexpstatone,
    #                                                         applyexpvalzero, applyexpvalone,
    #                                                         applyexpvalzeropred, applyexpvalonepred),
    #                                       output=os.path.join(dir_task_applymodels_exp, 'task_applymodels_exp.chk'))
    #
    # #
    # # END: major task apply gene expression models
    # # ==============================
    #
    # run_all = pipe.merge(task_func=touch_checkfile,
    #                      name='run_task_all',
    #                      input=output_from(run_task_sigcorr, run_task_sigmap,
    #                                        run_task_traindata_exp, run_task_testdata_exp,
    #                                        run_task_trainmodel_exp, run_task_applymodels_exp),
    #                      output=os.path.join(workbase, 'run_task_all.chk'))

    return pipe
