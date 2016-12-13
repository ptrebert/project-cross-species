# coding=utf-8

import os as os
import csv as csv

from ruffus import *

from pipelines.auxmods.auxiliary import collect_full_paths, touch_checkfile


def annotate_encode_files(fpaths, mddict, eiddict):
    """
    :param fpaths:
    :param mddict:
    :return:
    """
    annfiles = []
    narep = 1
    for fp in sorted(fpaths):
        facc = os.path.basename(fp).split('.')[0]
        # by construction, no missing keys
        infos = mddict[facc]
        k = infos['Assembly'], infos['Biosample term name'], infos['Biosample life stage'], infos['Lab']
        try:
            eid = eiddict[k]
        except KeyError:
            continue
        rep = infos['Biological replicate(s)']
        # manual fix
        if rep == '1, 2':  # not sure what this is, maybe 1 of 2? >> file ENCFF001KEH
            rep = 'n/a'
        if rep == 'n/a':
            rep = str(narep)
            narep += 1
        if isinstance(eid, list):
            # multi ID
            for i in eid:
                annfiles.append((fp, facc, rep, i, infos))
        else:
            annfiles.append((fp, facc, rep, eid, infos))
    assert annfiles, 'No ENCODE files annotated'
    return annfiles


def annotate_deep_files(fpaths, mddict, eiddict):
    """
    :param fpaths:
    :param mddict:
    :param eiddict:
    :return:
    """
    annfiles = []
    for fp in sorted(fpaths):
        fn = os.path.basename(fp)
        # by construction, no missing keys
        infos = mddict[fn]
        k = infos['assembly'], infos['biosample'], infos['lifestage'], infos['lab']
        try:
            eid = eiddict[k]
        except KeyError:
            continue
        rep = fn.split('.')[0].split('_')[6]
        _ = int(rep)  # just test, should always work
        if isinstance(eid, list):
            for i in eid:
                annfiles.append((fp, fn, rep, i, infos))
        else:
            annfiles.append((fp, fn, rep, eid, infos))
    assert annfiles, 'No DEEP files annotated'
    return annfiles


def link_encode_epigenomes(folder, mdfile, idfile, outfolder):
    """
    :param folder:
    :param mdfile:
    :param idfile:
    :param outfolder:
    :return:
    """
    os.makedirs(outfolder, exist_ok=True)
    bwfiles = collect_full_paths(folder, '*.bigWig', False)
    metadata = build_mddict(mdfile, 'File accession')
    eids = build_eiddict(idfile, 'ENCODE', 'epigenome')
    annfiles = annotate_encode_files(bwfiles, metadata, eids)
    linked_files = []
    for afp in annfiles:
        infos = afp[4]
        exptarget = infos['Experiment target']
        if exptarget == 'n/a':
            assert infos['Assay'] == 'DNase-seq', 'Unexpected assay type: {}'.format(infos['Assay'])
            exptarget = 'DNase'
        new_name = '_'.join([afp[3], infos['Assembly'], infos['Biosample term name'], exptarget, afp[2]])
        if afp[0].endswith('.bigWig'):
            fileext = '.bw'
        else:
            raise ValueError('Unexpected file type {}'.format(os.path.basename(afp[0])))
        new_path = os.path.join(outfolder, new_name + fileext)
        if not os.path.islink(new_path):
            os.link(afp[0], new_path)
        linked_files.append(new_path)
    assert linked_files, 'No ENCODE files linked to destination: {}'.format(outfolder)
    return linked_files


def link_deep_epigenomes(folder, mdfile, idfile, outfolder):
    """
    :param folder:
    :param mdfile:
    :param idfile:
    :param outfolder:
    :return:
    """
    os.makedirs(outfolder, exist_ok=True)
    bwfiles = collect_full_paths(folder, '*.bamcov.bw', False)
    metadata = build_mddict(mdfile, 'filename')
    eids = build_eiddict(idfile, 'DEEP', 'epigenome')
    annfiles = annotate_deep_files(bwfiles, metadata, eids)
    linked_files = []
    for afp in annfiles:
        infos = afp[4]
        exptarget = afp[1].split('_')[4]
        #new_name = '_'.join([afp[3], infos['assembly'], infos['biosample'], exptarget, afp[2]])
        new_name = '_'.join([afp[3], 'hs37d5', infos['biosample'], exptarget, afp[2]])
        fileext = '.bw'
        if not afp[0].endswith('.bw'):
            raise ValueError('Unexpected file type {}'.format(os.path.basename(afp[0])))
        new_path = os.path.join(outfolder, new_name + fileext)
        if not os.path.islink(new_path):
            os.link(afp[0], new_path)
        linked_files.append(new_path)
    assert linked_files, 'No DEEP files linked to destination: {}'.format(outfolder)
    return linked_files


def build_mddict(mdfile, key):
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


def build_eiddict(idfile, project, datatype):
    """
    :param idfile:
    :return:
    """
    eiddict = dict()
    with open(idfile, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        for row in reader:
            if row['project'] != project or row['type'] != datatype:
                continue
            k = row['assembly'], row['biosample'], row['lifestage'], row['lab']
            assert k not in eiddict, 'Duplicate key for IDs: {}'.format(k)
            if row['id'] == 'na':
                assert row['multid'] != 'na', 'Invalid ID for entry: {}'.format(row)
                eiddict[k] = row['multid'].split(',')
            else:
                eiddict[k] = row['id']
    assert eiddict, 'No metadata read from ID file {}'.format(idfile)
    return eiddict


def build_pipeline(args, config, sci_obj):
    """
    :param args:
    :param config:
    :param sci_obj:
    :return:
    """

    pipe = Pipeline(name=config.get('Pipeline', 'name'))

    # Main folders
    workbase = os.path.join(config.get('EnvPaths', 'workdir'), 'rawdata')

    # Special files
    enc_metadata = config.get('SpecialFiles', 'encmetadata')
    deep_metadata = config.get('SpecialFiles', 'deepmetadata')
    dataset_ids = config.get('SpecialFiles', 'datasetids')

    # Init task: link input epigenomes with readable names
    dlfolder = os.path.join(workbase, 'downloads')
    linkfolder = os.path.join(workbase, 'normlinks')

    encepi_init = pipe.originate(lambda x: x,
                                 link_encode_epigenomes(dlfolder,
                                                        enc_metadata,
                                                        dataset_ids,
                                                        linkfolder),
                                 name='encepi_init')

    deepepi_init = pipe.originate(lambda x: x,
                                  link_deep_epigenomes(dlfolder,
                                                       deep_metadata,
                                                       dataset_ids,
                                                       linkfolder))

    # ===========
    # Major task: convert all epigenomes from init tasks to HDF
    #
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        runjob = sci_obj.ruffus_gridjob()
    else:
        runjob = sci_obj.ruffus_localjob()

    bgout = os.path.join(workbase, 'conv', 'bedgraph')
    cmd = config.get('Pipeline', 'bwtobg')
    bwtobg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='bwtobg',
                            input=output_from(encepi_init, deepepi_init),
                            filter=suffix('.bw'),
                            output='.bg.gz',
                            output_dir=bgout,
                            extras=[cmd, runjob]).mkdir(bgout)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        runjob = sci_obj.ruffus_gridjob()
    else:
        runjob = sci_obj.ruffus_localjob()

    hdfout = os.path.join(workbase, 'conv', 'hdf')
    cmd = config.get('Pipeline', 'bgtohdfenc').replace('\n', ' ')
    re_filter = '(?P<EID>E[0-9]+)_(?P<ASSM>(hg19|mm9))_(?P<CELL>\w+)_(?P<MARK>\w+)_[0-9]+\.bg\.gz'
    bgtohdfenc = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                              name='bgtohdfenc',
                              input=output_from(bwtobg),
                              filter=formatter(re_filter),
                              output=os.path.join(hdfout, '{EID[0]}_{ASSM[0]}_{CELL[0]}_{MARK[0]}.h5'),
                              extras=[cmd, runjob]).mkdir(hdfout)

    cmd = config.get('Pipeline', 'bgtohdfdeep').replace('\n', ' ')
    re_filter = '(?P<EID>E[0-9]+)_(?P<ASSM>hs37d5)_(?P<CELL>\w+)_(?P<MARK>\w+)_[0-9]+\.bg\.gz'
    bgtohdfdeep = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                               name='bgtohdfdeep',
                               input=output_from(bwtobg),
                               filter=formatter(re_filter),
                               output=os.path.join(hdfout, '{EID[0]}_hg19_{CELL[0]}_{MARK[0]}.h5'),
                               extras=[cmd, runjob]).mkdir(hdfout)

    run_task_convepi = pipe.merge(task_func=touch_checkfile,
                                  name='task_convepi',
                                  input=output_from(bgtohdfenc, bgtohdfdeep),
                                  output=os.path.join(workbase, 'run_task_convepi.chk'))
    #
    # End of major task: convert epigenomes to HDF
    # ================


    return pipe

