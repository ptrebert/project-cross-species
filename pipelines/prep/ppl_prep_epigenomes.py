# coding=utf-8

import os as os
import csv as csv

from ruffus import *

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
            if row['project'] != 'ENCODE' or row['type'] != 'epigenome':
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
    metadata = build_encode_mddict(mdfile)
    eids = build_encode_eiddict(idfile)
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

    # Major task: convert all epigenomes from init tasks to HDF

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        runjob = sci_obj.ruffus_gridjob()
    else:
        runjob = sci_obj.ruffus_localjob()

    bgout = os.path.join(workbase, 'conv', 'bedgraph')
    cmd = config.get('Pipeline', 'bwtobg')
    bwtobg = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='bwtobg',
                            input=output_from(encepi_init),
                            filter=suffix('.bw'),
                            output='.bg.gz',
                            output_dir=bgout,
                            extras=[cmd, runjob]).mkdir(bgout)

    return pipe

