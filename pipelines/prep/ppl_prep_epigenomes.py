# coding=utf-8

import os as os
import re as re
import csv as csv
import collections as col

from ruffus import *

from pipelines.auxmods.auxiliary import collect_full_paths, \
    touch_checkfile, prep_dset_ids, prep_metadata


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
        rep = fn.split('.')[0].split('_')[1][3]
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
    metadata = prep_metadata(mdfile, 'File accession')
    eids = prep_dset_ids(idfile, 'ENCODE', 'epigenome')
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
    metadata = prep_metadata(mdfile, 'filename')
    eids = prep_dset_ids(idfile, 'DEEP', 'epigenome')
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


def make_norm_calls(inputfiles, outdir, components, cmd, jobcall):
    """
    :param inputfiles:
    :param outdir:
    :param components:
    :param cmd:
    :param jobcall:
    :return:
    """
    match_filename = re.compile(components)
    collect_input = col.defaultdict(list)
    collect_output = col.defaultdict(list)
    for fp in inputfiles:
        filename = os.path.basename(fp)
        mobj = match_filename.match(filename)
        if mobj is None:
            continue
        outputfile = os.path.join(outdir, filename)
        assm, mark = mobj.group('ASSM'), mobj.group('MARK')
        collect_input[(assm, mark)].append(fp)
        collect_output[(assm, mark)].append(outputfile)
    arglist = []
    expected_output = []
    for (assm, mark), infiles in collect_input.items():
        tmp = cmd.format(**{'assm': assm, 'inputfiles': ' '.join(infiles)})
        outfiles = collect_output[(assm, mark)]
        assert len(outfiles) == len(infiles), 'Sample mismatch: {} vs {}'.format(infiles, outfiles)
        arglist.append([tmp, jobcall])
        expected_output.extend([os.path.basename(f) for f in outfiles])
    return arglist, expected_output


def parse_bp_sra_acc(acc_file, url_file, outbase, cmd, jobcall):
    """
    :param acc_file:
    :param url_file:
    :param outbase:
    :param cmd:
    :param jobcall:
    :return:
    """
    accessions = []
    with open(acc_file, 'r') as accfile:
        for line in accfile:
            if not line.strip():
                continue
            accessions.append(line.strip())
    args = []
    with open(url_file, 'r', newline='') as ena_file:
        rows = csv.DictReader(ena_file, delimiter='\t')
        for entry in rows:
            acc = entry['run_accession']
            if acc in accessions:
                for rn, path in enumerate(entry['fastq_ftp'].split(';'), start=1):
                    tmp = cmd.format(**{'outputdir': outbase, 'ftp_url': path,
                                        'repnum': str(rn), 'sra_acc': acc})
                    outfile = os.path.join(outbase, acc + '_r' + str(rn) + '_md.csv')
                    args.append([acc_file, outfile, tmp, jobcall])
    assert args, 'No arguments created for BP download'
    return args


def create_blueprint_annotation(metadata_files, outputfile):
    """
    :param metadata_files:
    :return:
    """
    # this annotation is based on manual inspection of the annotation
    # in EBI/ENA
    # examples: SAMEA3928330 and SAMEA3855949
    tcell = {'full_biosample': 'naive_CD4p_Tcells', 'short_biosample': 'CD4pTN'}
    bcell = {'full_biosample': 'mature_resting_Bcells', 'short_biosample': 'Bcells'}
    collection = []
    delete_keys = set()
    for md in metadata_files:
        with open(md, 'r', newline='') as mdf:
            rows = csv.DictReader(mdf, delimiter=',')
            for r in rows:
                path, fn = os.path.split(md)
                fid, rn, _ = fn.split('_')
                if r['LibraryLayout'] == 'SINGLE':
                    fqn = os.path.join(path, fid + '.fastq.gz')
                else:
                    if int(rn.strip('r')) == 1:
                        fqn = os.path.join(path, fid + '_1.fastq.gz')
                    elif int(rn.strip('r')) == 2:
                        fqn = os.path.join(path, fid + '_2.fastq.gz')
                    else:
                        raise ValueError('Unexpected read number: {}'.format(rn))
                r['fastq_path'] = fqn
                _, _, _, _, ab_trg, ct, _, sex, brep = r['Subject_ID'].split('_')
                if ct == 'T':
                    r.update(tcell)
                elif ct == 'B':
                    r.update(bcell)
                else:
                    raise ValueError('Unexpected bio sample type: {}'.format(r['Subject_ID']))
                r['Experiment_target'] = ab_trg
                r['Sex'] = sex
                r['Biological_replicate'] = brep
                for k, v in r.items():
                    if not v.strip():
                        r[k] = 'n/a'
                        delete_keys.add(k)
                collection.append(r)
    collection = sorted(collection, key=lambda x: x['Subject_ID'])
    select_keys = set(collection[0].keys()) - delete_keys
    with open(outputfile, 'w', newline='') as dump:
        writer = csv.DictWriter(dump, fieldnames=list(select_keys),
                                restval='n/a', delimiter='\t',
                                extrasaction='ignore')
        writer.writeheader()
        writer.writerows(collection)
    return outputfile


def build_bp_srm_params(mdfile, se_cmd, pe_cmd, outdir, jobcall):
    """
    :param mdfile:
    :param se_cmd:
    :param pe_cmd:
    :param outdir:
    :return:
    """
    args = []
    if not os.path.isfile(mdfile):
        return args
    done = set()
    with open(mdfile, 'r', newline='') as mdf:
        rows = csv.DictReader(mdf, delimiter='\t')
        for r in rows:
            outname = 'bp_mm9_{}_s{}b{}_{}_{}.bam'.format(r['Experiment'],
                                                          r['Sex'][0],
                                                          r['Biological_replicate'],
                                                          r['short_biosample'],
                                                          r['Experiment_target'])
            outpath = os.path.join(outdir, outname)
            metfile = outpath.replace('.bam', '.txt')
            if r['LibraryLayout'] == 'SINGLE':
                inpath = r['fastq_path']
                assert outpath not in done, 'Duplicate creation: {}'.format(outname)
                tmp = se_cmd.format(**{'read1': inpath, 'metricsfile': metfile})
                args.append([inpath, outpath, tmp, jobcall])
                done.add(outpath)
            else:
                if outpath in done:
                    # pair already handled
                    continue
                if '_1.fastq.gz' in r['fastq_path']:
                    r1 = r['fastq_path']
                    r2 = r1.replace('_1.fastq.gz', '_2.fastq.gz')
                else:
                    r2 = r['fastq_path']
                    r1 = r2.replace('_2.fastq.gz', '_1.fastq.gz')
                tmp = pe_cmd.format(**{'read1': r1, 'read2': r2,
                                       'metricsfile': metfile})
                args.append([r1, outpath, tmp, jobcall])
                done.add(outpath)
    return args


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
    bp_metadata = config.get('SpecialFiles', 'bpmetadata')
    bp_sra_acc = config.get('SpecialFiles', 'bp_sra_acc')
    bp_ena_url = config.get('SpecialFiles', 'bp_ena_url')
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
                                                       linkfolder),
                                  name='deepepi_init')

    # ==========================
    # Download Blueprint data
    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        runjob = sci_obj.ruffus_gridjob()
    else:
        runjob = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'bp_fqdl').replace('\n', ' ')
    params_bp_sra_dl = parse_bp_sra_acc(bp_sra_acc, bp_ena_url, dlfolder, cmd, runjob)
    bpdl = pipe.files(sci_obj.get_jobf('in_out'),
                      params_bp_sra_dl,
                      name='bpdl')
    bpdl = bpdl.active_if(os.path.isfile(bp_sra_acc))

    bpann = pipe.merge(task_func=create_blueprint_annotation,
                       name='bpann',
                       input=output_from(bpdl),
                       output=bp_metadata)

    sci_obj.set_config_env(dict(config.items('NodeJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        runjob = sci_obj.ruffus_gridjob()
    else:
        runjob = sci_obj.ruffus_localjob()

    cmd_se = config.get('Pipeline', 'bp_map_se').replace('\n', ' ')
    cmd_pe = config.get('Pipeline', 'bp_map_pe').replace('\n', ' ')
    sr_map_dir = os.path.join(workbase, 'srmap')
    bp_map_params = build_bp_srm_params(bp_metadata, cmd_se, cmd_pe, sr_map_dir, runjob)
    bpmap = pipe.files(sci_obj.get_jobf('in_out'),
                       bp_map_params,
                       name='bpmap')
    bpmap = bpmap.mkdir(sr_map_dir)
    bpmap = bpmap.active_if(os.path.isfile(bp_metadata))
    bpmap = bpmap.follows(bpann)

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

    sigraw_out = os.path.join(workbase, 'conv', 'sigraw')
    cmd = config.get('Pipeline', 'bgtohdfenc').replace('\n', ' ')
    re_filter = '(?P<EID>EE[0-9]+)_(?P<ASSM>(hg19|mm9))_(?P<CELL>\w+)_(?P<MARK>\w+)_[0-9]+\.bg\.gz'
    bgtohdfenc = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                              name='bgtohdfenc',
                              input=output_from(bwtobg),
                              filter=formatter(re_filter),
                              output=os.path.join(sigraw_out, '{EID[0]}_{ASSM[0]}_{CELL[0]}_{MARK[0]}.h5'),
                              extras=[cmd, runjob]).mkdir(sigraw_out)

    sighdf_out = os.path.join(workbase, 'conv', 'hdf')
    cmd = config.get('Pipeline', 'normenc').replace('\n', ' ')
    re_filter = '(?P<EID>EE[0-9]+)_(?P<ASSM>(mm9|hg19))_(?P<CELL>\w+)_(?P<MARK>\w+)\.h5'
    syscalls, exp_output = make_norm_calls(collect_full_paths(sigraw_out, '*EE*.h5'),
                                           sighdf_out, re_filter, cmd, runjob)
    existing_output = os.listdir(sighdf_out)
    normenc = pipe.parallel(sci_obj.get_jobf('raw'),
                            syscalls,
                            name='normenc').active_if(any([fn not in existing_output for fn in exp_output]))

    cmd = config.get('Pipeline', 'bgtohdfdeep').replace('\n', ' ')
    re_filter = '(?P<EID>ED[0-9]+)_(?P<ASSM>hs37d5)_(?P<CELL>\w+)_(?P<MARK>\w+)_[0-9]+\.bg\.gz'
    bgtohdfdeep = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                               name='bgtohdfdeep',
                               input=output_from(bwtobg),
                               filter=formatter(re_filter),
                               output=os.path.join(sighdf_out, '{EID[0]}_hg19_{CELL[0]}_{MARK[0]}.h5'),
                               extras=[cmd, runjob]).mkdir(sighdf_out)

    run_task_convepi = pipe.merge(task_func=touch_checkfile,
                                  name='task_convepi',
                                  input=output_from(bgtohdfenc, bgtohdfdeep,
                                                    normenc),
                                  output=os.path.join(workbase, 'run_task_convepi.chk'))
    #
    # End of major task: convert epigenomes to HDF
    # ================

    return pipe
