# coding=utf-8

import os as os
import csv as csv
import gzip as gz
import io as io
import fnmatch as fnm
import json as js
import re as re

from ruffus import *

from pipelines.auxmods.auxiliary import collect_full_paths, prep_dset_ids, \
    prep_metadata, touch_checkfile


def annotate_encode_files(fpaths, metadata, dsetids):
    """
    :param fpaths:
    :param metadata:
    :param dsetids:
    :return:
    """
    annfiles = []
    groups = dict()
    gcount = 1
    uniq = set()
    for fp in sorted(fpaths):
        facc = os.path.basename(fp).split('.')[0]
        try:
            infos = metadata[facc]
        except KeyError:
            # Non-ENCODE file
            continue
        assert infos['Assembly'] == 'n/a', 'Transcriptome with assembly: {}'.format(infos)
        assert int(infos['Read length']) > 0, 'No read length specified for Fastq {}'.format(facc)
        k = infos['Biosample term name'], infos['Biosample life stage'], infos['Lab'], infos['Experiment accession']
        try:
            tid, assm = dsetids[k]
        except KeyError:
            continue
        infos['Assembly'] = assm
        if infos['Paired with'] != 'n/a':
            partner = infos['Paired with']
            try:
                group = groups[partner]
            except KeyError:
                group = gcount
                groups[facc] = group
                gcount += 1
        else:
            group = 0
        brep = infos['Biological replicate(s)']
        trep = infos['Technical replicate']
        read = infos['Paired end']
        if read == 'n/a':
            read = '0'
        if brep == 'n/a':
            brep = '0'
        if trep == 'n/a':
            trep = '0'
        num_id = 'b{}t{}g{}r{}'.format(brep, trep, group, read)
        if group != 0:
            assert num_id not in uniq, 'Duplicate identifier: {} -- {}'.format(num_id, infos)
        uniq.add(num_id)
        annfiles.append((fp, facc, num_id, tid, infos))
    assert annfiles, 'No ENCODE files annotated'
    return annfiles


def link_encode_transcriptomes(folder, mdfile, idfile, outfolder):
    """
    :param folder:
    :param mdfile:
    :param idfile:
    :param outfolder:
    :return:
    """
    os.makedirs(outfolder, exist_ok=True)
    fqfiles = collect_full_paths(folder, '*.fastq.gz', False)
    metadata = prep_metadata(mdfile, 'File accession')
    tids = prep_dset_ids(idfile, 'ENCODE', 'transcriptome')
    annfiles = annotate_encode_files(fqfiles, metadata, tids)
    linked_files = []
    for afp in annfiles:
        infos = afp[4]
        exptarget = 'mRNA'
        libtype = infos['Run type']
        readlength = infos['Read length']
        new_name = '_'.join([afp[3], infos['Assembly'], infos['Biosample term name'],
                             exptarget, afp[2], libtype + '-' + readlength])
        new_path = os.path.join(outfolder, new_name + '.fastq.gz')
        if not os.path.islink(new_path):
            os.link(afp[0], new_path)
        linked_files.append(new_path)
    assert linked_files, 'No ENCODE files linked to destination: {}'.format(outfolder)
    return linked_files


def annotate_deep_files(fpaths, metadata, dsetids):
    """
    :param fpaths:
    :param metadata:
    :param dsetids:
    :return:
    """
    annfiles = []
    for fp in sorted(fpaths):
        fn = os.path.basename(fp)
        try:
            infos = metadata[fn]
        except KeyError:
            # Non-DEEP file
            continue
        k = infos['biosample'], infos['lifestage'], infos['lab'], 'na'
        try:
            tid, assm = dsetids[k]
        except KeyError:
            continue
        infos['assembly'] = assm
        brep = fn.split('_')[1][3]
        if '_R1_' in fn:
            read = 1
        else:
            read = 2
        num_id = 'b{}t1g1r{}'.format(brep, read)
        annfiles.append((fp, fn, num_id, tid, infos))
    assert annfiles, 'No DEEP files annotated'
    return annfiles


def link_deep_transcriptomes(folder, mdfile, idfile, outfolder):
    """
    :param folder:
    :param mdfile:
    :param idfile:
    :param outfolder:
    :return:
    """
    os.makedirs(outfolder, exist_ok=True)
    fqfiles = collect_full_paths(folder, '*.fastq.gz', False)
    metadata = prep_metadata(mdfile, 'filename')
    tids = prep_dset_ids(idfile, 'DEEP', 'transcriptome')
    annfiles = annotate_deep_files(fqfiles, metadata, tids)
    linked_files = []
    for afp in annfiles:
        infos = afp[4]
        exptarget = 'mRNA'
        libtype = 'pe'
        readlength = '101'
        new_name = '_'.join([afp[3], infos['assembly'], infos['biosample'],
                             exptarget, afp[2], libtype + '-' + readlength])
        new_path = os.path.join(outfolder, new_name + '.fastq.gz')
        if not os.path.islink(new_path):
            os.link(afp[0], new_path)
        linked_files.append(new_path)
    assert linked_files, 'No DEEP files linked to destination: {}'.format(outfolder)
    return linked_files


def annotate_sra_files(fpaths, metadata, dsetids):
    """
    :param fpaths:
    :param metadata:
    :param dsetids:
    :return:
    """
    annfiles = []
    for fp in sorted(fpaths):
        fn, read = os.path.basename(fp).split('.')[0].split('_')
        try:
            infos = metadata[fn]
        except KeyError:
            # non-SRA file
            continue
        if int(infos['use']) < 1:
            continue
        k = infos['cell'], 'un', infos['lab'], infos['experiment']
        try:
            tid, assm = dsetids[k]
        except KeyError:
            continue
        brep = infos['rep']
        num_id = 'b{}t1g0r{}'.format(brep, read)
        annfiles.append((fp, fn, num_id, tid, infos))
    assert annfiles, 'No files annotated for SRA'
    return annfiles


def link_sra_transcriptomes(folder, mdfile, idfile, outfolder):
    """
    :param folder:
    :param mdfile:
    :param idfile:
    :param outfolder:
    :return:
    """
    fqfiles = collect_full_paths(folder, '*.fastq.gz')
    metadata = prep_metadata(mdfile, 'sraid')
    tids = prep_dset_ids(idfile, 'SRA', 'transcriptome')
    annfiles = annotate_sra_files(fqfiles, metadata, tids)
    linked_files = []
    for afp in annfiles:
        infos = afp[4]
        exptarget = 'mRNA'
        libtype = infos['runtype']
        readlength = infos['readlength']
        new_name = '_'.join([afp[3], infos['assembly'], infos['cell'],
                             exptarget, afp[2], libtype + '-' + readlength])
        new_path = os.path.join(outfolder, new_name + '.fastq.gz')
        if not os.path.islink(new_path):
            os.link(afp[0], new_path)
        linked_files.append(new_path)
    assert linked_files, 'No SRA files linked to destination: {}'.format(outfolder)
    return linked_files


def collect_qindices(inpath):
    """
    :param inpath:
    :return:
    """
    # NB: Salmon expects indices to be specified as respective directory
    indexdirs = os.listdir(inpath)
    lookup = dict()
    for idxdir in indexdirs:
        genemodel, km, _ = idxdir.split('.')
        spec, assm, _, version = genemodel.split('_')
        lookup[(assm, km)] = os.path.join(inpath, idxdir), genemodel
    return lookup


def make_qall_calls(fwreads, rvreads, outpath, idxpath, cmd, jobcall):
    """
    :param fwreads:
    :param rvreads:
    :param outpath:
    :param idxpath:
    :param cmd:
    :param jobcall:
    :return:
    """
    infos = re.compile('(?P<SAMPLE>\w+_mRNA)_(?P<NUMID>\w+)_(?P<LIB>[\w\-]+)\.fastq\.gz')
    assert len(fwreads) == len(rvreads), 'Some read files are missing a partner'
    indices = collect_qindices(idxpath)
    current = None
    readlength = 0
    files1, files2 = [], []
    arglist = []
    for read1, read2 in zip(sorted(fwreads), sorted(rvreads)):
        m1, m2 = infos.match(os.path.basename(read1)), infos.match(os.path.basename(read2))
        assert m1.group('SAMPLE') == m2.group('SAMPLE'), 'Sample mismatch: {} vs {}'.format(read1, read2)
        if m1.group('SAMPLE') != current and current is not None:
            assert readlength > 0, 'No read length extracted from filenames'
            assm = current.split('_')[1]
            if readlength < 70:
                idxdir, genemodel = indices[(assm, 'k19')]
            else:
                idxdir, genemodel = indices[(assm, 'k31')]
            genemap = genemodel + '.map.tsv'
            outdir = os.path.join(outpath, current)
            params = {'index': idxdir, 'outpath': outdir, 'genemap': genemap,
                      'reads1': ' '.join(files1), 'reads2': ' '.join(files2)}
            tmp = cmd.format(**params)
            joined_reads = files1 + files2
            arglist.append([joined_reads, os.path.join(outdir, 'quant.genes.sf'), tmp, jobcall])
            files1, files2 = [read1], [read2]
            current = m1.group('SAMPLE')
            readlength = int(m1.group('LIB').split('-')[1])
        elif current is None:
            files1, files2 = [read1], [read2]
            current = m1.group('SAMPLE')
            readlength = int(m1.group('LIB').split('-')[1])
        else:
            files1.append(read1)
            files2.append(read2)
    assert readlength > 0, 'No read length extracted from filenames'
    assm = current.split('_')[1]
    if readlength < 70:
        idxdir, genemodel = indices[(assm, 'k19')]
    else:
        idxdir, genemodel = indices[(assm, 'k31')]
    genemap = genemodel + '.map.tsv'
    outdir = os.path.join(outpath, current)
    params = {'index': idxdir, 'outpath': outdir, 'genemap': genemap,
              'reads1': ' '.join(files1), 'reads2': ' '.join(files2)}
    tmp = cmd.format(**params)
    joined_reads = files1 + files2
    arglist.append([joined_reads, os.path.join(outdir, 'quant.genes.sf'), tmp, jobcall])
    if fwreads:
        assert arglist, 'No argument list created for paired-end quantification'
    return arglist


def convert_salmon_quant(inputfile, outputfile, genemodels):
    """
    :param inputfile:
    :param outputfile:
    :param genemodels:
    :return:
    """
    assert not inputfile.endswith('.gz'), 'Expect gzip input file: {}'.format(inputfile)
    assembly = os.path.basename(outputfile).split('_')[1]
    genemodel = fnm.filter(genemodels, '*_{}_*'.format(assembly))
    assert len(genemodel) == 1, 'Non-unique gene model for assembly {} selected: {}'.format(assembly, genemodel)
    with gz.open(genemodel[0], 'rt', newline='') as gm:
        rows = csv.DictReader(gm, delimiter='\t')
        lut_gene = dict([(r['name'], r) for r in rows])
    outrows = []
    with open(inputfile, 'rt', newline='') as quant:
        rows = csv.DictReader(quant, delimiter='\t')
        for r in rows:
            # by construction, must not throw
            gene = lut_gene[r['Name']]
            gene.update(r)
            outrows.append(gene)
    assert outrows, 'No rows for output created'
    outrows = sorted(outrows, key=lambda d: (d['#chrom'], int(d['start']), int(d['end'])))
    with open(outputfile, 'w') as outf:
        writer = csv.DictWriter(outf, fieldnames=['#chrom', 'start', 'end', 'name', 'TPM', 'strand', 'symbol'],
                                extrasaction='ignore', delimiter='\t')
        writer.writeheader()
        writer.writerows(outrows)
    return outputfile


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
    sra_metadata = config.get('SpecialFiles', 'srametadata')
    dataset_ids = config.get('SpecialFiles', 'datasetids')

    # Init task: link input epigenomes with readable names
    dlfolder = os.path.join(workbase, 'downloads')
    linkfolder = os.path.join(workbase, 'normlinks')
    tempfolder = os.path.join(workbase, 'temp')
    qidxfolder = os.path.join(config.get('Pipeline', 'refdatabase'), 'transcriptome', 'qindex')
    genemodels = os.path.join(config.get('Pipeline', 'refdatabase'), 'genemodel', 'subsets', 'protein_coding')

    enctrans_init = pipe.originate(lambda x: x,
                                   link_encode_transcriptomes(dlfolder, enc_metadata,
                                                              dataset_ids, linkfolder),
                                   name='enctrans_init')

    deeptrans_init = pipe.originate(lambda x: x,
                                    link_deep_transcriptomes(dlfolder, deep_metadata,
                                                             dataset_ids, linkfolder),
                                    name='deeptrans_init')

    sratrans_init = pipe.originate(lambda x: x, collect_full_paths(dlfolder, '*.sra'),
                                   name='sratrans_init')

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'sradump').replace('\n', ' ')
    sradump = pipe.subdivide(task_func=sci_obj.get_jobf('in_pat'),
                             name='sradump',
                             input=output_from(sratrans_init),
                             filter=formatter('(?P<SRARUN>SRR[0-9]+)\.sra'),
                             output=os.path.join(tempfolder, '*{SRARUN[0]}*.fastq.gz'),
                             extras=[tempfolder, '*{SRARUN[0]}*.fastq.gz', cmd, jobcall])
    sradump = sradump.mkdir(tempfolder)

    sratrans_fq = pipe.originate(lambda x: x,
                                 link_sra_transcriptomes(tempfolder, sra_metadata,
                                                         dataset_ids, linkfolder),
                                 name='sratrans_fq')
    sratrans_fq = sratrans_fq.follows(sradump)
    sratrans_fq = sratrans_fq.active_if(len(collect_full_paths(tempfolder, '*.fastq.gz', False)) > 0)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'fastqc_raw')
    reports_raw = os.path.join(workbase, 'reports', 'raw')
    fastqc_raw = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                                name='fastqc_raw',
                                input=output_from(sratrans_fq, enctrans_init, deeptrans_init),
                                filter=formatter('(?P<FILEID>[\w\-]+)\.fastq\.gz'),
                                output=os.path.join(reports_raw, '{FILEID[0]}_fastqc.html'),
                                extras=[cmd, jobcall])
    fastqc_raw = fastqc_raw.mkdir(reports_raw)
    fastqc_raw = fastqc_raw.follows(sratrans_fq)

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    tmpquant = os.path.join(tempfolder, 'quant')

    cmd = config.get('Pipeline', 'qallpe').replace('\n', ' ')
    qallpe = pipe.files(sci_obj.get_jobf('ins_out'),
                        make_qall_calls(collect_full_paths(linkfolder, '*r1_pe*.gz', False),
                                        collect_full_paths(linkfolder, '*r2_pe*.gz', False),
                                        tmpquant, qidxfolder, cmd, jobcall),
                        name='qallpe')
    qallpe = qallpe.mkdir(tmpquant)
    qallpe = qallpe.follows(fastqc_raw)

    sci_obj.set_config_env(dict(config.items('JobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    aggout = os.path.join(workbase, 'conv', 'agg')
    cmd = config.get('Pipeline', 'hdfagg').replace('\n', ' ')
    hdfagg = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                          name='hdfagg',
                          input=output_from(qallpe),
                          filter=formatter('T(S|D|E)[0-9]+_(?P<ASSM>[0-9a-zA-Z]+)_\w+/quant.genes.sf$'),
                          output=os.path.join(aggout, '{ASSM[0]}_agg_exp.genes.h5'),
                          extras=[cmd, jobcall]).mkdir(aggout)

    dir_ortho_pred = os.path.join(config.get('EnvPaths', 'workdir'), 'processing', 'norm', 'task_ortho_pred')
    cmd = config.get('Pipeline', 'odb_pred').replace('\n', ' ')
    odb_pred = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                          name='odb_pred',
                          input=output_from(hdfagg),
                          output=os.path.join(dir_ortho_pred, 'orthopred_odb_6species.h5'),
                          extras=[cmd, jobcall])
    odb_pred = odb_pred.mkdir(dir_ortho_pred)
    odb_orth_file = os.path.join(config.get('Pipeline', 'orthdir'), 'hdf', 'odb9_6species.h5')
    odb_pred = odb_pred.active_if(os.path.isfile(odb_orth_file))

    cmd = config.get('Pipeline', 'hcop_pred').replace('\n', ' ')
    hcop_pred = pipe.merge(task_func=sci_obj.get_jobf('ins_out'),
                           name='hcop_pred',
                           input=output_from(hdfagg),
                           output=os.path.join(dir_ortho_pred, 'orthopred_hcop_6species.h5'),
                           extras=[cmd, jobcall])
    hcop_pred = hcop_pred.mkdir(dir_ortho_pred)
    hcop_orth_file = os.path.join(config.get('Pipeline', 'orthdir'), 'hdf', 'hcop_6species.h5')
    hcop_pred = hcop_pred.active_if(os.path.isfile(hcop_orth_file))
    # use ODB annotation as main
    hcop_pred = hcop_pred.active_if(False)

    bedout = os.path.join(workbase, 'conv', 'bed')
    exp_bed_files = collect_full_paths(bedout, '*.bed.gz', allow_none=True)
    trans_bed_init = pipe.originate(lambda x: x,
                                    exp_bed_files,
                                    name='trans_bed_init')
    trans_bed_init = trans_bed_init.follows(hdfagg)
    trans_bed_init = trans_bed_init.active_if(len(exp_bed_files) > 0)

    hdfout = os.path.join(workbase, 'conv', 'hdf')
    cmd = config.get('Pipeline', 'hdfconv').replace('\n', ' ')
    hdfconv = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                             name='hdfconv',
                             input=output_from(trans_bed_init),
                             filter=suffix('.bed.gz'),
                             output='.h5',
                             output_dir=hdfout,
                             extras=[cmd, jobcall])

    task_preptr = pipe.merge(task_func=touch_checkfile,
                             name='task_preptr',
                             input=output_from(fastqc_raw, qallpe, hdfagg,
                                               trans_bed_init, hdfconv,
                                               odb_pred, hcop_pred),
                             output=os.path.join(workbase, 'run_task_preptr.chk'))

    return pipe
