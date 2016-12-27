# coding=utf-8

import os as os
import csv as csv
import gzip as gz
import io as io
import fnmatch as fnm
import json as js
import sys as sys
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
        if '_Ct1_' in fn:
            brep = 1
        else:
            brep = 2
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
            if readlength > 52:
                idxdir, genemodel = indices[(assm, 'k31')]
            else:
                idxdir, genemodel = indices[(assm, 'k19')]
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
    if fwreads:
        assert arglist, 'No argument list created for paired-end quantification'
    return arglist

# the following functions may be outdated


def salmon_to_bed(inputfile, outputfile, genemodel, datadir, expid):
    """
    :param inputfile:
    :param outputfile:
    :param genemodel:
    :return:
    """
    ra_genemodel = dict([(d['id'], d) for d in js.load(open(genemodel, 'r'))])
    with open(inputfile, 'r', newline='') as quantfile:
        rows = csv.DictReader(quantfile, delimiter='\t')
        for r in rows:
            this_gene = ra_genemodel[r['Name']]
            tpm = float(r['TPM'])
            this_gene['tpm'] = tpm

    fp, fn = os.path.split(genemodel)
    parts = fn.split('_')
    species, assembly, auth, ver = parts[0], parts[1], parts[2], parts[3].split('.')[0]

    assoc_files = os.listdir(datadir)
    assoc_files = fnm.filter(assoc_files, expid + '*' + '.fastq.gz')
    assert assoc_files, 'No files for experiment ID: {}'.format(expid)
    common_name = ''
    for af in assoc_files:
        parts = af.split('_')
        new_name = '_'.join([parts[0], assembly, parts[3], parts[4], parts[5].split('.')[0]])
        if common_name:
            assert new_name == common_name, 'Information mismatch: {} \n {}'.format(common_name, assoc_files)
        common_name = new_name

    outpath = os.path.join(datadir, 'tmp', common_name + '.' + auth + '-' + ver + '.bed')
    genes = sorted(ra_genemodel.values(), key=lambda d: (d['chrom'], d['start'], d['end'], d['id']))
    with open(outpath, 'w') as out:
        writer = csv.DictWriter(out, fieldnames=['chrom', 'start', 'end', 'id', 'tpm', 'strand', 'symbol'],
                                extrasaction='ignore', delimiter='\t')
        writer.writeheader()
        writer.writerows(genes)
    return outputfile


def convert_hcop_table(inputfile, outputfile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    fn = os.path.basename(inputfile)
    from_species = fn.split('_')[0]
    to_species = fn.split('_')[1]
    outbuffer = io.StringIO()
    _ = outbuffer.write('\t'.join([from_species, to_species, 'support', 'assert_ids']) + '\n')
    with gz.open(inputfile, 'rb') as inf:
        # skip header
        _ = inf.readline()
        for line in inf:
            line = line.decode('ascii')
            if not line:
                continue
            _, from_ens, from_assert, _, to_ens, to_assert = line.strip().split()
            if not (from_ens.startswith('ENS') and to_ens.startswith('ENS')):
                continue
            support = str(min(len(from_assert.split(',')), len(to_assert.split(','))))
            _ = outbuffer.write('\t'.join([from_ens, to_ens, support, from_assert + '@' + to_assert]) + '\n')
    with open(outputfile, 'w') as outf:
        _ = outf.write(outbuffer.getvalue())
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
                             extras=[tempfolder, '*{SRARUN[0]}*.fastq.gz', cmd, jobcall]).mkdir(tempfolder)

    sratrans_fq = pipe.originate(lambda x: x,
                                 link_sra_transcriptomes(tempfolder, sra_metadata,
                                                         dataset_ids, linkfolder),
                                 name='sratrans_fq').follows(sradump).active_if(len(collect_full_paths(tempfolder, '*.fastq.gz', False)) > 0)

    sci_obj.set_config_env(dict(config.items('MemJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'fastqc')
    fastqc = pipe.transform(task_func=sci_obj.get_jobf('in_out'),
                            name='fastqc',
                            input=output_from(sratrans_fq, enctrans_init, deeptrans_init),
                            filter=formatter('(?P<FILEID>[\w\-]+)\.fastq\.gz'),
                            output=os.path.join(workbase, 'qcrep', '{FILEID[0]}_fastqc.zip'),
                            extras=[cmd, jobcall])

    sci_obj.set_config_env(dict(config.items('ParallelJobConfig')), dict(config.items('CondaPPLCS')))
    if args.gridmode:
        jobcall = sci_obj.ruffus_gridjob()
    else:
        jobcall = sci_obj.ruffus_localjob()

    cmd = config.get('Pipeline', 'qmmuk13').replace('\n', ' ')
    tmpquant = os.path.join(tempfolder, 'quant')
    qmmuk13 = pipe.collate(task_func=sci_obj.get_jobf('ins_out'),
                           name='qmmuk13',
                           input=output_from(enctrans_init),
                           filter=formatter('(?P<SAMPLE>\w+_mRNA)\w+(?P<LIB>se\-[0-9]+)\.fastq\.gz'),
                           output=os.path.join(tmpquant, '{SAMPLE[0]}', 'quant.genes.sf'),
                           extras=[cmd, jobcall]).mkdir(tmpquant)

    cmd = config.get('Pipeline', 'qallpe').replace('\n', ' ')
    qallpe = pipe.files(sci_obj.get_jobf('ins_out'),
                        make_qall_calls(collect_full_paths(linkfolder, '*r1_pe*.gz', False),
                                        collect_full_paths(linkfolder, '*r2_pe*.gz', False),
                                        tmpquant, qidxfolder, cmd, jobcall))

    task_preptr = pipe.merge(task_func=touch_checkfile,
                             name='task_preptr',
                             input=output_from(fastqc, qmmuk13, qallpe),
                             output=os.path.join(workbase, 'run_task_preptr.chk'))

    # orth_files = os.listdir(orthdir)
    # orth_files = [os.path.join(orthdir, f) for f in orth_files]
    #
    # initorth = pipe.originate(lambda x: x, orth_files, name='initorth')
    #
    # convorth = pipe.transform(task_func=convert_hcop_table,
    #                           name='convorth',
    #                           input=output_from(initorth),
    #                           filter=suffix('_six_column.txt.gz'),
    #                           output='_orthologs.tsv',
    #                           output_dir=orthdir).jobs_limit(2)

    return pipe
