
include: 'aux_utilities.smk'

localrules: compute_insert_size_estimates

SINGLE_EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.single.ok')
PAIRED_EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.paired1.ok')

EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{runid}.{layout,(single|paired)}{mate,[12]{0,1}}.ok')
HISTONE_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark,H3K[a-z0-9]+}_{sample,ET\-[A-Z0-9]+}_{runid}.{layout}.ok')

rule preprocess_epigenomes_master:
    input:
        expand('input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.read-length',
                zip,
                layout=EPIGENOME_WILDCARDS.layout,
                species=EPIGENOME_WILDCARDS.species,
                tissue=EPIGENOME_WILDCARDS.tissue,
                mark=EPIGENOME_WILDCARDS.mark,
                sample=EPIGENOME_WILDCARDS.sample),

        expand('input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.bg',
                zip,
                species=EPIGENOME_WILDCARDS.species,
                tissue=EPIGENOME_WILDCARDS.tissue,
                mark=EPIGENOME_WILDCARDS.mark,
                sample=EPIGENOME_WILDCARDS.sample),

        expand('input/bed/epigenome/{species}_{tissue}_{mark}_{sample}.peaks.bed',
                zip,
                species=HISTONE_WILDCARDS.species,
                tissue=HISTONE_WILDCARDS.tissue,
                mark=HISTONE_WILDCARDS.mark,
                sample=HISTONE_WILDCARDS.sample),

        expand('input/bedgraph/epigenome/peak_signal/{species}_{tissue}_{mark}_{sample}.cov.bg',
                zip,
                species=HISTONE_WILDCARDS.species,
                tissue=HISTONE_WILDCARDS.tissue,
                mark=HISTONE_WILDCARDS.mark,
                sample=HISTONE_WILDCARDS.sample)


def collect_single_end_epigenomes(wildcards):

    def filter_fun(filename):
        keep = wildcards.species in filename
        keep &= wildcards.tissue in filename
        keep &= wildcards.mark in filename
        keep &= wildcards.sample in filename
        keep &= filename.endswith('.single.ok')
        return keep

    epigenome_path = 'input/fastq/epigenome'
    process_fastq = []
    for root, dirs, files in os.walk(epigenome_path):
        if not files:
            continue
        fastq_files = list(filter(filter_fun, files))
        fastq_files = [os.path.join(root, f.replace('.ok', '.fastq.gz')) for f in fastq_files]
        process_fastq.extend(fastq_files)
    if not process_fastq:
        raise ValueError('No fastq/OK files for wildcards: {}'.format(wildcards))
    return sorted(process_fastq)


rule align_single_end_epigenome:
    input:
        fastq = collect_single_end_epigenomes,
        index = expand('references/indices/bowtie/{{species}}/{{species}}.{idxnum}.bt2',
                        idxnum=[1, 2, 3, 4]),
        autosome_regions = 'references/chromosomes/autosomes/{species}.auto.bed',
    output:
        bam = 'input/bam/epigenome/temp/single/{species}_{tissue}_{mark}_{sample}.filt.bam',
        metrics = 'input/bam/epigenome/metrics/{species}_{tissue}_{mark}_{sample}.bowtie2.metrics'
    log:
        bowtie = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.bowtie.log',
        samtools = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.samtools.log'
    benchmark: 'run/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.aln.rsrc'
    threads: 16
    params:
        index_loc = lambda wildcards: 'references/indices/bowtie/{}/{}'.format(wildcards.species, wildcards.species),
        min_mapq = config['min_mapq']
    run:
        if isinstance(input.fastq, list):
            input_files = ' -U ' + ' -U '.join(input.fastq)
        else:
            input_files = ' -U ' + ' -U '.join(list(input.fastq))

        exec = 'bowtie2 -q --very-sensitive-local'
        exec += ' --met-file {output.metrics}'
        exec += ' --threads {threads} --qc-filter'
        exec += ' -x {params.index_loc}'
        exec += input_files
        exec += ' -S /dev/stdout'
        exec += ' 2> {log.bowtie}'
        exec += ' | samtools view'
        exec += ' -b -q {params.min_mapq}'  # output BAM
        exec += ' -F 3844'  # exclude unmapped, secondary, failed QC, duplicate
        exec += ' -L {input.autosome_regions}'
        exec += ' -o {output.bam} -@ 2'
        exec += ' /dev/stdin'
        exec += ' &> {log.samtools}'
        shell(exec)


def collect_paired_end_epigenomes_mate1(wildcards):

    def filter_fun(filename):
        keep = wildcards.species in filename
        keep &= wildcards.tissue in filename
        keep &= wildcards.mark in filename
        keep &= wildcards.sample in filename
        keep &= filename.endswith('.paired1.ok')
        return keep

    epigenome_path = 'input/fastq/epigenome'
    process_fastq_pair1 = []
    for root, dirs, files in os.walk(epigenome_path):
        if not files:
            continue
        fastq_files = list(filter(filter_fun, files))
        for f in fastq_files:
            file_path = os.path.join(root, f.replace('.ok', '.fastq.gz'))
            if 'paired1' in file_path:
                process_fastq_pair1.append(file_path)
            else:
                raise ValueError('Unexpected mate number not detected: {}'.format(file_path))
    if not process_fastq_pair1:
        raise ValueError('No fastq/OK files for wildcards: {}'.format(wildcards))
    return sorted(process_fastq_pair1)


def collect_paired_end_epigenomes_mate2(wildcards):

    def filter_fun(filename):
        keep = wildcards.species in filename
        keep &= wildcards.tissue in filename
        keep &= wildcards.mark in filename
        keep &= wildcards.sample in filename
        keep &= filename.endswith('.paired2.ok')
        return keep

    epigenome_path = 'input/fastq/epigenome'
    process_fastq_pair2 = []
    for root, dirs, files in os.walk(epigenome_path):
        if not files:
            continue
        fastq_files = list(filter(filter_fun, files))
        for f in fastq_files:
            file_path = os.path.join(root, f.replace('.ok', '.fastq.gz'))
            if 'paired2' in file_path:
                process_fastq_pair2.append(file_path)
            else:
                raise ValueError('Unexpected mate number not detected: {}'.format(file_path))
    if not process_fastq_pair2:
        raise ValueError('No fastq/OK files for wildcards: {}'.format(wildcards))
    return sorted(process_fastq_pair2)


rule align_paired_end_epigenome:
    input:
        fastq_mate1 = collect_paired_end_epigenomes_mate1,
        fastq_mate2 = collect_paired_end_epigenomes_mate2,
        index = expand('references/indices/bowtie/{{species}}/{{species}}.{idxnum}.bt2',
                        idxnum=[1, 2, 3, 4]),
        autosome_regions = 'references/chromosomes/autosomes/{species}.auto.bed',
    output:
        bam = 'input/bam/epigenome/temp/paired/{species}_{tissue}_{mark}_{sample}.filt.bam',
        metrics = 'input/bam/epigenome/metrics/{species}_{tissue}_{mark}_{sample}.bowtie2.metrics'
    log:
        bowtie = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.bowtie.log',
        samtools = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.samtools.log'
    benchmark: 'run/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.aln.rsrc'
    threads: 16
    params:
        index_loc = lambda wildcards: 'references/indices/bowtie/{}/{}'.format(wildcards.species, wildcards.species),
        min_mapq = config['min_mapq']
    run:
        if isinstance(input.fastq_mate1, list):
            input_mate1 = ' -1 ' + ' -1 '.join(input.fastq_mate1)
            input_mate2 = ' -2 ' + ' -2 '.join(input.fastq_mate2)
        else:
            input_mate1 = ' -1 ' + ' -1 '.join(list(input.fastq_mate1))
            input_mate2 = ' -2 ' + ' -2 '.join(list(input.fastq_mate2))

        assert len(input_mate1) == len(input_mate2), 'Missing mate file:\n{}\n\n{}\n'.format(input_mate1, input_mate2)

        exec = 'bowtie2 -q --very-sensitive-local'
        exec += ' --met-file {output.metrics}'
        exec += ' --threads {threads} --qc-filter'
        exec += ' -x {params.index_loc}'
        exec += input_mate1
        exec += input_mate2
        exec += ' -S /dev/stdout'
        exec += ' 2> {log.bowtie}'
        exec += ' | samtools view'
        exec += ' -b -q {params.min_mapq}'  # output BAM
        exec += ' -f 2'  # include mapped in proper pair
        exec += ' -F 3844'  # exclude unmapped, secondary, failed QC, duplicate
        exec += ' -L {input.autosome_regions}'
        exec += ' -o {output.bam} -@ 2'
        exec += ' /dev/stdin'
        exec += ' &> {log.samtools}'
        shell(exec)


rule compute_bam_read_stats:
    input:
        'input/bam/epigenome/temp/{layout}/{filename}.filt.bam'
    output:
        'input/bam/epigenome/temp/{layout}/{filename}.read-stats'
    shell:
        'samtools stats {input} > {output}'


rule compute_insert_size_estimates:
    input:
        expand('input/bam/epigenome/temp/paired/{filename}.read-stats',
                filename=glob_wildcards('input/bam/epigenome/temp/paired/{filename}.read-stats').filename)
    output:
        'annotation/data_statistics/paired-end_insert_sizes.tsv'
    log:
        'log/annotation/data_statistics/paired-end_insert_sizes.log'
    run:
        import numpy as np

        collector = {
            'H3K4me3': {
                'weights': [],
                'sizes': []
                },
            'H3K27ac': {
                'weights': [],
                'sizes': []
                },
            'H3K36me3': {
                'weights': [],
                'sizes': []
                },
            'Input': {
                'weights': [],
                'sizes': []
            }
        }

        with open(log[0], 'w') as logfile:
            _ = logfile.write('Log started')
            for input_file in input:
                _ = logfile.write('Processing: {}\n'.format(input_file))
                mark = input_file.split('_')[2]
                _ = logfile.write('Mark identified as: {}\n'.format(mark))
                assert mark in collector, 'Unexpected data type: {}'.format(mark)
                with open(input_file, 'r') as stats:
                    for line in stats:
                        if line.startswith('IS'):
                            _ = logfile.write('Processing line...' + '\n')
                            columns = line.strip().split('\t')
                            _ = logfile.write('Columns extracted: {}'.format('_'.join(columns)) + '\n')
                            size = int(columns[1])
                            weight = int(columns[2])
                            _ = logfile.write('Storing values')
                            collector[mark]['sizes'].append(size)
                            collector[mark]['weights'].append(weight)

        with open(output[0], 'w') as table:
            for mark, stats in collector.items():
                sizes = np.array(stats['sizes'], dtype=np.int16)
                weights = np.array(stats['weights'], dtype=np.int32)
                size_estimate = np.average(sizes, weights=weights)
                size_estimate = np.floor(size_estimate)
                _ = table.write('{}\t{}\n'.format(mark, int(size_estimate)))
    # end of rule


rule extract_read_length:
    input:
        'input/bam/epigenome/temp/{layout}/{filename}.read-stats'
    output:
        'input/bam/epigenome/temp/{layout}/{filename}.read-length'
    shell:
        'egrep "^RL" {input} |  cut -f 2 > {output}'
    # end of rule


rule extract_old_bam_header:
    input:
        'input/bam/epigenome/temp/{layout}/{filename}.filt.bam'
    output:
        'input/bam/epigenome/temp/{layout}/{filename}.filt.old_header.sam'
    shell:
        'samtools view -H {input} > {output}'


rule adjust_old_header:
    input:
        header = 'input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.filt.old_header.sam',
        chrom = 'references/chromosomes/autosomes/{species}.auto.sizes'
    output:
        'input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.filt.new_header.sam'
    run:
        import io as io

        autosomes = set()
        with open(input.chrom, 'r') as chromosomes:
            [autosomes.add(l.split()[0]) for l in chromosomes.readlines()]

        new_header = io.StringIO()
        with open(input.header, 'r') as header:
            for line in header:
                if line.startswith('@'):
                    if line.startswith('@SQ'):
                        _, chrom, _ = line.strip().split()
                        chrom = chrom.replace('SN:', '')
                        if chrom in autosomes:
                            new_header.write(line)
                    else:
                        new_header.write(line)
                else:
                    new_header.write(line)

        with open(output[0], 'w') as dump:
            _ = dump.write(new_header.getvalue())
    # end of rule


rule reheader_bam_file:
    input:
        bam = 'input/bam/epigenome/temp/{layout}/{filename}.filt.bam',
        header = 'input/bam/epigenome/temp/{layout}/{filename}.filt.new_header.sam'
    output:
        'input/bam/epigenome/temp/{layout}/{filename}.filt.rhd.bam'
    shell:
        'samtools reheader {input.header} {input.bam} > {output}'


rule collect_mapped_read_statistics:
    input:
        expand('input/bam/epigenome/temp/{layout}/{filename}.read-stats',
                zip,
                layout=glob_wildcards('input/bam/epigenome/temp/{layout}/{filename}.read-stats').layout,
                filename=glob_wildcards('input/bam/epigenome/temp/{layout}/{filename}.read-stats').filename)
    output:
        'annotation/data_statistics/high-qual_mapped_reads.tsv'
    run:
        out_rows = []
        for statfile in input:
            species, tissue, mark, sample = os.path.basename(statfile).split('.')[0].split('_')
            layout = os.path.split(os.path.dirname(statfile))[1]
            with open(statfile, 'r') as stats:
                for line in stats:
                    if line.startswith('SN') and 'reads mapped' in line:
                        num_reads = int(line.split()[3])
                        assert num_reads > 0, 'Unexpected number of reads: {} - {}'.format(line.strip(), statfile)
                        out_rows.append((species, tissue, mark, layout, str(num_reads), sample))
                        break
        out_rows = sorted(out_rows, key=lambda x: (x[2], x[3], x[0]))
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(['\t'.join(r) for r in out_rows]) + '\n')
    # end of rule


rule subsample_bam_files:
    input:
        read_counts = 'annotation/data_statistics/high-qual_mapped_reads.tsv',
        bam = 'input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.filt.rhd.psort.bam',
        bam_idx = 'input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.filt.rhd.psort.bam.csi',
    output:
        'input/bam/epigenome/sampled/{layout}/{species}_{tissue}_{mark}_{sample}.sub.bam'
    log:
        'log/input/bam/epigenome/sampled/{layout}/{species}_{tissue}_{mark}_{sample}.sub.log'
    benchmark:
        'run/input/bam/epigenome/sampled/{layout}/{species}_{tissue}_{mark}_{sample}.sub.rsrc'
    params:
        limit_narrow = config['min_reads_narrow'],
        limit_broad = config['min_reads_broad']
    threads: 2
    run:
        import csv
        header = ['species', 'tissue', 'mark', 'layout', 'num_reads', 'sample']
        num_reads = 0
        with open(input.read_counts, 'r', newline='') as table:
            reader = csv.DictReader(table, delimiter='\t', fieldnames=header)
            for row in reader:
                if row['sample'] == wildcards.sample and row['mark'] == wildcards.mark:
                    num_reads = int(row['num_reads'])
                    break
        assert num_reads > 0, 'No read count for sample: {}'.format(wildcards)

        if wildcards.mark in ['H3K36me3', 'Input']:
            min_limit = params.limit_broad
        else:
            min_limit = params.limit_narrow

        if num_reads <= min_limit:
            subsample = ' '
        else:
            fraction = round(min_limit / num_reads, 3)
            assert 0 < fraction < 1, 'This is...unexpected: {} - {}'.format(fraction, wildcards)
            subsample = ' -s {} '.format(fraction)
        exec = 'samtools view -b '
        exec += subsample
        exec += ' -o {output} --output-fmt BAM '
        exec += ' --threads {threads}'
        exec += ' {input.bam} &>> {log}'
        with open(log[0], 'w') as logfile:
            _ = logfile.write(exec + '\n\n')
        shell(exec)


def find_signal_coverage_bam(wildcards):

    bam_root = 'input/bam/epigenome/temp'

    file_id = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, wildcards.sample])

    layout = None
    for root, dirs, files in os.walk('input/fastq/epigenome'):
        sample_fastq = list(filter(lambda x: file_id in x, files))
        for s in sample_fastq:
            if s.endswith('single.ok'):
                layout = 'single'
                break
            elif s.endswith('paired1.ok'):
                layout = 'paired'
                break
            else:
                continue
        if layout is not None:
            break

    if layout is None:
        raise ValueError('No lib layout for sample: {}'.format(wildcards))

    input_bam = os.path.join(bam_root, layout, file_id + '.filt.rhd.psort.bam')
    return input_bam


def find_signal_coverage_bam_index(wildcards):

    bam_root = 'input/bam/epigenome/temp'

    file_id = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, wildcards.sample])

    layout = None
    for root, dirs, files in os.walk('input/fastq/epigenome'):
        sample_fastq = list(filter(lambda x: file_id in x, files))
        for s in sample_fastq:
            if s.endswith('single.ok'):
                layout = 'single'
                break
            elif s.endswith('paired1.ok'):
                layout = 'paired'
                break
            else:
                continue
        if layout is not None:
            break

    if layout is None:
        raise ValueError('No lib layout for sample: {}'.format(wildcards))

    input_bai = os.path.join(bam_root, layout, file_id + '.filt.rhd.psort.bam.csi')
    return input_bai


def find_signal_coverage_gsize(wildcards):

    gsize = 'references/assemblies/autosomes/{species}.k{kmer}.stats'

    file_id = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, wildcards.sample])

    single_path = 'input/bam/epigenome/temp/single'
    paired_path = 'input/bam/epigenome/temp/paired'

    kmer = None
    for path in [single_path, paired_path]:

        if not os.path.isdir(path):
            raise ValueError('Checkpoint not triggered: {}'.format(wildcards))

        rlen_file = list(filter(lambda x: file_id in x and x.endswith('read-length'), os.listdir(path)))
        for rf in rlen_file:
            rlen_path = os.path.join(path, rf)
            with open(rlen_path, 'r') as dump:
                kmer = dump.readline().strip()
            break
        if kmer is not None:
            break

    if kmer is None:
        raise ValueError('No read length file for wildcards: {}'.format(wildcards))

    gsize = gsize.format(**{'species': wildcards.species,
                            'kmer': kmer})
    return gsize


def find_signal_coverage_kmer_length(wildcards):

    file_id = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, wildcards.sample])

    single_path = 'input/bam/epigenome/temp/single'
    paired_path = 'input/bam/epigenome/temp/paired'

    rlen_path = None
    for path in [single_path, paired_path]:

        if not os.path.isdir(path):
            raise ValueError('Checkpoint not triggered: {}'.format(wildcards))

        rlen_file = list(filter(lambda x: file_id in x and x.endswith('read-length'), os.listdir(path)))
        for r in rlen_file:
            rlen_path = os.path.join(path, r)
            break
        if rlen_path is not None:
            break
    return rlen_path


rule compute_signal_coverage:
    input:
        bam = find_signal_coverage_bam,
        bai = find_signal_coverage_bam_index,
        genome_size = find_signal_coverage_gsize,
        read_length = find_signal_coverage_kmer_length,
        insert_size_estimates = 'annotation/data_statistics/paired-end_insert_sizes.tsv',
    output:
        'input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.bg'
    threads: 6
    log:
        'log/input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.log'
    benchmark:
        'run/input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.rsrc'
    wildcard_constraints:
        species = '(' + '|'.join(config['species']) + ')'
    run:
        is_est = None
        with open(log[0], 'w') as logfile:
            _ = logfile.write('Processing file: {}'.format(input.bam) + '\n')
            effective_genome_size = None
            with open(input.genome_size, 'r') as counts:
                for line in counts:
                    if line.startswith('number of unique k-mers'):
                        effective_genome_size = line.split('\t')[1].strip()
                        try:
                            assert int(effective_genome_size), 'Non-integer genome size: {}'.format(input.genome_size)
                        except (AssertionError, ValueError):
                            _ = logfile.write('Invalid eff. genome size: {}'.format(line.strip()) + '\n')
                            raise
                        else:
                            break
            assert effective_genome_size is not None, 'No genome size for file: {}'.format(input.genome_size)
            _ = logfile.write('Effective genome size set to: {}'.format(effective_genome_size) + '\n')

            if 'epigenome/temp/single' in input.bam:
                mark = os.path.basename(input.bam).split('_')[2]
                _ = logfile.write('Single-end epigenome, mark is: {}'.format(mark) + '\n')
                with open(input.insert_size_estimates, 'r') as estimates:
                    for line in estimates:
                        if line.startswith(mark):
                            _, is_est = line.strip().split()
                            break
                if is_est is None:
                    _ = logfile.write('No insert size estimate extracted' + '\n')
                    raise ValueError
                else:
                    _ = logfile.write('Insert size estimate set to: {}'.format(is_est) + '\n')

            exec = 'bamCoverage --bam {input.bam}'
            exec += ' --outFileName {output}'  # for whatever reason, cannot write to /dev/stdout
            exec += ' --outFileFormat bedgraph'
            exec += ' --binSize 25 --numberOfProcessors {threads}'
            exec += ' --effectiveGenomeSize ' + effective_genome_size
            exec += ' --normalizeUsing RPGC' # according to docs, should be 1x normalization
            if is_est is not None:
                exec += ' --extendReads ' + is_est
            # Note to self:
            # putting all non-primary chromosomes into the ignore set
            # fails for certain assemblies such as macaque because
            # of the character limit for single command line calls
            #
            #exec += ' --ignoreForNormalization ' + ignore_set

            exec += ' --verbose &>> {log}'
            _ = logfile.write('\n' + exec + '\n\n')

        shell(exec)


def load_effective_genome_size(wildcards, input):
    """
    """
    effective_genome_size = None
    with open(input.genome_size, 'r') as counts:
        for line in counts:
            if line.startswith('number of unique k-mers'):
                effective_genome_size = line.split('\t')[1].strip()
                try:
                    assert int(effective_genome_size), 'Non-integer genome size: {}'.format(input.genome_size)
                except (AssertionError, ValueError):
                    _ = logfile.write('Invalid eff. genome size: {}'.format(line.strip()) + '\n')
                    raise
                else:
                    break
    assert effective_genome_size is not None, 'No genome size for file: {}'.format(input.genome_size)
    return effective_genome_size


def load_extension_size(wildcards, input):
    """
    """
    insert_size = None
    with open(input.insert_size_estimates, 'r') as table:
        for line in table:
            if line.startswith(wildcards.mark):
                _, insert_size = line.strip().split()
                break
    assert insert_size is not None, 'Could not load insert size estimate: {}'.format(wildcards)
    return insert_size


def select_paired_input_control(wildcards):

    import glob

    file_pattern = '_'.join([wildcards.species, wildcards.tissue, 'Input', wildcards.sample]) + '.read-stats'

    paired_path = 'input/bam/epigenome/temp/paired'
    search_pattern = os.path.join(paired_path, file_pattern)

    matched_input = glob.glob(search_pattern, recursive=False)
    if len(matched_input) == 0:
        file_pattern = '_'.join([wildcards.species, wildcards.tissue, 'Input', 'ET-*']) + '.read-stats'
        search_pattern = os.path.join(paired_path, file_pattern)
        unmatched_input = glob.glob(search_pattern, recursive=False)
        if len(unmatched_input) == 0:
            raise ValueError('No input found: {}'.format(wildcards))
        elif len(unmatched_input) == 1:
            selected_input = unmatched_input[0]
        else:
            file_pattern = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, '*']) + '.read-stats'
            search_pattern = os.path.join(paired_path, file_pattern)
            all_mark_files = glob.glob(search_pattern, recursive=False)
            if not len(all_mark_files) == len(unmatched_input):
                raise ValueError('PE: Unequal set size marks and inputs: {} v {}'.format(all_mark_files, unmatched_input))
            for m, i in zip(sorted(all_mark_files), sorted(unmatched_input)):
                if wildcards.sample in m:
                    selected_input = i
                    break

    else:
        selected_input = matched_input[0]

    selected_input = selected_input.replace('/temp/', '/sampled/').replace('read-stats', 'sub.nsort.bam')

    return selected_input


rule call_paired_end_peaks_macs:
    input:
        signal = 'input/bam/epigenome/sampled/paired/{species}_{tissue}_{mark}_{sample}.sub.nsort.bam',
        chromosomes = 'references/chromosomes/autosomes/{species}.auto.sizes',
        background = select_paired_input_control,
        genome_size = find_signal_coverage_gsize,
        read_length = 'input/bam/epigenome/temp/paired/{species}_{tissue}_{mark}_{sample}.read-length',
        insert_size_estimates = 'annotation/data_statistics/paired-end_insert_sizes.tsv',
    output:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.pe_peaks.xls'
    log:
        'log/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.pe_peaks.log'
    benchmark:
        'run/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.pe_peaks.rsrc'
    wildcard_constraints:
        mark = '(H3K4me3|H3K36me3|H3K27ac)',
        sample = 'ET\-[0-9A-Z]+'
    params:
        effgensize = load_effective_genome_size,
        extsize = load_extension_size,
        is_broad = lambda wildcards: ' --broad --broad-cutoff 0.1 ' if wildcards.mark == 'H3K36me3' else '',
        out_dir = lambda wildcards, output: os.path.dirname(output[0]),
        out_name = lambda wildcards, output: os.path.basename(output[0]).rsplit('_', 1)[0]
    conda: '../environment/conda/conda_cse_py2.txt'
    shell:
        'macs2 callpeak --treatment {input.signal} --control {input.background}' \
           ' --name {params.out_name} --outdir {params.out_dir}' \
           ' --format BAMPE --gsize {params.effgensize}' \
           ' --qvalue 0.05 --keep-dup 1 --nomodel' \
           ' --extsize {params.extsize} {params.is_broad}' \
           ' --scale-to small &> {log}'


rule call_paired_end_peaks_zerone:
    input:
        signal = 'input/bam/epigenome/sampled/paired/{species}_{tissue}_{mark}_{sample}.sub.nsort.bam',
        chromosomes = 'references/chromosomes/autosomes/{species}.auto.sizes',
        background = select_paired_input_control
    output:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.pe_peaks.bed'
    log:
        'log/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.pe_peaks.log'
    benchmark:
        'run/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.pe_peaks.rsrc'
    wildcard_constraints:
        mark = '(H3K4me3|H3K36me3|H3K27ac)',
        sample = 'ET\-[0-9A-Z]+'
    shell:
        'zerone --chip {input.signal} --mock {input.background} --confidence 0.95 1> {output} 2> {log}'


def select_single_input_control(wildcards):

    import glob

    file_pattern = '_'.join([wildcards.species, wildcards.tissue, 'Input', wildcards.sample]) + '.read-stats'

    paired_path = 'input/bam/epigenome/temp/single'
    search_pattern = os.path.join(paired_path, file_pattern)

    matched_input = glob.glob(search_pattern, recursive=False)
    if len(matched_input) == 0:
        file_pattern = '_'.join([wildcards.species, wildcards.tissue, 'Input', 'ET-*']) + '.read-stats'
        search_pattern = os.path.join(paired_path, file_pattern)
        unmatched_input = glob.glob(search_pattern, recursive=False)
        if len(unmatched_input) == 0:
            raise ValueError('No input found: {}'.format(wildcards))
        elif len(unmatched_input) == 1:
            selected_input = unmatched_input[0]
        else:
            file_pattern = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, '*']) + '.read-stats'
            search_pattern = os.path.join(paired_path, file_pattern)
            all_mark_files = glob.glob(search_pattern, recursive=False)
            if not len(all_mark_files) == len(unmatched_input):
                raise ValueError('SE: Unequal set size marks and inputs: {} v {}'.format(all_mark_files, unmatched_input))
            for m, i in zip(sorted(all_mark_files), sorted(unmatched_input)):
                if wildcards.sample in m:
                    selected_input = i
                    break

    else:
        selected_input = matched_input[0]

    selected_input = selected_input.replace('/temp/', '/sampled/').replace('read-stats', 'sub.psort.bam')

    return selected_input


rule call_single_end_peaks_macs:
    input:
        signal = 'input/bam/epigenome/sampled/single/{species}_{tissue}_{mark}_{sample}.sub.psort.bam',
        chromosomes = 'references/chromosomes/autosomes/{species}.auto.sizes',
        background = select_single_input_control,
        genome_size = find_signal_coverage_gsize,
        read_length = 'input/bam/epigenome/temp/single/{species}_{tissue}_{mark}_{sample}.read-length',
        insert_size_estimates = 'annotation/data_statistics/paired-end_insert_sizes.tsv',
    output:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.se_peaks.xls'
    log:
        'log/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.se_peaks.log'
    benchmark:
        'run/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.se_peaks.rsrc'
    wildcard_constraints:
        mark = '(H3K4me3|H3K36me3|H3K27ac)',
        sample = 'ET\-[0-9A-Z]+'
    params:
        effgensize = load_effective_genome_size,
        extsize = load_extension_size,
        is_broad = lambda wildcards: ' --broad --broad-cutoff 0.1 ' if wildcards.mark == 'H3K36me3' else '',
        out_dir = lambda wildcards, output: os.path.dirname(output[0]),
        out_name = lambda wildcards, output: os.path.basename(output[0]).rsplit('_', 1)[0]
    conda: '../environment/conda/conda_cse_py2.txt'
    shell:
        'macs2 callpeak --treatment {input.signal} --control {input.background}' \
           ' --name {params.out_name} --outdir {params.out_dir}' \
           ' --format BAM --gsize {params.effgensize}' \
           ' --qvalue 0.05 --keep-dup 1 --nomodel' \
           ' --extsize {params.extsize} {params.is_broad}' \
           ' --scale-to small &> {log}'


rule call_single_end_peaks_zerone:
    input:
        signal = 'input/bam/epigenome/sampled/single/{species}_{tissue}_{mark}_{sample}.sub.psort.bam',
        background = select_single_input_control,
    output:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.se_peaks.bed'
    log:
        'log/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.se_peaks.log'
    benchmark:
        'run/input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.se_peaks.rsrc'
    wildcard_constraints:
        mark = '(H3K4me3|H3K36me3|H3K27ac)',
        sample = 'ET\-[0-9A-Z]+'
    shell:
        'zerone --chip {input.signal} --mock {input.background} --confidence 0.95 1> {output} 2> {log}'


rule copy_peak_file:
    input:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.{layout_id}_peaks.xls'
    output:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.{layout_id}_peaks.bed',
    run:
        alt_input = input[0]
        if wildcards.mark == 'H3K36me3':
            alt_input = alt_input.replace('.xls', '.broadPeak')
        else:
            alt_input = alt_input.replace('.xls', '.narrowPeak')

        shell('cp {} {{output}}'.format(alt_input))


rule intersect_peak_calls:
    input:
        macs = 'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/macs2.{layout_id}_peaks.bed',
        zerone = 'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/zerone.{layout_id}_peaks.bed'
    output:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/shared.{layout_id}_peaks.tsv'
    shell:
        "bedtools intersect -wo -a {input.macs} -b {input.zerone} > {output}"


rule concat_shared_peaks:
    input:
        'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/shared.{layout_id}_peaks.tsv'
    output:
        concat = temp('input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/concat.{layout_id}_peaks.bed'),
        sort = 'input/bed/epigenome/temp/{species}_{tissue}_{mark}_{sample}/sorted.{layout_id}_peaks.bed',
    params:
        cut_cols = lambda wildcards: '11,12,13' if wildcards.mark in ['H3K4me3', 'H3K27ac'] else '10,11,12'
    shell:
        'cut -f 1,2,3 {input} > {output.concat} && cut -f {params.cut_cols} {input} >> {output.concat} '\
            ' && cat {output.concat} | sort -V -k1 -k2n,3n > {output.sort}'


def locate_peak_file(wildcards):
    """
    """
    import glob
    single_path = 'input/bam/epigenome/temp/single'
    paired_path = 'input/bam/epigenome/temp/paired'

    file_pattern = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, wildcards.sample])

    sample_peak_path = os.path.join('input/bed/epigenome/temp', file_pattern)

    lookup_pattern = file_pattern + '.*'
    folder_pattern = os.path.join(single_path, lookup_pattern)

    matched_bam = glob.glob(folder_pattern, recursive=False)
    if len(matched_bam) > 0:
        return os.path.join(sample_peak_path, 'sorted.se_peaks.bed')

    folder_pattern = os.path.join(paired_path, lookup_pattern)
    matched_bam = glob.glob(folder_pattern, recursive=False)

    if len(matched_bam) > 0:
        return os.path.join(sample_peak_path, 'sorted.pe_peaks.bed')
    raise ValueError('No peak file for wildcards {}'.format(wildcards))


rule merge_shared_peaks:
    input:
        locate_peak_file
    output:
        'input/bed/epigenome/{species}_{tissue}_{mark}_{sample}.peaks.bed'
    shell:
        'bedtools merge -d 100 -i {input} > {output}'


rule reduce_signal_to_peaks:
    input:
        peaks = 'input/bed/epigenome/{species}_{tissue}_{mark}_{sample}.peaks.bed',
        signal = 'input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.bg'
    output:
        'input/bedgraph/epigenome/peak_signal/{species}_{tissue}_{mark}_{sample}.cov.bg'
    wildcard_constraints:
        mark = '(H3K4me3|H3K36me3|H3K27ac)',
        sample = 'ET\-[0-9A-Z]+'
    shell:
        'bedtools intersect -u -a {input.signal} -b {input.peaks} > {output}'