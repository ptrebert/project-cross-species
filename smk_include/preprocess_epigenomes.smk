
localrules: compute_insert_size_estimates

SINGLE_EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.single.ok')
PAIRED_EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.paired1.ok')

EPIGENOME_WILDCARDS = glob_wildcards('input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.filt.bam')

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
                sample=EPIGENOME_WILDCARDS.sample)


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
        index_loc = lambda wildcards: 'references/indices/bowtie/{}/{}'.format(wildcards.species, wildcards.species)
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
        exec += ' -b -q 5'  # output BAM, MAPQ >= 5
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
        index_loc = lambda wildcards: 'references/indices/bowtie/{}/{}'.format(wildcards.species, wildcards.species)
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
        exec += ' -b -q 5'  # output BAM, MAPQ >= 5
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

    gsize = 'references/assemblies/whole-genome/{species}.k{kmer}.stats'

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
