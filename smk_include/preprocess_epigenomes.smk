

#def determine_single_combinations(species, tissues, marks, samples):
#
#    valid_combinations = []
#    root_folder = 'input/fastq/epigenome'
#    existing_combinations = set()
#    for root, dirs, files in os.walk(root_folder):
#        if files:
#            valid_files = list(filter(lambda x: x.endswith('.single.ok'), files))
#            [existing_combinations.add(vf.rsplit('_', 1)[0]) for vf in valid_files]
#
#    for spec, tissue, mark, sample in zip(species, tissues, marks, samples):
#        s = spec[1]
#        t = tissue[1]
#        m = mark[1]
#        p = sample[1]
#        combo = '_'.join([s, t, m, p])
#        if combo in existing_combinations:
#            valid_combinations.append({'species': s,
#                                        'tissue': t,
#                                        'mark': m,
#                                        'sample': p})
#    return valid_combinations
#
#
#def determine_paired_combinations(species, tissues, marks, samples):
#
#    valid_combinations = []
#    root_folder = 'input/fastq/epigenome'
#    existing_combinations = set()
#    for root, dirs, files in os.walk(root_folder):
#        if files:
#            valid_files = list(filter(lambda x: x.endswith('.paired1.ok'), files))
#            [existing_combinations.add(vf.rsplit('_', 1)[0]) for vf in valid_files]
#
#    for spec, tissue, mark, sample in zip(species, tissues, marks, samples):
#        s = spec[1]
#        t = tissue[1]
#        m = mark[1]
#        p = sample[1]
#        combo = '_'.join([s, t, m, p])
#        if combo in existing_combinations:
#            valid_combinations.append({'species': s,
#                                        'tissue': t,
#                                        'mark': m,
#                                        'sample': p})
#    return valid_combinations


SINGLE_EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.single.ok')
PAIRED_EPIGENOME_WILDCARDS = glob_wildcards('input/fastq/epigenome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.paired1.ok')

EPIGENOME_WILDCARDS = glob_wildcards('input/bam/epigenome/temp/{layout}/{species}_{tissue}_{mark}_{sample}.filt.bam')

rule preprocess_epigenomes_master:
    input:
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
                        idxnum=[1, 2, 3, 4])
    output:
        bam = 'input/bam/epigenome/temp/single/{species}_{tissue}_{mark}_{sample}.filt.bam',
        metrics = 'input/bam/epigenome/metrics/{species}_{tissue}_{mark}_{sample}.bowtie2.metrics'
    log:
        bowtie = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.bowtie.log',
        samtools = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.samtools.log'
    benchmark: 'run/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.aln.rsrc'
    wildcard_constraints:
        species = '[a-z]+'
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
                        idxnum=[1, 2, 3, 4])
    output:
        bam = 'input/bam/epigenome/temp/paired/{species}_{tissue}_{mark}_{sample}.filt.bam',
        metrics = 'input/bam/epigenome/metrics/{species}_{tissue}_{mark}_{sample}.bowtie2.metrics'
    log:
        bowtie = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.bowtie.log',
        samtools = 'log/input/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.samtools.log'
    benchmark: 'run/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.aln.rsrc'
    wildcard_constraints:
        species = '[a-z]+'
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


checkpoint collect_read_length_statistics:
    input:
        expand('input/bam/epigenome/temp/{layout}/{{filename}}.read-stats',
                layout=['single', 'paired'])
    output:
        directory('input/bam/epigenome/temp/{layout}/read-length')
    run:
        for stat_file in input:
            rlen = None
            with open(stat_file, 'r') as stats:
                for line in stats:
                    if line.startswith('RL'):
                        rlen = line.read().split()[1]
                        assert int(rlen), 'Non-integer read length for file {}: {}'.format(stat_file, rlen)
                        break
                    continue
            assert rlen is not None, 'No read length for file {}'.format(stat_file)
            output_path = os.path.join(output[0], 'read-length', '{wildcards.filename}.k{}'.format(rlen))
            with open(output_path, 'w') as dump:
                pass


def collect_signal_coverage_input(wildcards):

    layout_single = checkpoints.collect_read_length_statistics.get(layout='single').output[0]
    layout_paired = checkpoints.collect_read_length_statistics.get(layout='paired').output[0]

    file_id = '_'.join([wildcards.species, wildcards.tissue, wildcards.mark, wildcards.sample])

    single_bam = os.path.join(os.path.split(layout_single)[0], file_id + '.filt.bam')
    paired_bam = os.path.join(os.path.split(layout_paired)[0], file_id + '.filt.bam')

    if os.path.isfile(single_bam):
        search_path = layout_single
        input_bam = single_bam
    elif os.path.isfile(paired_bam):
        search_path = layout_paired
        input_bam = paired_bam
    else:
        return dict()

    rlen_file = list(filter(lambda x: file_id in x, os.listdir(search_path)))
    if len(rlen_file) < 1:
        raise ValueError('No read length file for file {}'.format(input_bam))

    rlen_path = os.path.join(search_path, rlen_file)
    kmer = rlen_file.split('.')[-1].strip('k')

    gstat_file = 'references/assemblies/whole-genome/{}.k{}.stats'.format(wildcards.species, kmer)

    assert os.path.isfile(rlen_path), 'rlen_path invalid'
    assert os.path.isfile(gstat_file), 'gstat_file invalid'

    file_combo = {'bam': input_bam,
                  'genome_size': gstat_file,
                  'kmer_size': rlen_path}

    return file_combo


rule compute_signal_coverage:
    input:
        unpack(collect_signal_coverage_input)
    output:
        'input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.bg'
    threads: 6
    log:
        'log/input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.bg'
    benchmark:
        'run/input/bedgraph/epigenome/{species}_{tissue}_{mark}_{sample}.cov.bg'
    run:
        effective_genome_size = None
        with open(input.genome_size, 'r') as counts:
            for line in counts:
                if line.startswith('number of unique k-mers'):
                    effective_genome_size = line.split()[1]
                    assert int(effective_genome_size), 'Non-integer genome size: {}'.format(input.genome_size)
                    break
        assert effective_genome_size is not None, 'No genome size for file: {}'.format(input.genome_size)

        exec = 'bamCoverage --bam {input.bam}'
        exec += ' --outFileName {output}'
        exec += ' --outFileFormat bedgraph'
        exec += ' --binSize 25 --numberOfProcessors {threads}'
        exec += ' --effectiveGenomeSize ' + effective_genome_size
        exec += ' --normalizeUsing RPGC' # according to docs, should be 1x normalization
        exec += ' --verbose &> {log}'
        shell(exec)
