



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