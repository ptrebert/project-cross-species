

RNA_READS_SE_WILDCARDS = glob_wildcards('input/fastq/transcriptome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.single.ok')

RNA_READS_PE_WILDCARDS = glob_wildcards('input/fastq/transcriptome/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.paired1.ok')


def transcriptome_single_end_combinations(species, tissues, marks, samples, kmers):

    valid_combinations = []
    root_folder = 'input/fastq/transcriptome'
    existing_combinations = set()
    for root, dirs, files in os.walk(root_folder):
        if files:
            valid_files = list(filter(lambda x: x.endswith('.single.ok'), files))
            [existing_combinations.add(vf.rsplit('_', 1)[0]) for vf in valid_files]

    for spec, tissue, mark, sample in zip(species, tissues, marks, samples):
        s = spec[1]
        t = tissue[1]
        m = mark[1]
        assert m == 'rna', 'Invalid mark {} for transcriptome {}'.format(m, sample)
        p = sample[1]
        combo = '_'.join([s, t, m, p])
        if combo in existing_combinations:
            for _, k in kmers:
                valid_combinations.append({'species': s,
                                            'tissue': t,
                                            'mark': m,
                                            'sample': p,
                                            'kmer': k})
    return valid_combinations


def transcriptome_paired_end_combinations(species, tissues, marks, samples, kmers):

    valid_combinations = []
    root_folder = 'input/fastq/transcriptome'
    existing_combinations = set()
    for root, dirs, files in os.walk(root_folder):
        if files:
            valid_files = list(filter(lambda x: x.endswith('.paired1.ok'), files))
            [existing_combinations.add(vf.rsplit('_', 1)[0]) for vf in valid_files]

    for spec, tissue, mark, sample in zip(species, tissues, marks, samples):
        s = spec[1]
        t = tissue[1]
        m = mark[1]
        assert m == 'rna', 'Invalid mark {} for transcriptome {}'.format(m, sample)
        p = sample[1]
        combo = '_'.join([s, t, m, p])
        if combo in existing_combinations:
            for _, k in kmers:
                valid_combinations.append({'species': s,
                                            'tissue': t,
                                            'mark': m,
                                            'sample': p,
                                            'kmer': k})
    return valid_combinations


rule preprocess_transcriptomes_master:
    input:
        expand('references/indices/salmon/whole-genome/protein-coding/{species}.k{kmer}.idx/sa.bin',
                species=sorted(config['species']),
                kmer=config['rna_kmer_seeds']),

        expand('input/tabular/transcriptome/temp/single/{species}_{tissue}_{mark}_{sample}.k{kmer}/quant.genes.sf',
                transcriptome_single_end_combinations,
                species=RNA_READS_SE_WILDCARDS.species,
                tissue=RNA_READS_SE_WILDCARDS.tissue,
                mark=RNA_READS_SE_WILDCARDS.mark,
                sample=RNA_READS_SE_WILDCARDS.sample,
                kmer=config['rna_kmer_seeds']),

        expand('input/tabular/transcriptome/temp/paired/{species}_{tissue}_{mark}_{sample}.k{kmer}/quant.genes.sf',
                transcriptome_paired_end_combinations,
                species=RNA_READS_PE_WILDCARDS.species,
                tissue=RNA_READS_PE_WILDCARDS.tissue,
                mark=RNA_READS_PE_WILDCARDS.mark,
                sample=RNA_READS_PE_WILDCARDS.sample,
                kmer=config['rna_kmer_seeds'])


rule build_salmon_index:
    input:
        'references/gene-models/whole-genome/protein-coding/{species}.wg.transcripts.fa'
    output:
        'references/indices/salmon/whole-genome/protein-coding/{species}.k{kmer}.idx/sa.bin'
    log:
        'log/references/indices/salmon/whole-genome/protein-coding/{species}.k{kmer}.idx.log'
    benchmark:
        'run/references/indices/salmon/whole-genome/protein-coding/{species}.k{kmer}.idx.rsrc'
    threads: 8
    params:
        index_dir = lambda wildcards, output: os.path.dirname(output[0])
    run:
        exec = 'salmon index --type quasi'
        exec += ' --transcripts {input}'
        exec += ' --kmerLen {wildcards.kmer}'
        exec += ' --keepDuplicates'  # should not have any effect
        exec += ' --threads {threads}'
        exec += ' --perfectHash'
        exec += ' --index {params.index_dir}'
        exec += ' &> {log}'
        shell(exec)


#rule build_star_index:
#    input:
#        'references/gene-models/whole-genome/protein-coding/{species}.wg.transcripts.fa'
#    output:
#        'references/indices/star/whole-genome/protein-coding/{species}.k{kmer}.idx/SA'
#    log:
#        'log/references/indices/star/whole-genome/protein-coding/{species}.k{kmer}.idx.log'
#    benchmark:
#        'run/references/indices/star/whole-genome/protein-coding/{species}.k{kmer}.idx.rsrc'
#    threads: 8
#    params:
#        index_dir = lambda wildcards, output: os.path.dirname(output[0])
#    run:
#        exec = 'salmon index --type quasi'
#        exec += ' --transcripts {input}'
#        exec += ' --kmerLen {wildcards.kmer}'
#        exec += ' --keepDuplicates'  # should not have any effect
#        exec += ' --threads {threads}'
#        exec += ' --perfectHash'
#        exec += ' --index {params.index_dir}'
#        exec += ' &> {log}'
#        shell(exec)


def collect_single_end_transcriptomes(wildcards):

    def filter_fun(filename):
        keep = wildcards.species in filename
        keep &= wildcards.tissue in filename
        keep &= wildcards.mark in filename
        keep &= wildcards.sample in filename
        keep &= filename.endswith('.single.ok')
        return keep

    transcriptome_path = 'input/fastq/transcriptome'
    process_fastq = []
    for root, dirs, files in os.walk(transcriptome_path):
        if not files:
            continue
        fastq_files = list(filter(filter_fun, files))
        fastq_files = [os.path.join(root, f.replace('.ok', '.fastq.gz')) for f in fastq_files]
        process_fastq.extend(fastq_files)
    return sorted(process_fastq)


rule quantify_single_end_transcriptome:
    input:
        fastq = collect_single_end_transcriptomes,
        index = 'references/indices/salmon/whole-genome/protein-coding/{species}.k{kmer}.idx/sa.bin',
        gtmap = 'references/gene-models/whole-genome/protein-coding/{species}.wg.gt-map.tsv'
    output:
        quant = 'input/tabular/transcriptome/temp/single/{species}_{tissue}_{mark}_{sample}.k{kmer}/quant.genes.sf',
        metadata = 'input/tabular/transcriptome/temp/single/{species}_{tissue}_{mark}_{sample}.k{kmer}/aux_info/meta_info.json',
    log:
        salmon = 'log/input/tabular/transcriptome/temp/{species}_{tissue}_{mark}_{sample}.k{kmer}.se.salmon.log',
    benchmark: 'run/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.k{kmer}.se.quant.rsrc'
    wildcard_constraints:
        species = '[a-z]+'
    threads: 16
    params:
        index_loc = lambda wildcards, input: os.path.dirname(input.index),
        output_dir = lambda wildcards, output: os.path.dirname(output.quant)
    run:
        exec = 'salmon quant --libType A'
        exec += ' --index {params.index_loc}'
        exec += ' --unmatedReads {input.fastq}'
        exec += ' --validateMappings'
        exec += ' --threads {threads}'
        exec += ' --seqBias --gcBias'
        exec += ' --geneMap {input.gtmap}'
        exec += ' --forgettingFactor 0.8'
        exec += ' --output {params.output_dir}'
        exec += ' &> {log}'
        shell(exec)


def collect_paired_end_transcriptomes_mate1(wildcards):

    def filter_fun(filename):
        keep = wildcards.species in filename
        keep &= wildcards.tissue in filename
        keep &= wildcards.mark in filename
        keep &= wildcards.sample in filename
        keep &= filename.endswith('.paired1.ok')
        return keep

    transcriptome_path = 'input/fastq/transcriptome'
    process_fastq_pair1 = []
    for root, dirs, files in os.walk(transcriptome_path):
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


def collect_paired_end_transcriptomes_mate2(wildcards):

    def filter_fun(filename):
        keep = wildcards.species in filename
        keep &= wildcards.tissue in filename
        keep &= wildcards.mark in filename
        keep &= wildcards.sample in filename
        keep &= filename.endswith('.paired2.ok')
        return keep

    transcriptome_path = 'input/fastq/transcriptome'
    process_fastq_pair2 = []
    for root, dirs, files in os.walk(transcriptome_path):
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


rule quantify_paired_end_transcriptome:
    input:
        fastq1 = collect_paired_end_transcriptomes_mate1,
        fastq2 = collect_paired_end_transcriptomes_mate2,
        index = 'references/indices/salmon/whole-genome/protein-coding/{species}.k{kmer}.idx/sa.bin',
        gtmap = 'references/gene-models/whole-genome/protein-coding/{species}.wg.gt-map.tsv'
    output:
        quant = 'input/tabular/transcriptome/temp/paired/{species}_{tissue}_{mark}_{sample}.k{kmer}/quant.genes.sf',
        metadata = 'input/tabular/transcriptome/temp/paired/{species}_{tissue}_{mark}_{sample}.k{kmer}/aux_info/meta_info.json',
    log:
        salmon = 'log/input/tabular/transcriptome/temp/{species}_{tissue}_{mark}_{sample}.k{kmer}.pe.salmon.log',
    benchmark: 'run/bam/epigenome/temp/{species}_{tissue}_{mark}_{sample}.k{kmer}.pe.quant.rsrc'
    wildcard_constraints:
        species = '[a-z]+'
    threads: 16
    params:
        index_loc = lambda wildcards, input: os.path.dirname(input.index),
        output_dir = lambda wildcards, output: os.path.dirname(output.quant)
    run:
        assert len(input.fastq1) == len(input.fastq2), \
            'Missing read files: {} v {}'.format(input.fastq1, input.fastq2)
        exec = 'salmon quant --libType A'
        exec += ' --index {params.index_loc}'
        exec += ' --mates1 {input.fastq1}'
        exec += ' --mates2 {input.fastq2}'
        exec += ' --validateMappings'
        exec += ' --threads {threads}'
        exec += ' --seqBias --gcBias'
        exec += ' --geneMap {input.gtmap}'
        exec += ' --forgettingFactor 0.8'
        exec += ' --output {params.output_dir}'
        exec += ' &> {log}'
        shell(exec)
