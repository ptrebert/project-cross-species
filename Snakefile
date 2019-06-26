
include: 'smk_include/aux_utilities.smk'
include: 'smk_include/query_ensembl_rest.smk'
include: 'smk_include/load_ensembl_references.smk'
include: 'smk_include/load_ena_data.smk'
include: 'smk_include/load_deep_data.smk'
include: 'smk_include/preprocess_epigenomes.smk'

onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] SUCCESS - CSE" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
        shell('mail -s "[Snakemake] ERROR - CSE" {} < {{log}}'.format(config['notify_email']))


SPECIES = sorted(ENSEMBL_REFERENCES_ASSEMBLIES.keys())

TRANSCRIPTOME_PROJECTS = [k for k, v in config['bioproject_label_map'].items() if v.startswith('T')]
EPIGENOME_PROJECTS = [k for k, v in config['bioproject_label_map'].items() if v.startswith('E')]

READSET_WILDCARDS = glob_wildcards('input/fastq/{datatype}/{bioproject}/{readset}.request')

def determine_readset_combinations(datatypes, bioprojects, readsets):
    valid_combinations = []
    root_folder = 'input/fastq'
    for datatype, bioproject, readset in zip(datatypes, bioprojects, readsets):
        dt = datatype[1]
        bp = bioproject[1]
        rs = readset[1]
        path = os.path.join(root_folder, dt, bp, rs + '.request')
        if os.path.isfile(path):
            valid_combinations.append({'datatype': dt,
                                        'bioproject': bp,
                                        'readset': rs})
    return valid_combinations

SINGLE_WILDCARDS = glob_wildcards('input/fastq/{datatype}/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.single.ok')

def determine_single_combinations(species, tissues, marks, samples):

    valid_combinations = []
    root_folder = 'input/fastq/epigenome'
    existing_combinations = set()
    for root, dirs, files in os.walk(root_folder):
        if files:
            valid_files = list(filter(lambda x: x.endswith('.single.ok'), files))
            [existing_combinations.add(vf.rsplit('_', 1)[0]) for vf in valid_files]

    for spec, tissue, mark, sample in zip(species, tissues, marks, samples):
        s = spec[1]
        t = tissue[1]
        m = mark[1]
        p = sample[1]
        combo = '_'.join([s, t, m, p])
        if combo in existing_combinations:
            valid_combinations.append({'species': s,
                                        'tissue': t,
                                        'mark': m,
                                        'sample': p})
    return valid_combinations


PAIRED_WILDCARDS = glob_wildcards('input/fastq/{datatype}/{bioproject}/{species}_{tissue}_{mark}_{sample}_{run}.paired1.ok')

def determine_paired_combinations(species, tissues, marks, samples):

    valid_combinations = []
    root_folder = 'input/fastq/epigenome'
    existing_combinations = set()
    for root, dirs, files in os.walk(root_folder):
        if files:
            valid_files = list(filter(lambda x: x.endswith('.paired1.ok'), files))
            [existing_combinations.add(vf.rsplit('_', 1)[0]) for vf in valid_files]

    for spec, tissue, mark, sample in zip(species, tissues, marks, samples):
        s = spec[1]
        t = tissue[1]
        m = mark[1]
        p = sample[1]
        combo = '_'.join([s, t, m, p])
        if combo in existing_combinations:
            valid_combinations.append({'species': s,
                                        'tissue': t,
                                        'mark': m,
                                        'sample': p})
    return valid_combinations


rule master:
    input:
        'input/checkpoints/EGAC00001000331.epigenome.folder',
        'input/checkpoints/EGAC00001000331.transcriptome.folder',
        expand('references/gene-models/whole-genome/protein-coding/{species}.wg.transcripts.fa',
                species=SPECIES),

        expand('input/fastq/{datatype}/EGAC00001000331.link.ok',
                datatype=['epigenome', 'transcriptome']),

        expand(rules.handle_deep_batch_download.output,
                datatype=['epigenome', 'transcriptome']),
        expand(rules.create_ena_request_files.output,
                datatype='epigenome',
                bioproject=EPIGENOME_PROJECTS),
        expand(rules.create_ena_request_files.output,
                datatype='transcriptome',
                bioproject=TRANSCRIPTOME_PROJECTS),
        expand('input/fastq/{datatype}/{bioproject}/{readset}.ok',
                determine_readset_combinations,
                datatype=READSET_WILDCARDS.datatype,
                bioproject=READSET_WILDCARDS.bioproject,
                readset=READSET_WILDCARDS.readset),

        expand('input/bam/epigenome/temp/single/{species}_{tissue}_{mark}_{sample}.filt.bam',
                determine_single_combinations,
                species=SINGLE_WILDCARDS.species,
                tissue=SINGLE_WILDCARDS.tissue,
                mark=SINGLE_WILDCARDS.mark,
                sample=SINGLE_WILDCARDS.sample),

        expand('input/bam/epigenome/temp/paired/{species}_{tissue}_{mark}_{sample}.filt.bam',
                determine_paired_combinations,
                species=PAIRED_WILDCARDS.species,
                tissue=PAIRED_WILDCARDS.tissue,
                mark=PAIRED_WILDCARDS.mark,
                sample=PAIRED_WILDCARDS.sample),

    message: 'Executing MASTER rule'


