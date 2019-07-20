
include: 'smk_include/aux_utilities.smk'
include: 'smk_include/query_ensembl_rest.smk'
include: 'smk_include/load_ensembl_references.smk'
include: 'smk_include/load_ena_data.smk'
include: 'smk_include/load_deep_data.smk'
include: 'smk_include/preprocess_epigenomes.smk'
include: 'smk_include/preprocess_transcriptomes.smk'

wildcard_constraints:
    species='[a-z]+',
    tissue='[a-z4]+'


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


rule master:
    input:
        rules.preprocess_transcriptomes_master.input,
        rules.preprocess_epigenomes_master.input,
        expand('references/alignments/whole-genome/{reference}_vs_{target}.lastz-net.tar.gz',
                reference='human',
                target=list(set(SPECIES) - set(['human']))),
        expand('references/alignments/whole-genome/{reference}_vs_{target}.lastz-net.tar.gz',
                reference='mouse',
                target=['chicken', 'cow', 'dog', 'opossum', 'pig', 'platypus', 'rat']),

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

    message: 'Executing MASTER rule'


