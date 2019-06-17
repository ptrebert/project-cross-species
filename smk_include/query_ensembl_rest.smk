
localrules: query_ensembl_species_list, process_ensembl_species_list

rule query_ensembl_species_list:
    output:
        'annotation/ensembl_json_dump/ensembl_species.json'
    log:
        'log/ensembl_json_dump/ensembl_species.log'
    params:
        script_dir = config['script_dir'],
        log_config = config['script_log_config']
    run:
        exec = '{params.script_dir}/query_ensembl_rest.py'
        exec += ' --log-config {params.log_config}'
        exec += ' --query-type species'
        exec += ' --output {output}'
        exec += ' &> {log}'
        shell(exec)


rule process_ensembl_species_list:
    input:
        dump = 'annotation/ensembl_json_dump/ensembl_species.json',
        bioprojects = config['bioproject_table']
    output:
        #'annotation/species/ensembl_species.tsv',
        dynamic('annotation/species/{species}.info')
    params:
        script_dir = config['script_dir'],
        log_config = config['script_log_config']
    run:
        exec = '{params.script_dir}/ensembl_dumps/process_species_dump.py'
        exec += ' --log-config {params.log_config}'
        exec += ' --json-dump {input.dump}'
        exec += ' --species-table {input.bioprojects}'
        exec += ' --create-info-files'
        exec += ' --output {output}'
#        exec += ' &> {log}'
        shell(exec)
