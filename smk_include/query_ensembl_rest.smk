
localrules: query_ensembl_species_list

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
        exec += ' 2> {log}'
        shell(exec)