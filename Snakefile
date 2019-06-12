
include: 'smk_include/aux_utilities.smk'
include: 'smk_include/query_ensembl_rest.smk'

onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] SUCCESS - CSE" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
        shell('mail -s "[Snakemake] ERROR - CSE" {} < {{log}}'.format(config['notify_email']))

rule master:
    input:
        'annotation/ensembl_json_dumps/ensembl_species.json'

    message: 'Executing MASTER rule'

