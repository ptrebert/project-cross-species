
ENSEMBL_REFERENCES_ASSEMBLIES = {
    'bonobo': 'fasta/pan_paniscus/dna/Pan_paniscus.panpan1.1.dna_sm.toplevel.fa.gz',
    'cat': 'fasta/felis_catus/dna/Felis_catus.Felis_catus_9.0.dna_sm.toplevel.fa.gz',
    'chicken': 'fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna_sm.toplevel.fa.gz',
    'chimpanzee': 'fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa.gz',
    'cow': 'fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa.gz',
    'dog': 'fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna_sm.toplevel.fa.gz',
    'goat': 'fasta/capra_hircus/dna/Capra_hircus.ARS1.dna_sm.toplevel.fa.gz',
    'gorilla': 'fasta/gorilla_gorilla/dna/Gorilla_gorilla.gorGor4.dna_sm.toplevel.fa.gz',
    'horse': 'fasta/equus_caballus/dna/Equus_caballus.EquCab3.0.dna_sm.toplevel.fa.gz',
    'human': 'fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz',
    'macaque': 'fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_8.0.1.dna_sm.toplevel.fa.gz',
    'mouse': 'fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz',
    'opossum': 'fasta/monodelphis_domestica/dna/Monodelphis_domestica.monDom5.dna_sm.toplevel.fa.gz',
    'orangutan': 'fasta/pongo_abelii/dna/Pongo_abelii.PPYG2.dna_sm.toplevel.fa.gz',
    'pig': 'fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz',
    'platypus': 'fasta/ornithorhynchus_anatinus/dna/Ornithorhynchus_anatinus.OANA5.dna_sm.toplevel.fa.gz',
    'rabbit': 'fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_sm.toplevel.fa.gz',
    'rat': 'fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna_sm.toplevel.fa.gz',
    'sheep': 'fasta/ovis_aries/dna/Ovis_aries.Oar_v3.1.dna_sm.toplevel.fa.gz'
}


rule download_genome_assemblies:
    input:
        'annotation/species/{species_name}.info'
    output:
        'references/assemblies/raw_download/{species_name}.fna.gz'
    run:
        base_url = config['ensembl_release_url']
        load_path = ENSEMBL_REFERENCES_ASSEMBLIES[wildcards.species_name]
        source_url = os.path.join(base_url, load_path)
        exec = 'wget --quiet -O {output} ' + source_url
        shell(exec)


rule filter_sort_genome_assemblies:
    input:
        'references/assemblies/raw_download/{species_name}.fna.gz'
    output:
        fasta = 'references/assemblies/whole-genome/{species_name}.wg.fa',
        table = 'references/chromosomes/whole-genome/{species_name}.wg.sizes'
    log: 'log/references/filter_sort/whole-genome/{species_name}.log'
    params:
        script_dir = config['script_dir'],
        log_config = config['script_log_config']
    run:
        exec = '{params.script_dir}/ensembl_references/process_ensembl_assembly.py'
        exec += ' --log-config {params.log_config}'
        exec += ' --fasta-in {input} --debug'
        exec += ' --fasta-out {output.fasta}'
        exec += ' --chromosome-sizes {output.table}'
        exec += ' --debug'
        exec += ' &> {log}'
        shell(exec)


rule build_bowtie_index:
    input:
        'references/assemblies/whole-genome/{species_name}.wg.fa'
    output:
        expand('references/indices/bowtie/{{species_name}}/{{species_name}}.{idxnum}.bt2',
                idxnum=[1, 2, 3, 4])
    log: 'log/references/indices/bowtie/{species_name}.log'
    benchmark: 'run/references/indices/bowtie/{species_name}.rsrc'
    threads: 16
    run:
        exec = 'bowtie2-build --threads {threads}'
        exec += ' {input}'
        exec += ' references/indices/bowtie/{wildcards.species_name}/{wildcards.species_name}'
        exec += ' &> {log}'
        shell(exec)

