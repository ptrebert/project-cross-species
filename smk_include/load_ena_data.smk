
localrules: load_bioproject_file_reports, parse_ena_file_report, create_ena_request_files, check_file_integrity

rule load_bioproject_file_reports:
    output:
        'annotation/ena_project_reports/{bioproject}.tsv'
    run:
        load_url = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession='
        load_url += '{accession}'
        load_url += '&result=read_run&fields='
        load_url += 'study_accession,sample_accession,secondary_sample_accession'
        load_url += ',experiment_accession,run_accession,submission_accession'
        load_url += ',tax_id,scientific_name,instrument_platform,instrument_model'
        load_url += ',library_name,library_layout,library_strategy,library_source'
        load_url += ',library_selection,read_count,center_name,study_title,fastq_md5'
        load_url += ',study_alias,experiment_alias,run_alias,fastq_ftp,submitted_ftp'
        load_url += ',sample_alias,sample_title'
        load_url += '&download=txt'

        tmp = load_url.format(**{'accession': bp})
        if os.path.isfile(output[0]):
            exec = 'touch {output}'
        else:
            exec = 'wget --quiet -O {{output}} "{}"'.format(tmp)
        shell(exec)


rule parse_ena_file_report:
    input:
        species_table = 'annotation/species/ensembl_species.tsv',
        file_report = 'annotation/ena_project_reports/{bioproject}.tsv'
    output:
        'input/checkpoints/{bioproject}.{datatype}.download'
    run:
        bioproject_parser = get_file_report_parser(wildcards.bioproject)
        existing_downloads = set()
        if os.path.isfile(output[0]):
            with open(output[0], 'r') as listing:
                [existing_downloads.add(tuple(record.strip().split('\t'))) for record in listing]

        requested_downloads = bioproject_parser(input.species_table, input.file_report)

        assert len(requested_downloads) == len(set(requested_downloads)), \
            'Duplicates read for project: {}'.format(input.file_report)

        rewrite = False
        if any([r not in existing_downloads for r in requested_downloads]):
            rewrite = True
        if rewrite:
            with open(output[0], 'w') as listing:
                for req in sorted(requested_downloads):
                    # write triples:
                    # local file name - remote URL - remote_md5
                    _ = listing.write('\t'.join(req) + '\n')
        # end of rule


rule create_ena_request_files:
    input:
        'input/checkpoints/{bioproject}.{datatype}.download'
    output:
        request_files = touch('input/fastq/{datatype}/{bioproject}/{bioproject}.chk')
        #request_files = directory('input/fastq/{datatype}/{bioproject}')
    wildcard_constraints:
        datatype = '(transcriptome|epigenome)',
        bioproject = 'PR[A-Z0-9]+'
    run:
        #output_dir = output.request_files
        output_dir = os.path.dirname(output.request_files)
        os.makedirs(output_dir, exist_ok=True)
        assert os.path.isdir(output_dir), 'No output folder created: {}'.format(output_dir)
        with open(input[0], 'r') as listing:
            for record in listing:
                local_file, remote_path, remote_md5 = record.strip().split('\t')
                req_file = local_file.replace('.fastq.gz', '.request')
                req_file_path = os.path.join(output_dir, req_file)
                local_file_path = os.path.join(output_dir, local_file)
                with open(req_file_path, 'w') as dump:
                    _ = dump.write('\t'.join([local_file_path, remote_path, remote_md5]) + '\n')
        # end of rule


rule handle_request_file:
    input:
        rules.create_ena_request_files.output[0],
        request_file = 'input/fastq/{datatype}/{bioproject}/{readset}.request',
    output:
        'input/fastq/{datatype}/{bioproject}/{readset}.fastq.gz'
    log:
        'log/input/fastq/{datatype}/{bioproject}/{readset}.log'
    run:
        with open(input.request_file, 'r') as info:
            local_path, remote_path, _ = info.read().strip().split('\t')
        if os.path.isfile(local_path):
            # since rule is re-executed, just touch the file
            exec = 'touch {}'.format(local_path)
        else:
            if remote_path.startswith('ftp.'):
                remote_path = 'ftp://' + remote_path
            else:
                raise ValueError('Unsupported transfer protocol: {}'.format(input.request_file))
            exec = 'wget --timeout=30 --tries=2 -O {} "{}" &> {{log}}'.format(local_path, remote_path)
        shell(exec)


rule compute_md5_checksum:
    input:
        'input/fastq/{datatype}/{bioproject}/{readset}.fastq.gz'
    output:
        'input/fastq/{datatype}/{bioproject}/{readset}.md5'
    shell:
        'md5sum {input} > {output}'


rule check_file_integrity:
    input:
        request = 'input/fastq/{datatype}/{bioproject}/{readset}.request',
        md5 = 'input/fastq/{datatype}/{bioproject}/{readset}.md5'
    output:
        'input/fastq/{datatype}/{bioproject}/{readset}.ok'
    run:
        with open(input.request, 'r') as info:
            local_path, _, remote_md5 = info.read().strip().split('\t')
        with open(input.md5, 'r') as info:
            local_md5, _ = info.read().strip().split()
        if not local_md5 == remote_md5:
            os.unlink(local_path)
        else:
            with open(output[0], 'w') as valid:
                pass
    # end of rule


def get_file_report_parser(bioproject):

    mapping = {
        'PRJEB6906': parser_villar_epigenomes,
        'PRJNA143627': parser_brawand_transcriptomes,
        'PRJNA164819': parser_broad_transcriptomes,
        'PRJNA184055': parser_fushan_transcriptomes,
        'PRJNA244152': parser_jiang_transcriptomes,
        'PRJNA30709': parser_human_encode,
        'PRJNA371787': parser_mouse_blueprint,
        'PRJNA373805': parser_mouse_blueprint,
        'PRJNA63443': parser_human_encode,
        'PRJNA63471': parser_mouse_encode,
        'PRJNA66167': parser_mouse_encode
        }

    parser_function = mapping.get(bioproject, None)
    if parser_function is None:
        raise ValueError('No parser defined for bioproject {}'.format(bioproject))

    return parser_function


def make_local_file_name(species, tissue, mark, accessions, bioproject):

    bp_label = config['bioproject_label_map'][bioproject]

    filename = '{}_{}_{}_{}-{}_{}'.format(species,
                                          tissue,
                                          mark,
                                          bp_label,
                                          accessions.sample_accession,
                                          accessions.run_accession
                                          )
    return filename


def combine_local_remote(local_name, record):

    is_paired = ';' in record.fastq_ftp
    load_info = []
    for url, md5 in zip(record.fastq_ftp.split(';'), record.fastq_md5.split(';')):
        if is_paired:
            mate_num = url.split('_')[-1].strip('.fastq.gz')
            assert int(mate_num), 'Unexpected mate number: {}'.format(record)
            tmp = local_name + '.paired{}.fastq.gz'.format(mate_num)
            if mate_num == '1' and 'mouse_ncd4_Input_ET-SAMN06312021_SRR5237812' in local_name:
                # for the ncd4 H3K27ac sample, create a matched Input
                # by selecting mate 1 reads of the paired Input
                tmp2 = local_name.replace('ET-SAMN06312021', 'ET-SAMN06312017') + '.single.fastq.gz'
                load_info.append((tmp2, url, md5))
            if mate_num == '1' and 'mouse_ncd4_Input_ET-SAMN06312022_SRR5237810' in local_name:
                tmp2 = local_name.replace('ET-SAMN06312022', 'ET-SAMN06312018') + '.single.fastq.gz'
                load_info.append((tmp2, url, md5))
        else:
            tmp = local_name + '.single.fastq.gz'
        load_info.append((tmp, url, md5))
    return load_info


def parser_mouse_encode(species, file_report):

    import pandas as pd
    import collections as col

    mouse_tissues = ['Liver', 'Heart', 'Kidney', 'ES-E14', 'stem cell line E14']
    histone_marks = config['histone_marks']
    exclude_age = [
        'E0', '10.5', '11.5', '12.5',
        '13.5', '14.5', '15.5', '16.5',
        'superseded', '_0', '0 day'
    ]

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    table = pd.read_csv(file_report, sep='\t')
    for row in table.itertuples():
        encode_sample = row.sample_title
        species = 'mouse'
        tissue, mark = None, None
        for t in mouse_tissues:
            if t in encode_sample:
                tissue = t.lower()
                if tissue == 'es-e14' or tissue == 'stem cell line e14':
                    tissue = 'esc'
                break
        if tissue is None:
            continue
        for m in histone_marks:
            if m in encode_sample:
                mark = m
        if mark is None:
            if 'rna-seq' in row.library_strategy.lower():
                mark = 'rna'
                if 'single' in row.library_layout.lower():
                    continue
            else:
                continue
        assert mark is not None, 'No target identified: {}'.format(row)

        if tissue != 'esc':
            if any([ea in encode_sample for ea in exclude_age]):
                continue

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )
        if 'SAMN01174211' in local_name and mark == 'H3K4me3':
            # Mouse Encode Chip-Seq Stanford experiment
            # has no K36me3 and K27ac, so pointless of
            # including only K4me3
            continue
        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    if 'rna' in selected_data_files[0][0]:
        # use only samples that were replicated at least once
        rep_counter = col.Counter()
        for name, _, _ in selected_data_files:
            if 'paired1' in name:
                sample = name.split('_')[3]
                rep_counter[sample] += 1

        rep_samples = []
        for local, remote, md5 in selected_data_files:
            sample = local.split('_')[3]
            if rep_counter[sample] > 1:
                rep_samples.append((local, remote, md5))
        selected_data_files = rep_samples

    return selected_data_files


def parser_villar_epigenomes(species, file_report):

    import pandas as pd

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    for row in metadata.itertuples():
        assert row.library_layout == 'SINGLE', 'Library could be paired: {}'.format(row)
        row_taxon = row.tax_id
        if row_taxon not in species_table['taxon_id'].values:
            if row_taxon == 10092:
                # for whatever reason, the species taxon id
                # is given for house mouse [10092] instead
                # of the generic mus musculus [10090]
                row_taxon = 10090
            else:
                continue

        species = species_table.loc[species_table['taxon_id'] == row_taxon, 'name'].values
        assert len(species) == 1, 'No species selected: {}'.format(row)
        species = species[0]
        tissue = 'liver'
        mark = 'none'
        if 'H3K27Ac' in row.submitted_ftp:
            mark = 'H3K27ac'
        elif 'H3K4me3' in row.submitted_ftp:
            mark = 'H3K4me3'
        elif '_input_' in row.submitted_ftp or '_Input_' in row.submitted_ftp:
            mark = 'Input'
        else:
            # For 4 files, annotation is incomplete/wrong:
            # do2602_unk_liver_unk_rn5Rattus_norvegicusRat9_SAN01.fq.gz
            # do2601_unk_liver_unk_rn5Rattus_norvegicusRat7_SAN01.fq.gz
            # do2603_unk_liver_unk_rn5Rattus_norvegicusRat10_SAN01.fq.gz
            # do2592_unk_liver_unk_canFam2Canis_familiarisDog5_SAN01.fq.gz
            # Based available data, these seem to be Input samples
            mark = 'Input'

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    return selected_data_files


def parser_brawand_transcriptomes(species, file_report):

    import pandas as pd

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    for row in metadata.itertuples():
        if row.tax_id not in species_table['taxon_id'].values:
            continue
        tissue = ''
        if ' lv ' in row.sample_title:
            tissue = 'liver'
        elif ' ht ' in row.sample_title:
            tissue = 'heart'
        elif ' kd ' in row.sample_title:
            tissue = 'kidney'
        else:
            continue
        assert row.library_layout == 'SINGLE', 'Library could be paired: {}'.format(row)
        species = species_table.loc[species_table['taxon_id'] == row.tax_id, 'name'].values
        assert len(species) == 1, 'No species selected: {}'.format(row)
        species = species[0]
        mark = 'rna'

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    return selected_data_files


def parser_broad_transcriptomes(species, file_report):

    import pandas as pd

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]
    for row in metadata.itertuples():
        species = 'opossum'
        mark = 'rna'
        tissue = ''
        if '_Liver' in row.sample_alias:
            tissue = 'liver'
        elif '_Kidney' in row.sample_alias:
            tissue = 'kidney'
        elif '_Blood' in row.sample_alias:
            tissue = 'blood'
        elif '_Heart' in row.sample_alias:
            tissue = 'heart'
        else:
            continue
        assert row.library_layout == 'PAIRED', 'Library could be single: {}'.format(row)

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    return selected_data_files


def parser_fushan_transcriptomes(species, file_report):

    import pandas as pd

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    for row in metadata.itertuples():
        assert row.library_layout == 'PAIRED', 'Library could be single: {}'.format(row)
        row_taxon = row.tax_id
        if row_taxon not in species_table['taxon_id'].values:
            continue

        species = species_table.loc[species_table['taxon_id'] == row_taxon, 'name'].values
        assert len(species) == 1, 'No species selected: {}'.format(row)
        species = species[0]
        tissue = None
        mark = 'rna'
        if '.lv.' in row.sample_title:
            tissue = 'liver'
        elif '.kd.' in row.sample_title:
            tissue = 'kidney'
        else:
            continue
        assert tissue is not None, 'No tissue detected: {}'.format(row)

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    return selected_data_files


def parser_jiang_transcriptomes(species, file_report):

    import pandas as pd

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    for row in metadata.itertuples():
        species = 'sheep'
        mark = 'rna'
        tissue = ''
        if 'Liver' in row.sample_title:
            tissue = 'liver'
        elif 'Kidney' in row.sample_title:
            tissue = 'kidney'
        elif 'Heart' in row.sample_title:
            tissue = 'heart'
        else:
            continue
        assert row.library_layout == 'PAIRED', 'Library could be single: {}'.format(row)

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    return selected_data_files


def parser_mouse_blueprint(species, file_report):

    import pandas as pd

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    for row in metadata.itertuples():
        species = 'mouse'
        mark = None
        tissue = 'ncd4'
        if 'NaiveCD4' not in row.sample_title:
            continue
        if 'H3K4me3' in row.sample_title:
            mark = 'H3K4me3'
        elif 'H3K27ac' in row.sample_title:
            mark = 'H3K27ac'
        elif 'H3K36me3' in row.sample_title:
            mark = 'H3K36me3'
        elif 'INPUT' in row.sample_title:
            mark = 'Input'
        elif row.library_strategy == 'RNA-Seq':
            mark = 'rna'
        else:
            continue
        assert mark is not None, 'No target detected: {}'.format(row)

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    return selected_data_files


def parser_human_encode(species, file_report):

    import pandas as pd
    import collections as col

    species_table = pd.read_csv(species, sep='\t')
    metadata = pd.read_csv(file_report, sep='\t')

    selected_data_files = []
    bioproject = os.path.basename(file_report).split('.')[0]

    for row in metadata.itertuples():
        species = 'human'
        mark = None
        tissue = None
        if row.library_strategy not in ['ChIP-Seq', 'RNA-Seq']:
            continue
        if 'H1-hESC' not in row.sample_title:
            continue
        tissue = 'esc'

        if 'H3K4me3' in row.sample_title:
            mark = 'H3K4me3'
        elif 'H3K27ac' in row.sample_title:
            mark = 'H3K27ac'
        elif 'H3K36me3' in row.sample_title:
            mark = 'H3K36me3'
        elif 'Input' in row.sample_title or 'Control' in row.sample_title:
            mark = 'Input'
        elif row.library_strategy == 'RNA-Seq':
            mark = 'rna'
        else:
            continue
        assert mark is not None, 'No target detected: {}'.format(row)

        if mark == 'rna':
            if row.library_strategy == 'SINGLE':
                continue
            if 'cytosol' in row.library_name or 'nucleus' in row.library_name:
                continue

        local_name = make_local_file_name(species,
                                          tissue,
                                          mark,
                                          row,
                                          bioproject
                                          )

        load_info = combine_local_remote(local_name, row)
        selected_data_files.extend(load_info)

    if 'rna' in selected_data_files[0][0]:
        # use only samples that were replicated at least once
        rep_counter = col.Counter()
        for name, _, _ in selected_data_files:
            if 'paired1' in name:
                sample = name.split('_')[3]
                rep_counter[sample] += 1

        rep_samples = []
        for local, remote, md5 in selected_data_files:
            sample = local.split('_')[3]
            if rep_counter[sample] > 1:
                rep_samples.append((local, remote, md5))
        selected_data_files = rep_samples

    return selected_data_files
