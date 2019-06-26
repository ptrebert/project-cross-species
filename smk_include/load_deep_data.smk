
localrules: create_deep_checkpoint_file, link_deep_data

DEEP_ASPERA_SOURCES = [
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf02_LiHe_Ct > replicate1-H3K27ac > paired > run131023_SN7001180_0077_AC2TKLACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf02_LiHe_Ct > replicate1-H3K36me3 > paired > run131023_SN7001180_0077_AC2TKLACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf02_LiHe_Ct > replicate1-H3K4me3 > paired > run131023_SN7001180_0077_AC2TKLACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf02_LiHe_Ct > replicate1-Input > paired > run131023_SN7001180_0077_AC2TKLACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf03_LiHe_Ct > replicate1-H3K27ac > paired > run140918_SN7001180_0146_C4CFEACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf03_LiHe_Ct > replicate1-H3K36me3 > paired > run140918_SN7001180_0146_C4CFEACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf03_LiHe_Ct > replicate1-H3K4me3 > paired > run140918_SN7001180_0146_C4CFEACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hf03_LiHe_Ct > replicate1-Input > paired > run140918_SN7001180_0146_C4CFEACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm09_LiHe_Ct > replicate1-H3K27ac > paired > run160805_SN7001180_0278_BC9LH5ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm09_LiHe_Ct > replicate1-H3K36me3 > paired > run160805_SN7001180_0278_BC9LH5ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm09_LiHe_Ct > replicate1-H3K4me3 > paired > run160805_SN7001180_0278_BC9LH5ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm09_LiHe_Ct > replicate1-Input > paired > run160805_SN7001180_0278_BC9LH5ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm16_LiHe_Ct > replicate1-H3K27ac > paired > run160719_SN7001180_0276_BC975GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm16_LiHe_Ct > replicate1-H3K36me3 > paired > run160719_SN7001180_0276_BC975GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm16_LiHe_Ct > replicate1-H3K4me3 > paired > run160719_SN7001180_0276_BC975GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 41_Hm16_LiHe_Ct > replicate1-Input > paired > run160719_SN7001180_0276_BC975GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf03_BlTN_Ct > replicate1-H3K27ac > paired > run141031_SN106_0104_A_C376GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf03_BlTN_Ct > replicate1-H3K36me3 > paired > run141031_SN106_0104_A_C376GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf03_BlTN_Ct > replicate1-H3K4me3 > paired > run141031_SN106_0104_A_C376GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf03_BlTN_Ct > replicate1-Input > paired > run141031_SN106_0104_A_C376GACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf04_BlTN_Ct > replicate1-H3K27ac > paired > run150327_SN106_0122_A_C7228ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf04_BlTN_Ct > replicate1-H3K36me3 > paired > run150327_SN106_0122_A_C7228ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf04_BlTN_Ct > replicate1-H3K4me3 > paired > run150327_SN106_0122_A_C7228ACXX > sequence",
"download > sequencing > chip_seq_sequencing > view-by-pid > 51_Hf04_BlTN_Ct > replicate1-Input > paired > run150327_SN106_0122_A_C7228ACXX > sequence",
"download > sequencing > strand_specific_mrna_sequencing > view-by-pid > 41_Hf02_LiHe_Ct > replicate1 > paired > run140212_SN758_0152_BC3ETGACXX > sequence",
"download > sequencing > strand_specific_mrna_sequencing > view-by-pid > 41_Hf03_LiHe_Ct > replicate1 > paired > run140827_SN751_0200_BC533PACXX > sequence",
"download > sequencing > strand_specific_mrna_sequencing > view-by-pid > 41_Hm09_LiHe_Ct > replicate1 > paired > run161128_D00695_0001_ACADVUANXX > sequence",
"download > sequencing > strand_specific_mrna_sequencing > view-by-pid > 41_Hm16_LiHe_Ct > replicate1 > paired > run160609_D00695_0055_AC8TCWANXX > sequence",
"download > sequencing > strand_specific_mrna_sequencing > view-by-pid > 51_Hf03_BlTN_Ct > replicate1 > paired > run141203_SN471_0187_B_C5T86ACXX > sequence",
"download > sequencing > strand_specific_mrna_sequencing > view-by-pid > 51_Hf04_BlTN_Ct > replicate1 > paired > run150402_SN541_0306_B_C6B4MACX > sequence"
]



rule create_deep_checkpoint_file:
    output:
        epigenome = 'input/checkpoints/EGAC00001000331.epigenome.folder',
        transcriptome = 'input/checkpoints/EGAC00001000331.transcriptome.folder'
    run:
        sample_map = {
            '41_Hf02_LiHe_Ct': 'SAM41FM0209CT',
            '41_Hf03_LiHe_Ct': 'SAM41FM0316CT',
            '41_Hm09_LiHe_Ct': 'SAM41FM0209CT',
            '41_Hm16_LiHe_Ct': 'SAM41FM0316CT',
            '51_Hf03_BlTN_Ct': 'SAM51FF0304CT',
            '51_Hf04_BlTN_Ct': 'SAM51FF0304CT'
            }
        epigenome_folders = []
        transcriptome_folders = []
        for source in DEEP_ASPERA_SOURCES:
            components = source.split(' > ')
            if 'chip_seq' in components[2]:
                datatype = 'epigenome'
            elif 'mrna' in components[2]:
                datatype = 'transcriptome'
            else:
                raise ValueError('Unexpected data type: {}'.format(components[2]))
            sample = sample_map[components[4]]
            if datatype == 'epigenome':
                sample = 'ET-' + sample
            else:
                sample = 'TT-' + sample
            run = components[7].split('_', 1)[1].replace('_', '')
            mark = None
            if datatype == 'transcriptome':
                mark = 'rna'
            else:
                for m in config['histone_marks']:
                    if m in components[5]:
                        mark = m
                        break
            assert mark is not None, 'No experiment target detected: {}'.format(source)
            if 'BlTN' in components[4]:
                tissue = 'ncd4'
            elif 'LiHe' in components[4]:
                tissue = 'hepa'
            else:
                raise ValueError('Unexpected tissue: {}'.format(source))
            local_name = '_'.join(['human', tissue, mark, sample, run])
            local_folder = os.path.join('input', 'batch', datatype, local_name)
            remote_path = '/'.join(components)
            if datatype == 'epigenome':
                epigenome_folders.append((local_folder, remote_path))
            else:
                transcriptome_folders.append((local_folder, remote_path))

        with open(output.epigenome, 'w') as dump:
            for record in sorted(epigenome_folders):
                _ = dump.write('\t'.join(record) + '\n')

        with open(output.transcriptome, 'w') as dump:
            for record in sorted(transcriptome_folders):
                _ = dump.write('\t'.join(record) + '\n')
    # end of rule


rule handle_deep_batch_download:
    input:
        batch_file = 'input/checkpoints/EGAC00001000331.{datatype}.folder',
        aspera_credentials = config['aspera_credentials']
    output:
        touch('input/batch/{datatype}/EGAC00001000331.batch.ok')
    log: 'log/input/batch/deep_download_{datatype}.log'
    run:
        with open(input.aspera_credentials, 'r') as cred:
            user = cred.readline().strip()
            pw = cred.readline().strip()

        log_dir = 'log/input/batch/' + wildcards.datatype
        os.makedirs(log_dir, exist_ok=True)

        with open(input.batch_file, 'r') as batch:
            for line in batch:
                local_path, remote_path = line.strip().split()
                os.makedirs(local_path, exist_ok=True)
                exec = 'ASPERA_SCP_PASS="' + pw + '"'
                exec += ' ascp -q -T -k1 -Q -P 33001 -l800M'
                exec += ' --overwrite=diff --symbolic-links=follow'
                exec += ' --src-base=' + remote_path
                exec += ' -L ' + log_dir
                exec += ' ' + user + '@dkfzaspera.dkfz-heidelberg.de:' + remote_path
                exec += ' ' + local_path
                exec += ' &> {log}'
                shell(exec)


rule link_deep_data:
    input:
        check_file = 'input/batch/{datatype}/EGAC00001000331.batch.ok',
        batch_file = 'input/checkpoints/EGAC00001000331.{datatype}.folder'
    output:
        touch('input/fastq/{datatype}/EGAC00001000331.link.ok')
    log: 'log/input/fastq/{datatype}/EGAC00001000331.link.log'
    run:
        import re

        output_dir = output[0].split('.')[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(log[0], 'w') as logfile:
            with open(input.batch_file, 'r') as listing:
                for line in listing:
                    local_path, _ = line.strip().split()
                    sample_prefix = os.path.split(local_path)[1]
                    for fastq in os.listdir(local_path):
                        src = os.path.join(local_path, fastq)
                        mobj = re.search('L00[1-9]', fastq)
                        if mobj is None:
                            raise ValueError('No lane identifier: {}'.format(fastq))
                        lane = mobj.group(0)
                        if '_R1' in fastq:
                            mate = 1
                        elif '_R2' in fastq:
                            mate = 2
                        else:
                            raise ValueError('Unexpected mate number: {}'.format(src))
                        output_fastq = sample_prefix + lane + '.paired{}.fastq.gz'.format(mate)
                        output_valid = sample_prefix + lane + '.paired{}.ok'.format(mate)

                        dst = os.path.join(output_dir, output_fastq)
                        _ = logfile.write('{}\tfrom\t{}\n'.format(dst, src))
                        os.link(src, dst)

                        ok_path = os.path.join(output_dir, output_valid)
                        with open(ok_path, 'w') as ok_file:
                            pass
        # end of rule
