
rule samtools_bai_index_bam:
    input:
        '{filepath}.bam'
    output:
        '{filepath}.bam.bai'
    threads: 4
    shell:
        'samtools index -b -@ {threads} {input}'


rule samtools_csi_index_bam:
    input:
        '{filepath}.bam'
    output:
        '{filepath}.bam.csi'
    threads: 4
    shell:
        'samtools index -c -@ {threads} {input}'


rule samtools_position_sort_bam:
    input:
        '{filepath}.bam'
    output:
        '{filepath}.psort.bam'
    log:
        'log/{filepath}.psort.log'
    benchmark:
        'run/{filepath}.psort.rsrc'
    threads: 6
    resources:
        mem_mb= 12 * 1024  # 2G per thread
    shell:
        'samtools sort -l 9 -m 2G -@ {threads} --output-fmt BAM -o {output} {input} &> {log}'