rule index_ref:
    input:
        ref = ref_dir + ref
    output:
        ref_dir + os.path.splitext(ref)[0] + '.dict',
        ref_dir + ref + '.fai'
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk CreateSequenceDictionary -R ../resources/{input.ref} && '
        'samtools faidx ../resources/{input.ref}'

rule mark_duplicates:
    input:
        ref_dir + os.path.splitext(ref)[0] + '.dict',
        ref_dir + ref + '.fai',
        bam = lambda wc: get_aligner_directory(wc) + '{sample}_sorted.bam',
        ref = ref_dir + ref
    params:
        out_dir = '../results/{sample}/{aligner}/bam_preprocessing'
    output:
        '../results/{sample}/{aligner}/bam_preprocessing/dedup_reads_{sample}.bam',
        '../results/{sample}/{aligner}/bam_preprocessing/dedup_reads_{sample}.bai'
    conda:
        "../envs/gatk.yaml"
    threads:
        config['cores']
    shell:
        'gatk MarkDuplicates -I {input.bam} '
        '-O {output[0]} '
        '-M {params.out_dir}/MarkedDupMetrics_{wildcards.sample}.txt '
        '-ASO coordinate '
        '--CREATE_INDEX true'

if mode == 'RNA':
    rule split_cigar_reads:
        input:
            ref = ref_dir + ref,
            bam = '../results/{sample}/{aligner}/bam_preprocessing/dedup_reads_{sample}.bam'
        output:
            '../results/{sample}/{aligner}/bam_preprocessing/dedup_split_reads_{sample}.bam'
        conda:
            "../envs/gatk.yaml"
        threads:
            config['cores']
        shell:
            'gatk SplitNCigarReads --reference {input.ref} '
            '--input {input.bam} '
            '--output {output}'

    rule add_or_replace_rgps:
        input:
            '../results/{sample}/{aligner}/bam_preprocessing/dedup_split_reads_{sample}.bam'
        params:
            id = '{sample}',
            lb = os.path.splitext(ref)[0],
            pl = config['pl'],
            pu = config['pu'],
            sm = '{sample}'
        output:
            '../results/{sample}/{aligner}/bam_preprocessing/dedup_rna_reads_{sample}.bam',
            '../results/{sample}/{aligner}/bam_preprocessing/dedup_rna_reads_{sample}.bai'
        conda:
            "../envs/gatk.yaml"
        threads:
            config['cores']
        shell:
            'gatk AddOrReplaceReadGroups -INPUT {input} '
            '-OUTPUT {output[0]} '
            '-RGID {params.id} -RGLB {params.lb} -RGPL {params.pl} -RGPU {params.pu} -RGSM {params.sm} '
            '-SORT_ORDER coordinate '
            '--CREATE_INDEX true'