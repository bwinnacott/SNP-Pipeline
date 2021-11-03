rule HaplotypeCaller:
    input:
        bam = get_bam_by_mode,
        ref = '../resources/' + ref
    output:
        '../results/{sample}/{aligner}/bqsr/raw_variants_{sample}.vcf'
    conda:
        "../envs/gatk.yaml"
    threads:
        config['cores']
    shell:
        'gatk HaplotypeCaller -R {input.ref} '
        '-I {input.bam} '
        '-O {output} '
        '--native-pair-hmm-threads {threads}'

rule SelectVariants_pass1:
    input:
        vcf = '../results/{sample}/{aligner}/bqsr/raw_variants_{sample}.vcf',
        ref = '../resources/' + ref
    output:
        '../results/{sample}/{aligner}/bqsr/raw_snps_{sample}.vcf',
        '../results/{sample}/{aligner}/bqsr/raw_indels_{sample}.vcf'
    conda:
        "../envs/gatk.yaml"
    threads:
        config['cores']
    shell:
        'gatk SelectVariants -R {input.ref} '
        '-V {input.vcf} '
        '--select-type-to-exclude SNP '
        '-O {output[0]} && '
        'gatk SelectVariants -R {input.ref} '
        '-V {input.vcf} '
        '--select-type-to-exclude INDEL '
        '-O {output[1]}'

rule VariantFiltration:
    input:
        snps_vcf = '../results/{sample}/{aligner}/bqsr/raw_snps_{sample}.vcf',
        indels_vcf = '../results/{sample}/{aligner}/bqsr/raw_indels_{sample}.vcf',
        ref = '../resources/' + ref
    params:
        QD_filter_snp = config['QD_filter_snp'],
        FS_filter_snp = config['FS_filter_snp'],
        MQ_filter_snp = config['MQ_filter_snp'],
        SOR_filter_snp = config['SOR_filter_snp'],
        MQRankSum_filter_snp = config['MQRankSum_filter_snp'],
        ReadPosRankSum_filter_snp = config['ReadPosRankSum_filter_snp'],
        QD_filter_indel = config['QD_filter_indel'],
        FS_filter_indel = config['FS_filter_indel'],
        SOR_filter_indel = config['SOR_filter_indel']
    output:
        '../results/{sample}/{aligner}/bqsr/filtered_snps_{sample}.vcf',
        '../results/{sample}/{aligner}/bqsr/filtered_indels_{sample}.vcf'
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk VariantFiltration -R {input.ref} '
        '-V {input.snps_vcf} '
        '-O {output[0]} '
        '-filter-name "QD_filter" -filter "{params.QD_filter_snp}" '
        '-filter-name "FS_filter" -filter "{params.FS_filter_snp}" '
        '-filter-name "MQ_filter" -filter "{params.MQ_filter_snp}" '
        '-filter-name "SOR_filter" -filter "{params.SOR_filter_snp}" '
        '-filter-name "MQRankSum_filter" -filter "{params.MQRankSum_filter_snp}" '
        '-filter-name "ReadPosRankSum_filter" -filter "{params.ReadPosRankSum_filter_snp}" && '
        'gatk VariantFiltration -R {input.ref} '
        '-V {input.indels_vcf} '
        '-O {output[1]} '
        '-filter-name "QD_filter" -filter "{params.QD_filter_indel}" '
        '-filter-name "FS_filter" -filter "{params.FS_filter_indel}" '
        '-filter-name "SOR_filter" -filter "{params.SOR_filter_indel}"'

rule SelectVariants_pass2:
    input:
        filtered_snps_vcf = '../results/{sample}/{aligner}/bqsr/filtered_snps_{sample}.vcf',
        filtered_indels_vcf = '../results/{sample}/{aligner}/bqsr/filtered_indels_{sample}.vcf'
    output:
        '../results/{sample}/{aligner}/bqsr/bqsr_snps_{sample}.vcf',
        '../results/{sample}/{aligner}/bqsr/bqsr_indels_{sample}.vcf'
    conda:
        "../envs/gatk.yaml"
    threads:
        config['cores']
    shell:
        'gatk SelectVariants --exclude-filtered true '
        '-V {input.filtered_snps_vcf} '
        '-O {output[0]} && '
        'gatk SelectVariants --exclude-filtered true '
        '-V {input.filtered_indels_vcf} '
        '-O {output[1]}'

rule BaseRecalibrator:
    input:
        bam = get_bam_by_mode,
        bqsr_snps = '../results/{sample}/{aligner}/bqsr/bqsr_snps_{sample}.vcf',
        bqsr_indels = '../results/{sample}/{aligner}/bqsr/bqsr_indels_{sample}.vcf',
        ref = '../resources/' + ref
    output:
        '../results/{sample}/{aligner}/bqsr/recal_data_{sample}.table'
    conda:
        "../envs/gatk.yaml"
    threads:
        config['cores']
    shell:
        'gatk BaseRecalibrator -R {input.ref} '
        '-I {input.bam} '
        '--known-sites {input.bqsr_snps} '
        '--known-sites {input.bqsr_indels} '
        '-O {output}'

rule ApplyBQSR:
    input:
        bam = get_bam_by_mode,
        bqsr = '../results/{sample}/{aligner}/bqsr/recal_data_{sample}.table',
        ref = '../resources/' + ref
    output:
        '../results/{sample}/{aligner}/bqsr/recal_reads_{sample}.bam'
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk ApplyBQSR -R {input.ref} '
        '-I {input.bam} '
        '-bqsr {input.bqsr} '
        '-O {output}'