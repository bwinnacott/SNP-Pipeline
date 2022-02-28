if config['use_haplotypecaller']:
    rule HaplotypeCaller:
        input:
            ref = ref_dir + ref,
            bam = lambda wc: get_input_bam(wc,calling=True)
        params:
            RNA_params = get_hapcaller_rna_params,
            ploidy = config['ploidy']
        output:
            '../results/{sample}/{aligner}/HaplotypeCaller/raw_variants_HaplotypeCaller_{sample}.vcf'
        conda:
            "../envs/gatk.yaml"
        threads: 
            max(config['cores'],config['hap_caller_cores'])
        shell:
            'gatk HaplotypeCaller -R {input.ref} '
            '-I {input.bam} '
            '-ploidy {params.ploidy} '
            '-O {output} '
            '--native-pair-hmm-threads {threads} '
            '{params.RNA_params}'

    rule FilterHapVariants:
        input:
            vcf = '../results/{sample}/{aligner}/HaplotypeCaller/raw_variants_HaplotypeCaller_{sample}.vcf',
            ref = ref_dir + ref
        params:
            qual = config['qual_filter_hapcaller'],
            QD_filter_snp = config['QD_filter_snp'],
            FS_filter_snp = config['FS_filter_snp'],
            MQ_filter_snp = config['MQ_filter_snp'],
            SOR_filter_snp = config['SOR_filter_snp'],
            MQRankSum_filter_snp = config['MQRankSum_filter_snp'],
            ReadPosRankSum_filter_snp = config['ReadPosRankSum_filter_snp']
        output:
            '../results/{sample}/{aligner}/HaplotypeCaller/filtered_HaplotypeCaller_variants_{sample}.vcf',
            '../results/{sample}/{aligner}/HaplotypeCaller/final_HaplotypeCaller_variants_{sample}.vcf'
        conda:
            "../envs/gatk.yaml"
        shell:
            'gatk VariantFiltration -R {input.ref} '
            '-V {input.vcf} '
            '-O {output[0]} '
            '-filter "{params.qual}" -filter-name "Quality" '
            '-filter "{params.QD_filter_snp}" -filter-name "QD_filter" '
            '-filter "{params.FS_filter_snp}" -filter-name "FS_filter" '
            '-filter "{params.MQ_filter_snp}" -filter-name "MQ_filter" '
            '-filter "{params.SOR_filter_snp}" -filter-name "SOR_filter" '
            '-filter "{params.MQRankSum_filter_snp}" -filter-name "MQRankSum_filter" '
            '-filter "{params.ReadPosRankSum_filter_snp}" -filter-name "ReadPosRankSum_filter" && '
            'gatk SelectVariants --exclude-filtered true '
            '-V {output[0]} '
            '-O {output[1]}'

if config['use_freebayes']:
    rule Freebayes:
        input:
            ref = ref_dir + ref,
            bam = lambda wc: get_input_bam(wc,calling=True),
            index = lambda wc: get_input_bam(wc,calling=True,ind=True)
        params:
            ploidy = config['ploidy'],
            min_count = config['min_alternate_count'],
            min_freq = config['min_alternate_freq'],
            min_cov = config['min_coverage'],
            min_map_qual = config['min_map_qual'],
            min_base_qual = config['min_base_qual'],
            chunksize = config['chunksize']
        output:
            '../results/{sample}/{aligner}/Freebayes/raw_variants_Freebayes_{sample}.vcf'
        conda:
            "../envs/freebayes.yaml"
        threads:
            max(config['cores'],config['free_cores'])
        script:
            '../scripts/freebayes.py'

    rule FilterFreebayesCalls:
        input:
            '../results/{sample}/{aligner}/Freebayes/raw_variants_Freebayes_{sample}.vcf'
        params:
            qual = config['qual_filter_free']
        output:
            '../results/{sample}/{aligner}/Freebayes/final_Freebayes_variants_{sample}.vcf'
        conda:
            "../envs/freebayes.yaml"
        shell:
            'vcffilter -f "QUAL > {params.qual}" {input} > {output}'

if config['use_bcftools']:
    rule BcftoolsCall:
        input:
            bam = lambda wc: get_input_bam(wc,calling=True),
            ref = ref_dir + ref
        params:
            ploidy = '--ploidy 1 ' if config['ploidy'] == 1 else '',
            max_dp = config['max_depth'],
            min_MQ = config['min_MQ'],
            min_BQ = config['min_BQ'],
            caller = config['caller']
        output:
            '../results/{sample}/{aligner}/Bcftools/raw_variants_Bcftools_{sample}.vcf'
        conda:
            "../envs/bcftools.yaml"
        threads:
            max(config['cores'],config['bcftools_cores'])
        shell:
            'bcftools mpileup -Ou '
            '--max-depth {params.max_dp} '
            '--min-MQ {params.min_MQ} '
            '--min-BQ {params.min_BQ} '
            '--threads {threads} '
            '-f {input.ref} {input.bam} | '
            'bcftools call -Ov -v {params.caller} {params.ploidy}> {output}'

    rule FilterBcftoolsCalls:
        input:
            '../results/{sample}/{aligner}/Bcftools/raw_variants_Bcftools_{sample}.vcf'
        params:
            qual_filter = config['qual_filter_bcf'],
            depth_filter = config['depth_filter_bcf']
        output:
            '../results/{sample}/{aligner}/Bcftools/filtered_Bcftools_variants_{sample}.vcf',
            '../results/{sample}/{aligner}/Bcftools/final_Bcftools_variants_{sample}.vcf'
        conda:
            "../envs/bcftools.yaml"
        shell:
            'bcftools filter -s LowQual -i "{params.qual_filter} & {params.depth_filter}" {input} > {output[0]} && '
            'bcftools view -f .,PASS {output[0]} > {output[1]}'