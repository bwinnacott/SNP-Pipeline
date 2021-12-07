if config['use_mutect2']:
    rule Mutect2:
        input:
            bam = lambda wc: get_input_bam(wc,calling=True),
            ref = ref_dir + ref
        output:
            '../results/{sample}/{aligner}/Mutect2/raw_variants_Mutect2_{sample}.vcf',
            '../results/{sample}/{aligner}/Mutect2/CollectDirectionalCounts_{sample}.tar.gz'
        conda:
            "../envs/gatk.yaml"
        threads:
            max(config['cores'],config['mutect_cores'])
        shell:
            'gatk Mutect2 -R {input.ref} '
            '-I {input.bam} '
            '-O {output[0]} '
            '--f1r2-tar-gz {output[1]} '
            '--native-pair-hmm-threads {threads}'

    #rule GetPileupSummaries:
        #input:
            #bam = lambda wc: get_input_bam(wc,calling=True),
            #raw_variants = '../results/{sample}/{aligner}/Mutect2/raw_variants_Mutect2_{sample}.vcf'
        #output:
            #'../results/{sample}/{aligner}/Mutect2/PileupTable_{sample}.tsv'
        #conda:
            #"../envs/gatk.yaml"
        #threads:
            #config['cores']
        #shell:
            #'gatk GetPileupSummaries -I {input.bam} '
            #'-V {input.raw_variants} '
            #'-L {input.raw_variants} '
            #'-O {output}'

    rule LearnReadOrientationModel:
        input:
            counts = '../results/{sample}/{aligner}/Mutect2/CollectDirectionalCounts_{sample}.tar.gz'
        output:
            '../results/{sample}/{aligner}/Mutect2/ReadOrientationModel_{sample}.tar.gz'
        conda:
            "../envs/gatk.yaml"
        shell:
            'gatk LearnReadOrientationModel -I {input.counts} '
            '-O {output}'

    rule FilterMutectCalls:
        input:
            read_model = '../results/{sample}/{aligner}/Mutect2/ReadOrientationModel_{sample}.tar.gz',
            raw_variants = '../results/{sample}/{aligner}/Mutect2/raw_variants_Mutect2_{sample}.vcf',
            ref = ref_dir + ref
        output:
            '../results/{sample}/{aligner}/Mutect2/final_Mutect2_variants_{sample}.vcf'
        conda:
            "../envs/gatk.yaml"
        shell:
            'gatk FilterMutectCalls -R {input.ref} '
            '-V {input.raw_variants} '
            '-ob-priors {input.read_model} '
            '-O {output}'

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
            qual = config['QUAL']
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
            max_dp = config['max_depth'],
            min_MQ = config['min_MQ'],
            min_BQ = config['min_BQ'],
            caller = config['caller'],
            qual_filter = config['qual_filter'],
            depth_filter = config['depth_filter']
        output:
            '../results/{sample}/{aligner}/Bcftools/final_Bcftools_variants_{sample}.vcf'
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
            'bcftools call -Ou -v {params.caller} | '
            'bcftools filter -s LowQual -i "{params.qual_filter} & {params.depth_filter}" > {output}'