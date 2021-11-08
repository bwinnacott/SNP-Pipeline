if config['use_mutect2']:
    rule Mutect2:
        input:
            bam = lambda wc: get_input_bam(wc,calling=True),
            ref = '../resources/' + ref
        output:
            '../results/{sample}/{aligner}/Mutect2/raw_variants_recal_{sample}.vcf',
            '../results/{sample}/{aligner}/Mutect2/CollectDirectionalCounts_{sample}.tar.gz'
        conda:
            "../envs/gatk.yaml"
        threads:
            config['cores']
        shell:
            'gatk Mutect2 -R {input.ref} '
            '-I {input.bam} '
            '-O {output[0]} '
            '--f1r2-tar-gz {output[1]} '
            '--native-pair-hmm-threads {threads}'

    rule LearnReadOrientationModel:
        input:
            counts = '../results/{sample}/{aligner}/Mutect2/CollectDirectionalCounts_{sample}.tar.gz'
        output:
            '../results/{sample}/{aligner}/Mutect2/ReadOrientationModel_{sample}.tar.gz'
        conda:
            "../envs/gatk.yaml"
        threads:
            config['cores']
        shell:
            'gatk LearnReadOrientationModel -I {input.counts} '
            '-O {output}'

    rule FilterMutectCalls:
        input:
            read_model = '../results/{sample}/{aligner}/Mutect2/ReadOrientationModel_{sample}.tar.gz',
            var_recal = '../results/{sample}/{aligner}/Mutect2/raw_variants_recal_{sample}.vcf',
            ref = '../resources/' + ref
        output:
            '../results/{sample}/{aligner}/Mutect2/final_Mutect2_variants_{sample}.vcf'
        conda:
            "../envs/gatk.yaml"
        threads:
            config['cores']
        shell:
            'gatk FilterMutectCalls -R {input.ref} '
            '-V {input.var_recal} '
            '-ob-priors {input.read_model} '
            '-O {output}'

if config['use_freebayes']:
    rule Freebayes:
        input:
            ref = '../resources/' + ref,
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
            '../results/{sample}/{aligner}/Freebayes/final_Freebayes_variants_{sample}.vcf'
        conda:
            "../envs/freebayes.yaml"
        threads: 1
        script:
            '../scripts/freebayes.py'