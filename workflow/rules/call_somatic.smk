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
            '../results/{sample}/{aligner}/Mutect2/filtered_Mutect2_variants_{sample}.vcf',
            '../results/{sample}/{aligner}/Mutect2/final_Mutect2_variants_{sample}.vcf'
        conda:
            "../envs/gatk.yaml"
        shell:
            'gatk FilterMutectCalls -R {input.ref} '
            '-V {input.raw_variants} '
            '-ob-priors {input.read_model} '
            '-O {output[0]} && '
            'bcftools view -Ov -f PASS {output[0]} > {output[1]}'

