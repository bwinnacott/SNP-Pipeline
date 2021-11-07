rule compress_and_index:
    input:
        '../results/{sample}/{aligner}/Mutect2/FilteredMutectVariants_{sample}.vcf',
        '../results/{sample}/{aligner}/Freebayes/freebayes_variants_{sample}.vcf'
    output:
        '../results/{sample}/{aligner}/Mutect2/FilteredMutectVariants_{sample}.vcf.gz',
        '../results/{sample}/{aligner}/Freebayes/freebayes_variants_{sample}.vcf.gz',
        '../results/{sample}/{aligner}/Mutect2/FilteredMutectVariants_{sample}.vcf.gz.tbi',
        '../results/{sample}/{aligner}/Freebayes/freebayes_variants_{sample}.vcf.gz.tbi'
    conda:
        "../envs/merge.yaml"
    threads: 1
    shell:
        'bgzip {input}[0] && '
        'bgzip {input}[1] && '
        'tabix {output}[0] && '
        'tabix {output}[1]'

rule intersection:
    input:
        '../results/{sample}/{aligner}/Mutect2/FilteredMutectVariants_{sample}.vcf.gz',
        '../results/{sample}/{aligner}/Freebayes/freebayes_variants_{sample}.vcf.gz'
    params:
        nfiles = config['nfiles'],
        output = get_intersection_output
    output:
        directory('../results/{sample}/{aligner}/final_calls/')
    conda:
        "../envs/merge.yaml"
    threads: 1
    shell:
        'bcftools isec -n +{params.nfiles} {params.output} {input}'