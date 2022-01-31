rule CompressAndIndex:
    input:
        '../results/{sample}/{aligner}/{caller}/final_{caller}_variants_{sample}.vcf'
    output:
        '../results/{sample}/{aligner}/{caller}/final_{caller}_variants_{sample}.vcf.gz',
        '../results/{sample}/{aligner}/{caller}/final_{caller}_variants_{sample}.vcf.gz.tbi'
    conda:
        "../envs/merge.yaml"
    shell:
        'bgzip {input} && '
        'tabix {output[0]}'

rule Intersection:
    input:
        expand('../results/{sample}/{aligner}/{caller}/final_{caller}_variants_{sample}.vcf.gz',aligner=aligner,caller=callers,allow_missing=True)
    params:
        nfiles = lambda wc: get_intersection_num(wc),
        out_dir = '../results/{sample}/final_calls/'
        #output = get_intersection_output
    output:
        '../results/{sample}/final_calls/final_calls_{sample}.vcf'
    conda:
        "../envs/merge.yaml"
    shell:
        'bcftools isec -Oz -n +{params.nfiles} -p {params.out_dir} {input} && '
        'bcftools concat -a -D -o {output} {params.out_dir}*.vcf.gz && '
        'rm {params.out_dir}*.vcf.gz*'

rule SeparateVariants:
    input:
        '../results/{sample}/final_calls/final_calls_{sample}.vcf'
    output:
        '../results/{sample}/final_calls/final_calls_{sample}_snps.vcf',
        '../results/{sample}/final_calls/final_calls_{sample}_indels.vcf'
    conda:
        "../envs/bcftools.yaml"
    shell:
        'bcftools view --type snps {input} > {output[0]} && '
        'bcftools view --type indels {input} > {output[1]}'