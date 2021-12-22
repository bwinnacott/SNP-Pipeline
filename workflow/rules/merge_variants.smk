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
        nfiles = config['nfiles'],
        out_dir = '../results/{sample}/final_calls/'
        #output = get_intersection_output
    output:
        '../results/{sample}/final_calls/final_calls_{sample}.vcf'
    conda:
        "../envs/merge.yaml"
    shell:
        'mkdir -p {params.out_dir} && '
        'bcftools isec -n +{params.nfiles} -w 1 -o {output} {input}'

rule FilterIntersection:
    input:
        '../results/{sample}/final_calls/final_calls_{sample}.vcf'
    output:
        '../results/{sample}/final_calls/final_calls_{sample}_snps.vcf',
        '../results/{sample}/final_calls/final_calls_{sample}_indels.vcf',
        '../results/{sample}/final_calls/final_calls_{sample}_snps_pass.vcf',
        '../results/{sample}/final_calls/final_calls_{sample}_indels_pass.vcf'
    conda:
        "../envs/bcftools.yaml"
    shell:
        'bcftools view --type snps {input} > {output[0]} && '
        'bcftools view --type indels {input} > {output[1]} && '
        'bcftools view -f .,PASS {output[0]} > {output[2]} && '
        'bcftools view -f .,PASS {output[1]} > {output[3]}'