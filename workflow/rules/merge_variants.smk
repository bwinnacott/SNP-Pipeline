rule compress_and_index:
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

rule intersection:
    input:
        expand('../results/{sample}/{aligner}/{caller}/final_{caller}_variants_{sample}.vcf.gz',aligner=aligner,caller=callers,allow_missing=True)
    params:
        nfiles = config['nfiles'],
        output = get_intersection_output
    output:
        directory('../results/{sample}/final_calls/')
    conda:
        "../envs/merge.yaml"
    shell:
        # might think about writing script to handle different scenarios for intersection output
        # i.e., renaming the output files if all files are written out...
        'mkdir -p {output} && '
        'bcftools isec -n +{params.nfiles} {params.output} {input}'
