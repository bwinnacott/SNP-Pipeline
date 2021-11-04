rule bwa_index:
    input:
        ref = '../resources/' + ref
    params:
        algorithm = config['algorithm']
    output:
        directory('../resources/' + os.path.splitext(ref)[0] + '_bwa_index')
    conda:
        "../envs/bwa.yaml"
    shell:
        'mkdir {output} && '
        'cd {output} && '
        'bwa index -a {params.algorithm} {input.ref}'

rule bwa_mem:
    input:
        refdir = '../resources/' + os.path.splitext(ref)[0] + '_bwa_index'
        ref = '../resources/' + ref
    params:
        sample = lambda wc: get_aligner_input(wc,aligner='bwa'),
        outdir = '../results/{sample}/bwa',
        RG = "@RG\tID:'{sample}'\tSM:'{sample}'\tPL:{config['pl']}\tLB:{os.path.splitext(ref)[0]}\tPU:{config['lb']}"
    output:
        '../results/{sample}/bwa/{sample}_sorted.bam'
    conda:
        "../envs/bwa.yaml"
    threads: 8
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'bwa mem -t {threads} '
        '-R {params.RG} '
        '../{input.ref} '
        '{params.sample} > {wildcards.sample}.sam && '
        'samtools view -bS {wildcards.sample}.sam > {wildcards.sample}.bam && '
        'samtools sort {wildcards.sample}.bam -o {wildcards.sample}_sorted.bam -T sort && '
        'rm {wildcards.sample}.bam && cd ../../workflow'