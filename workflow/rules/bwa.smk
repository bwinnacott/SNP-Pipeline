rule bwa_index:
    input:
        ref = '../resources/' + ref
    params: 
        algorithm = config['algorithm']
    output:
        '../resources/' + ref + '.amb',
        '../resources/' + ref + '.ann',
        '../resources/' + ref + '.bwt',
        '../resources/' + ref + '.pac',
        '../resources/' + ref + '.sa'
    conda:
        "../envs/bwa.yaml"
    shell:
        'bwa index -a {params.algorithm} {input.ref}'

rule bwa_mem:
    input:
        '../resources/' + ref + '.amb',
        '../resources/' + ref + '.ann',
        '../resources/' + ref + '.bwt',
        '../resources/' + ref + '.pac',
        '../resources/' + ref + '.sa',
        ref = '../resources/' + ref
    params:
        sample = lambda wc: get_aligner_input(wc,aligner='bwa'),
        outdir = '../results/{sample}/bwa',
        PL = config['pl'],
        LB = os.path.splitext(ref)[0],
        PU = config['pu']
    output:
        '../results/{sample}/bwa/{sample}_sorted.bam'
    conda:
        "../envs/bwa.yaml"
    threads: 8
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'bwa mem -t {threads} '
        '-R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{params.PL}\\tLB:{params.LB}\\tPU:{params.PU}" '
        '../../{input.ref} '
        '{params.sample} > {wildcards.sample}.sam && '
        'samtools view -bS {wildcards.sample}.sam > {wildcards.sample}.bam && '
        'samtools sort {wildcards.sample}.bam -o {wildcards.sample}_sorted.bam -T sort && '
        'rm {wildcards.sample}.bam && cd ../../../workflow'