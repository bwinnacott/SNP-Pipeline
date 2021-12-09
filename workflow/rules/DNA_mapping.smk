rule BWAIndex:
    input:
        ref_dir + ref
    params: 
        algorithm = config['algorithm']
    output:
        ref_dir + ref + '.amb',
        ref_dir + ref + '.ann',
        ref_dir + ref + '.bwt',
        ref_dir + ref + '.pac',
        ref_dir + ref + '.sa'
    conda:
        "../envs/bwa.yaml"
    shell:
        'bwa index -a {params.algorithm} {input}'

rule BWAMem:
    input:
        sample_dir + '/running_pipeline.txt',
        ref_dir + ref + '.amb',
        ref_dir + ref + '.ann',
        ref_dir + ref + '.bwt',
        ref_dir + ref + '.pac',
        ref_dir + ref + '.sa',
        ref = ref_dir + ref
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
    threads:
        max(config['cores'],config['bwa_mem_cores'])
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'bwa mem -t {threads} '
        '-R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{params.PL}\\tLB:{params.LB}\\tPU:{params.PU}" '
        '../../{input.ref} '
        '{params.sample} > {wildcards.sample}.sam && '
        'samtools view -bS {wildcards.sample}.sam > {wildcards.sample}.bam && '
        'rm {wildcards.sample}.sam && '
        'samtools sort {wildcards.sample}.bam -o {wildcards.sample}_sorted.bam -T sort && '
        'rm {wildcards.sample}.bam && cd ../../../workflow'