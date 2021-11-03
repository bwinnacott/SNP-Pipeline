rule hisat2_index:
    input:
        ref = '../resources/' + ref
    params:
        prefix = os.path.splitext(ref)[0]
    output: 
        directory('../resources/' + os.path.splitext(ref)[0] + '_hisat2_index')
    conda:
        "../envs/hisat2.yaml"
    threads: 8
    shell:
        'mkdir {output} && '
        'hisat2-build -p {threads} '
        '{input.ref} {output}/{params.prefix}'

rule hisat2_align:
    input:
        gtf = '../resources/' + gtf,
        refdir = '../resources/' + os.path.splitext(ref)[0] + '_hisat2_index'
    params:
        sample = lambda wc: get_aligner_input(wc,aligner='hisat2'),
        outdir = '../results/{sample}/hisat2/hisat2_output',
        prefix = os.path.splitext(ref)[0],
        minintronlen = config['minintronlen'],
        maxintronlen = config['maxintronlen']
    output:
        '../results/{sample}/hisat2/hisat2_output/{sample}_sorted.bam'
    conda:
        "../envs/hisat2.yaml"
    threads: 8
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'ss_script=$(which hisat2_extract_splice_sites.py) && '
        '$ss_script ../../../{input.gtf} > splicesites.txt && '
        'hisat2 -p {threads} '
        '-x ../../../{input.refdir}/{params.prefix} '
        '{params.sample} '
        '--min-intronlen {params.minintronlen} '
        '--max-intronlen {params.maxintronlen} '
        '--known-splicesite-infile splicesites.txt '
        '-S {wildcards.sample}.sam && '
        'samtools view -bS {wildcards.sample}.sam > {wildcards.sample}.bam && '
        'samtools sort {wildcards.sample}.bam -o {wildcards.sample}_sorted.bam && '
        'rm {wildcards.sample}.bam && cd ../../../../workflow'
