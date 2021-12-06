rule star_align:
    input:
        refdir = ref_dir + os.path.splitext(ref)[0] + '_star_index'
    params:
        sample = lambda wc: get_aligner_input(wc,aligner='star'),
        outdir = '../results/{sample}/star/star_output',
        Command = get_star_readfile_command,
        Type = config['outFilterType'],
        MultimapNmax = config['outFilterMultimapNmax'],
        overhangMin = config['alignSJoverhangMin'],
        DBoverhangMin = config['alignSJDBoverhangMin'],
        MismatchNmax = config['outFilterMismatchNmax'],
        MismatchNoverReadLmax = config['outFilterMismatchNoverReadLmax'],
        IntronMin = config['alignIntronMin'],
        IntronMax = config['alignIntronMax'],
        MatesGapMax = config['alignMatesGapMax'],
        quantMode = config['quantMode']
    output:
        '../results/{sample}/star/star_output/{sample}_sorted.bam'
    conda:
        "../envs/star.yaml"
    threads:
        max(config['cores'],config['star_align_cores'])
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'STAR --runThreadN {threads} '
        '{params.Command} '
        '--genomeDir ../../../{input.refdir} '
        '--readFilesIn {params.sample} '
        '--outFilterType {params.Type} '
        '--outFilterMultimapNmax {params.MultimapNmax} '
        '--alignSJoverhangMin {params.overhangMin} '
        '--alignSJDBoverhangMin {params.DBoverhangMin} '
        '--outFilterMismatchNmax {params.MismatchNmax} '
        '--outFilterMismatchNoverReadLmax {params.MismatchNoverReadLmax} '
        '--alignIntronMin {params.IntronMin} '
        '--alignIntronMax {params.IntronMax} '
        '--alignMatesGapMax {params.MatesGapMax} '
        '--twopassMode Basic '
        '--outSAMtype BAM SortedByCoordinate '
        '--quantMode {params.quantMode} && '
        'mv Aligned.sortedByCoord.out.bam {wildcards.sample}_sorted.bam && cd ../../../../workflow'

rule hisat2_align:
    input:
        annotation = ref_dir + annotation,
        refdir = ref_dir + os.path.splitext(ref)[0] + '_hisat2_index'
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
    threads:
        max(config['cores'],config['hisat_align_cores'])
    shell:
        'mkdir -p {params.outdir} && '
        'cd {params.outdir} && '
        'ss_script=$(which hisat2_extract_splice_sites.py) && '
        '$ss_script ../../../{input.annotation} > splicesites.txt && '
        'hisat2 -p {threads} '
        '-x ../../../{input.refdir}/{params.prefix} '
        '{params.sample} '
        '--min-intronlen {params.minintronlen} '
        '--max-intronlen {params.maxintronlen} '
        '--known-splicesite-infile splicesites.txt '
        '-S {wildcards.sample}.sam && '
        'samtools view -bS {wildcards.sample}.sam > {wildcards.sample}.bam && '
        'rm {wildcards.sample}.sam && '
        'samtools sort {wildcards.sample}.bam -o {wildcards.sample}_sorted.bam && '
        'rm {wildcards.sample}.bam && cd ../../../../workflow'
