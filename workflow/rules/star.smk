rule star_index:
    input:
        ref = '../resources/' + ref,
        gtf = '../resources/' + gtf
    params:
        sjdbOverhang = config['sjdbOverhang'],
        sjdbGTFtag = config['sjdbGTFtag'],
        genomeSAindexNbases = config['genomeSAindexNbases']
    output:
        directory(os.path.splitext(ref)[0] + '_star_index')
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        'mkdir ../resources/{output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir ../resources/{output} '
        '--genomeFastaFiles {input.ref} '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang {params.sjdbOverhang} '
        '--sjdbGTFtagExonParentTranscript {params.sjdbGTFtag} '
        '--genomeSAindexNbases {params.genomeSAindexNbases}'

rule star_align:
    input:
        directory(os.path.splitext(ref)[0] + '_star_index'),
        refdir = '../resources/' + os.path.splitext(ref)[0] + '_star_index'
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
    threads: 8
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
        '--outSAMtype BAM SortedByCoordinate '
        '--quantMode {params.quantMode} && '
        'mv Aligned.sortedByCoord.out.bam {wildcards.sample}_sorted.bam && cd ../../../../workflow'

