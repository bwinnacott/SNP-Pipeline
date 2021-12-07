if annotation.endswith(tuple(['.gff3','.gff'])):
    annotation_split = os.path.splitext(annotation)
    annotation = annotation_split[0] + '.gtf'
    rule ConvertGffToGtf:
        input:
            gff = ref_dir + ''.join(annotation_split)
        output:
            ref_dir + annotation
        conda:
            "../envs/star.yaml"
        shell:
            'gffread {input.gff} -T -o {output}'

rule StarIndex:
    input:
        ref = ref_dir + ref,
        annotation = ref_dir + annotation
    params:
        sjdbOverhang = config['sjdbOverhang'],
        sjdbGTFtag = config['sjdbGTFtag'],
        genomeSAindexNbases = config['genomeSAindexNbases']
    output:
        directory(ref_dir + os.path.splitext(ref)[0] + '_star_index')
    conda:
        "../envs/star.yaml"
    threads:
        max(config['cores'],config['star_index_cores'])
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.ref} '
        '--sjdbGTFfile {input.annotation} '
        '--sjdbOverhang {params.sjdbOverhang} '
        '--sjdbGTFtagExonParentTranscript {params.sjdbGTFtag} '
        '--genomeSAindexNbases {params.genomeSAindexNbases}'

rule Hisat2Index:
    input:
        ref = ref_dir + ref
    params:
        prefix = os.path.splitext(ref)[0]
    output: 
        directory(ref_dir + os.path.splitext(ref)[0] + '_hisat2_index')
    conda:
        "../envs/hisat2.yaml"
    threads:
        max(config['cores'],config['hisat_index_cores'])
    shell:
        'mkdir {output} && '
        'hisat2-build -p {threads} '
        '{input.ref} {output}/{params.prefix}'