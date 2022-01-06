if not config['database']:
    rule BuildDatabase:
        input:
            ref_file = ref_dir + ref,
            ann_file = ref_dir + annotation
        output:
            directory('../resources/snpeff/' + config['reference'])
        params:
            reference = config['reference'],
            ref_name = os.path.splitext(ref)[0],
            snpeff_config = '../config/snpEff.config'
        conda:
            "../envs/snpeff.yaml"
        shell:
            # create directory for which reference database files are to be stored
            'mkdir -p {output} && '
            # copy reference files to created directory
            'cp {input.ref_file} {output}/sequences.fa && '
            'cp {input.ann_file} {output}/genes.gtf && '
            # add the new database info to the SnpEff config file
            'echo -e "\n# {params.reference}\n{params.ref_name}.genome : {params.reference}" >> {params.snpeff_config} && '
            # build the SnpEff database
            'snpEff build -c {params.snpeff_config} -gtf22 -v {params.reference}'

rule SnpEff:
    input:
        calls = '../results/{sample}/final_calls/final_calls_{sample}_snps_pass.vcf',
        db_dir = '../resources/snpeff/' + config['reference']
    output:
        out_calls = '../results/{sample}/final_calls/final_calls_{sample}_snps_snpeff.vcf',
        stats = '../results/{sample}/final_calls/final_calls_{sample}_snps_snpeff.csv'
    params:
        snpeff_config = '../config/snpEff.config',
        ref = config['reference'] if not config['database'] else config['database']
    conda:
        "../envs/snpeff.yaml"
    threads:
        max(config['cores'],config['snpeff_cores'])
    shell:
        'snpEff -t {threads} '
        '-c {params.snpeff_config} '
        '-csvStats {output.stats} '
        '{params.ref} '
        '{input.calls} > {output.out_calls}'