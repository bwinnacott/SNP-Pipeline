rule GetContigInfo:
    input:
        '../results/{sample}/final_calls/final_calls_{sample}_snps_pass.vcf'
    output:
        '../results/{sample}/final_calls/contig_ids_{sample}.txt',
        '../results/{sample}/final_calls/contig_lens_{sample}.txt',
        '../results/{sample}/final_calls/final_calls_{sample}_snps_pass_noheader.vcf'
    shell:
        'cat {input} | grep "^##contig" | cut -d , -f 1 | cut -d = -f 3 > {output[0]} && '
        'cat {input} | grep "^##contig" | cut -d , -f 2 | cut -d = -f 2 | sed "s/>//" > {output[1]} && '
        'cat {input} | grep -v "^##" | sed "s/^#//" > {output[2]}'

rule CheckPloidy:
    input:
        '../results/{sample}/final_calls/final_calls_{sample}_snps_pass_noheader.vcf',
        '../results/{sample}/final_calls/contig_ids_{sample}.txt',
        '../results/{sample}/final_calls/contig_lens_{sample}.txt'
    output:
        '../results/{sample}/final_calls/ploidy_check_report.html'
    conda:
        "../envs/ploidy_check.yaml"
    script:
        '../scripts/ploidy_check_report.Rmd'