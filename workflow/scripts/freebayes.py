from snakemake.shell import shell

if snakemake.threads == 1:
    freebayes = "freebayes"
else:
    regions = (f"<(fasta_generate_regions.py {snakemake.input.ref} {snakemake.params.chunksize})")
    freebayes = (f"freebayes-parallel {regions} {snakemake.threads}")

shell(
    '{freebayes} -f {snakemake.input.ref} '
        '-p {snakemake.params.ploidy} '
        '-C {snakemake.params.min_count} '
        '-F {snakemake.params.min_freq} '
        '-m {snakemake.params.min_map_qual} '
        '-q {snakemake.params.min_base_qual} '
        '--min-coverage {snakemake.params.min_cov} '
        '{snakemake.input.bam} > {snakemake.output}'
)