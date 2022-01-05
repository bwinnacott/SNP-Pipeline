from snakemake.shell import shell

if not snakemake.params.db:
	shell(
		# create directory for which reference database files are to be stored
		'mkdir -p {snakemake.output}/{snakemake.params.ref_name} && '
		# copy reference files to created directory
		'cp {snakemake.input.ref_file} {snakemake.output}/{snakemake.params.ref_name}/sequences.fa && '
		'cp {snakemake.input.ann_file} {snakemake.output}/{snakemake.params.ref_name}/genes.gtf && '
		# add the new database info to the SnpEff config file
		'echo -e "\n# {snakemake.params.ref_name}\n{snakemake.params.ref_name}.genome : {snakemake.params.reference}" >> {snakemake.params.snpeff_config} && '
		# build the SnpEff database
		'snpEff build -c {snakemake.params.snpeff_config} -gtf22 -v {snakemake.params.ref_name}'
	)
else:
	shell(
		'snpEff download -dataDir {snakemake.output}/{snakemake.params.ref_name} {snakemake.params.db}'
	)