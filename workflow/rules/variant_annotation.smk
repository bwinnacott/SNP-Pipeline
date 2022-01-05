rule GetDatabase:
	input:
		ref_file = ref_dir + ref,
		ann_file = ref_dir + annotation
	output:
		directory('../resources/snpeff/')
	params:
		reference = config['reference'],
		ref_name = os.path.splitext(ref)[0],
		db = config['database'],
		snpeff_config = '../config/snpEff.config'
	log:
		"logs/GetDatabase/test.txt"
	conda:
		"../envs/snpeff.yaml"
	script:
		"../scripts/get_snpeff_database.py"