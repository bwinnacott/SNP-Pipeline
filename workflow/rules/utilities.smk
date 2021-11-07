def get_sample_basenames(sample_dir):
    samples = []
    paired = 0
    unpaired = 0
    for f in os.listdir(sample_dir):
        if re.search('_R[1-2]',f):
            paired += 1
            samples.append(re.split('_R[1-2]',f,maxsplit=1)[0])
        else:
            unpaired += 1
            samples.append(f.split('.',1)[0])
    
    assert (paired % 2 == 0), f"Detected paired end file with no match. Ensure all paired files are matched. If \
    not matched, remove '_R1' or '_R2' from file name."
    #print(f"Found total of {paired} paired read files and {unpaired} unpaired read files.")
    
    return samples

def get_sample_exts(sample_dir):
    exts = []
    for f in os.listdir(sample_dir):
        exts.append(f.split('.',1)[1])

    return set(exts)

def get_aligner_input(wildcards,aligner):
    if samples.count(wildcards.sample) == 2:
        sample = set([f for f in os.listdir(sample_dir) if re.split('_R[1-2]',f,maxsplit=1)[0] == wildcards.sample])
        for ext in exts:
            if os.path.exists(sample_dir + '/' + wildcards.sample + '_R1.' + ext):
                R1 = os.path.abspath(sample_dir + '/' + wildcards.sample + '_R1.' + ext)
                R2 = os.path.abspath(sample_dir + '/' + wildcards.sample + '_R2.' + ext)
                if aligner == 'star' or aligner == 'bwa':
                    return R1,R2
                else:
                    return f"-1 {R1} -2 {R2}"
    else:
        for ext in exts:
            if os.path.exists(sample_dir + '/' + wildcards.sample + '.' + ext):
                if aligner == 'star' or aligner == 'bwa':
                    return os.path.abspath(sample_dir + '/' + wildcards.sample + '.' + ext)
                else:
                    return f"-U {os.path.abspath(sample_dir + '/' + wildcards.sample + '.' + ext)}"

def get_star_readfile_command(wildcards):
    for ext in exts:
        if os.path.exists(sample_dir + '/' + wildcards.sample + '_R1.' + ext):
            if ext.endswith('.gz'):
                return '--readFilesCommand zcat'
            else:
                return ''
        elif os.path.exists(sample_dir + '/' + wildcards.sample + '.' + ext):
            if ext.endswith('.gz'):
                return '--readFilesCommand zcat'
            else:
                return ''

def get_aligner_directory(wildcards):
    if mode == 'RNA':
        return '../results/' + wildcards.sample + '/' + wildcards.aligner + '/' + wildcards.aligner + '_output' + '/'
    else:
        return '../results/' + wildcards.sample + '/' + wildcards.aligner + '/'

def get_bam_by_mode(wildcards,ind=False):
    ext = '.bai' if ind else '.bam'
    if mode == 'RNA':
        return '../results/' + wildcards.sample + '/' + wildcards.aligner + '/bam_preprocessing/dedup_rna_reads_' + wildcards.sample + ext
    else:
        return '../results/' + wildcards.sample + '/' + wildcards.aligner + '/bam_preprocessing/dedup_reads_' + wildcards.sample + ext

def get_intersection_output(wildcards):
    if config['output_combined']:
        return '-w1 -o ../results/' + wildcards.sample + '/' + wildcards.aligner + '/final_calls/final_calls.vcf'
    else:
        return '-p ../results/' + wildcards.sample + '/' + wildcards.aligner + '/final_calls'