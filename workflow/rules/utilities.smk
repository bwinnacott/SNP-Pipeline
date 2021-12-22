import os
import re
import gzip

#rule GetSampleInfo:
#    output:
#        temp(sample_dir + '/running_pipeline.txt')
#    run:
#        with open(sample_dir + '/running_pipeline.txt','w') as f:
#            pass
#        
#        samples_df = pd.read_csv(config['samples'],sep='\t').set_index('Sample_name')
#        samples = samples_df.index.tolist()
        #storage.store('samples',samples)
#        for sample in samples:
#            sample_info = samples_df.loc[sample,['R1','R2']]
#            for f in sample_info[~pd.isna(sample_info)]:
#                assert os.path.exists(sample_dir + '/' + f), (f'The file {f} associated with {sample} in "samples.tsv" is not \
#present in the specified directory ({sample_dir}). Please fix before proceeding. Exiting...')

#            storage.store(sample,[sample_info[0],sample_info[1]])

rule GetSampleInfo:
    output:
        temp(sample_dir + '/running_pipeline.txt')
    run:
        samples_df = pd.read_csv(config['samples'],sep='\t').set_index('Sample_name')
        samples = samples_df.index.tolist()
        for sample in samples:
            sample_info = samples_df.loc[sample,['R1','R2']]
            for f in sample_info[~pd.isna(sample_info)]:
                assert os.path.exists(sample_dir + '/' + f), (f'The file {f} associated with {sample} in "samples.tsv" is not \
present in the specified directory ({sample_dir}). Please fix before proceeding. Exiting...')

        samples_df.to_csv(sample_dir + '/running_pipeline.txt',sep = '\t')

def get_resource_file(ref_dir,type=None):
    ref_exts = tuple(['.fasta','.fa'])
    ann_exts = tuple(['.gff','.gff3','.gtf'])
    for f in os.listdir(ref_dir):
        if f.endswith(ref_exts) and type == 'fasta':
            return f
        elif f.endswith(ann_exts) and type == 'gtf':
            return f
        elif f.endswith('.gz') and type == 'fasta':
            sys.exit('Reference file is compressed. Decompress file before restarting pipeline. Exiting...')
        elif f.endswith('.gz') and type == 'gtf':
            sys.exit('Annotation file is compressed. Decompress file before restarting pipeline. Exiting...')
        else:
            continue
    
    if type == 'fasta':
        sys.exit('Reference fasta file not detected. Ensure file extension is either of ".fa" or ".fasta" and \
exists within the specified reference directory. Exiting...')
    elif type == 'gtf':
        sys.exit('Gene annotation file not detected. Ensure file extension is either of ".gff(3)" or ".gtf" and \
exists within the specified reference directory. Exiting...')

def get_callers():
    callers_set = []
    if config['use_mutect2']:
        callers_set.append('Mutect2')
    if config['use_freebayes']:
        callers_set.append('Freebayes')
    if config['use_bcftools']:
        callers_set.append('Bcftools')
    
    return callers_set

def is_paired(wildcards,R1,R2):
    if pd.isna(R1) | pd.isna(R2):
        return False
    else:
        return True

#def get_aligner_input(wildcards,aligner):
#    target_dir = f'../../{sample_dir}' if aligner == 'bwa' else f'../../../{sample_dir}'
#    R1,R2 = storage.fetch(wildcards.sample)
#    if is_paired(wildcards,R1,R2):
#        if aligner == 'star' or aligner == 'bwa':
#            return f"{target_dir}/{R1}",f"{target_dir}/{R2}"
#        else:
#            return f"-1 {target_dir}/{R1} -2 {target_dir}/{R2}"
#    else:
#        if aligner == 'star' or aligner == 'bwa':
#            return f"{target_dir}/{R1}"
#        else:
#            return f"-U {target_dir}/{R1}"

def get_aligner_input(wildcards,aligner):
    target_dir = f'../../{sample_dir}' if aligner == 'bwa' else f'../../../{sample_dir}'
    samples = pd.read_csv(sample_dir + '/running_pipeline.txt',sep = '\t',index_col=0)
    R1,R2 = samples.loc[wildcards.sample,'R1'],samples.loc[wildcards.sample,'R2']
    if is_paired(wildcards,R1,R2):
        if aligner == 'star' or aligner == 'bwa':
            return f"{target_dir}/{R1}",f"{target_dir}/{R2}"
        else:
            return f"-1 {target_dir}/{R1} -2 {target_dir}/{R2}"
    else:
        if aligner == 'star' or aligner == 'bwa':
            return f"{target_dir}/{R1}"
        else:
            return f"-U {target_dir}/{R1}"

def get_star_readfile_command(wildcards):
    R1 = storage.fetch(wildcards.sample)[0]
    if R1.endswith('.gz'):
        return '--readFilesCommand zcat'
    else:
        return ''

def get_aligner_directory(wildcards):
    prefix = '../results/' + wildcards.sample + '/' + wildcards.aligner + '/'
    if mode == 'RNA':
        return prefix + wildcards.aligner + '_output/'
    else:
        return prefix

def get_input_bam(wildcards,calling=False,ind=False):
    ext = '.bai' if ind else '.bam'
    prefix = '../results/' + wildcards.sample + '/' + wildcards.aligner
    suffix = wildcards.sample + ext
    if calling and config['apply_bqsr']:
        return prefix + '/bqsr/recal_reads_' + suffix
    elif mode == 'RNA':
        return prefix + '/bam_preprocessing/dedup_rna_reads_' + suffix 
    else:
        return prefix + '/bam_preprocessing/dedup_reads_' + suffix

#def get_intersection_output(wildcards):
    #if config['output_combined']:
        #return '-w1 -o ../results/' + wildcards.sample + '/final_calls/final_calls_' + wildcards.sample + '.vcf'
    #else:
        #return '-p ../results/' + wildcards.sample + '/final_calls'
