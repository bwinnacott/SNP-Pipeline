import pandas as pd
import os
import re

def get_samples():
    

def get_resource_file(ref_dir,type=None):
    ref_exts = tuple(['.fasta','.fa','.fasta.gz','.fa.gz'])
    ann_exts = tuple(['.gff','.gff3','.gtf'])
    for f in os.listdir(ref_dir):
        if f.endswith(ref_exts) and type == 'fasta':
            return f
        elif f.endswith(ann_exts) and type == 'gtf':
            return f
        else:
            continue
    
    if type == 'fasta':
        sys.exit('Reference fasta file not detected. Ensure file extension is either of ".fa" or ".fasta". Exiting...')
    elif type == 'gtf':
        sys.exit('Gene annotation file not detected. Ensure file extension is either of ".gff(3)" or ".gtf". Exiting...')

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

def get_callers():
    callers_set = []
    if config['use_mutect2']:
        callers_set.append('Mutect2')
    if config['use_freebayes']:
        callers_set.append('Freebayes')
    
    return callers_set

def get_aligner_input(wildcards,aligner):
    prefix = sample_dir + '/' + wildcards.sample
    if samples.count(wildcards.sample) == 2:
        sample = set([f for f in os.listdir(sample_dir) if re.split('_R[1-2]',f,maxsplit=1)[0] == wildcards.sample])
        for ext in exts:
            if os.path.exists(prefix + '_R1.' + ext):
                R1 = os.path.abspath(prefix + '_R1.' + ext)
                R2 = os.path.abspath(prefix + '_R2.' + ext)
                if aligner == 'star' or aligner == 'bwa':
                    return R1,R2
                else:
                    return f"-1 {R1} -2 {R2}"
    else:
        for ext in exts:
            if os.path.exists(prefix + '.' + ext):
                if aligner == 'star' or aligner == 'bwa':
                    return os.path.abspath(prefix + '.' + ext)
                else:
                    return f"-U {os.path.abspath(prefix + '.' + ext)}"

def get_star_readfile_command(wildcards):
    prefix = sample_dir + '/' + wildcards.sample
    for ext in exts:
        if os.path.exists(prefix + '_R1.' + ext):
            if ext.endswith('.gz'):
                return '--readFilesCommand zcat'
            else:
                return ''
        elif os.path.exists(prefix + '.' + ext):
            if ext.endswith('.gz'):
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

def get_intersection_output(wildcards):
    if config['output_combined']:
        return '-w1 -o ../results/' + wildcards.sample + '/final_calls/final_calls.vcf'
    else:
        return '-p ../results/' + wildcards.sample + '/final_calls'