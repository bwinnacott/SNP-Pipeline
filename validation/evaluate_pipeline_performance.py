#!/usr/bin/env python
import pandas as pd
import argparse
import allel
import numpy as np
import os

def get_args():
	# get command line arguments
	parser = argparse.ArgumentParser(description='Computes performance metrics for simulated SNV datasets post variant calling. \
	Note: when providing list of filenames for vcf and simulation files, the name should be an absolute path if the file(s) are \
	present in another directory.')
	parser.add_argument('-i','--in_file',action='store',type=str,default=None,help='Input vcf file.')
	parser.add_argument('-s','--in_sim',action='store',type=str,default=None,help='Input simulation file. Contains info for simulated (known) variants.')
	parser.add_argument('-o','--out_file',action='store',type=str,help='Output file.')
	parser.add_argument('--caller',action='store',type=str,default="Intersection",help='Specifies variant calling tool used to generate input vcf file(s).')
	parser.add_argument('--metric',action='store',type=str,default='recall',choices=['precision','recall'],help='Specifies metric to calculate (i.e., precision or recall).')
	parser.add_argument('--vaf_low',action='store',type=float,default=0,help='Specifies lower bound of VAF range.')
	parser.add_argument('--vaf_high',action='store',type=float,default=0.5,help='Specifies upper bound of VAF range.')
	parser.add_argument('--step',action='store',type=float,default=0.05,help='Specifies step size for which vaf should increase.')
	parser.add_argument('-I','--in_file_list',action='store',type=str,default=[],help='Input list (txt file) of vcf filenames.')
	parser.add_argument('-S','--in_sim_list',action='store',type=str,default=[],help='Input list (txt file) of simulation filenames. Required \
						if -I is specified. Order of filenames in list must match that of input vcf filename list.')

	return parser.parse_args()

def get_file_lists(vcf_filename_list,sim_filename_list):
	# store filenames when list of files is provided
	vcf_filenames = []
	sim_filenames = []
	# open each file and store filenames in above objects
	with open(vcf_filename_list,'r') as vcf,open(sim_filename_list,'r') as sim:
		for vcf_filename in vcf:
			vcf_filename = vcf_filename.strip()
			vcf_filenames.append(vcf_filename)
		for sim_filename in sim:
			sim_filename = sim_filename.strip()
			sim_filenames.append(sim_filename)
	
	return vcf_filenames,sim_filenames

def extract_vaf(file,caller):
	calls = allel.read_vcf(file,fields="*",alt_number=1)
	if caller.lower() == "mutect2":
		return [AF[0] for AF in calls['calldata/AF']]
	elif caller.lower() == "freebayes":
		DP = calls['calldata/DP']
		AO = calls['calldata/AO']
		return [round(x[0]/y[0],2) for x,y in zip([z for z in AO],[w for w in DP])]
	elif caller.lower() == "bcftools":
		DP4 = calls['variants/DP4']
		return [round((x[2]+x[3])/(x[0]+x[1]+x[2]+x[3]),2) for x in DP4]
	else:
		return [None] * len(calls['calldata/DP'])

def get_simulated_calls(vcf_file,sim_file,caller):
	# read in input vcf file to dataframe
	callset = allel.vcf_to_dataframe(vcf_file,fields="*",alt_number=1)
	# add column with info specific to each variant
	callset['Variant'] = callset['CHROM'] + ',' + callset['POS'].astype(str)
	# read in simulated output file
	sim_info = pd.read_csv(sim_file,header=0,delim_whitespace=True)
	# repeat column addition (different column names to combine)
	sim_info['Variant'] = sim_info['chromosome'] + ',' + sim_info['position.1'].astype(str)
	# merge the callset and simulated variant summary file
	simulated_calls = pd.merge(callset,sim_info,how='inner',on='Variant')
	# add VAF values to original dataframe for precision calculations at VAF ranges
	#callset['output_VAF'] = extract_vaf(vcf_file,caller)

	return callset,sim_info,simulated_calls

def discordant_alt_cnt(sub_calls,check_alt=False):
	# keep track of how many alternate alleles called don't match simulation record at given position
	discordant_alts = 0
	if not check_alt:
		return discordant_alts
	# convert calls in VAF range to dictionary
	calls_dict = sub_calls.to_dict('records')
	# iterate over each called variant and check if called allele matches simulated allele
	for row in calls_dict:
		if row['ALT'] != row['alt_allele']:
			discordant_alts += 1

	return discordant_alts

def get_vaf_counts(variant_df,vaf_low,vaf_high,step_size,calls=False):
	# store count values for each vaf range
	vaf_counts = {}
	# loop over defined vaf range intervals (default 0.05) and get counts
	for x in np.arange(vaf_low,vaf_high,step_size):
		upper_vaf = x + step_size
		vaf_range = variant_df[(variant_df['output_VAF'] > x) & (variant_df['output_VAF'] < upper_vaf)]
		vaf_counts[str(round(x,2))+'-'+str(round(upper_vaf,2))] = len(vaf_range.index) - discordant_alt_cnt(vaf_range,check_alt=calls)
	
	return vaf_counts

def calculate_metric(vaf_counts,vaf_totals):
	# store metric values by vaf range
	metric_by_vaf = {}
	# loop over each range and get calculate metric
	for vaf_range,count in vaf_counts.items():
		if vaf_totals[vaf_range] == 0:
			vaf_totals[vaf_range] == None
		else:
			metric_by_vaf[vaf_range] = round(count/vaf_totals[vaf_range],3)
	
	return metric_by_vaf

def get_metric_df(sim_calls,total_variants,sample_index,vaf_low,vaf_high,step_size):
	# get number of simulated variants called per vaf range
	called_sim_counts = get_vaf_counts(sim_calls,vaf_low,vaf_high,step_size,calls=True)
	# get total number of variants simulated per vaf range
	total_counts = get_vaf_counts(total_variants,vaf_low,vaf_high,step_size)
	print(total_counts)
	# calculate sensitivity for each vaf range
	metric = calculate_metric(called_sim_counts,total_counts)
	# convert to dataframe
	metric_df = pd.DataFrame(metric,index=[sample_index])

	return metric_df

def main():
	# store command line arguments in variable
	args = get_args()
	# if vcf filename list is provided, enter if
	if args.in_file_list:
		# ensure the corresponding variant simulation summary files are specified in list
		assert args.in_sim_list != [], "List of simulated variant summary filenames is absent, despite list of \
input list of vcf filenames present. Provide filename list of simulated variant summaries. Exiting..."
		# get list of filenames for vcf and simulation files
		vcf_filenames,sim_filenames = get_file_lists(args.in_file_list,args.in_sim_list)
		# initiate dataframe to store output
		final_df = pd.DataFrame(columns=[str(round(x,2))+'-'+str(round(x+args.step,2)) for x in np.arange(args.vaf_low,args.vaf_high,args.step)])
		# loop over each pair of filenames
		for ind,(vcf,sim) in enumerate(zip(vcf_filenames,sim_filenames)):
			# if absolute path is provided, update sample name
			if os.path.isabs(vcf):
				file_name = os.path.split(vcf)[1]
				sample_index = os.path.splitext(file_name)[0]
			# else just split the extension
			else:
				sample_index = os.path.splitext(vcf_filenames[ind])[0]
			# get dataframes for original callset, simulation info, and simulated callset
			callset,sim_summary,sim_calls = get_simulated_calls(vcf,sim,args.caller)
			if args.metric == 'recall':
				# get dataframe of sensitivity values for each vaf range
				metric_df = get_metric_df(sim_calls,sim_summary,sample_index,args.vaf_low,args.vaf_high,args.step)
			elif args.metric == 'precision':
				# get dataframe of precision values for each vaf range
				metric_df = get_metric_df(sim_calls,callset,sample_index,args.vaf_low,args.vaf_high,args.step)
			else:
				sys.exit("Metric specified at command line does not match either 'precision' or 'recall'.")
			# update final df
			final_df = final_df.append(metric_df)
		# write out to output file
		final_df.to_csv(args.out_file,sep='\t')
	# else only a single vcf file is provided as input
	else:
		# make sure files are specified
		assert args.in_file != None, ("Missing input vcf file. Exiting...")
		assert args.in_sim != None, ("Missing simulated variant summary information file. Exiting...")
		# get sample name and dataframes for simulated variants
		sample_index = os.path.splitext(args.in_file)
		callset,sim_summary,sim_calls = get_simulated_calls(args.in_file,args.in_sim,args.caller)
		if args.metric == 'recall':
			# get dataframe of sensitivity values for each vaf range
			metric_df = get_metric_df(sim_calls,sim_summary,sample_index,args.vaf_low,args.vaf_high,args.step)
		elif args.metric == 'precision':
			# get dataframe of precision values for each vaf range
			metric_df = get_metric_df(sim_calls,callset,sample_index,args.vaf_low,args.vaf_high,args.step)
		else:
			sys.exit("Metric specified at command line does not match either 'precision' or 'recall'.")
		# write out final sensitivity df to file
		metric_df.to_csv(args.out_file,sep='\t')

if __name__ == "__main__":
	main()