#!/usr/bin/python

# USAGE: py comborate.py -a jessica/donor_jessica.txt -r jessica/response_jessica.txt -p jessica/iedb_netmhciipan_el_result.csv -o jessica/output

import sys
import os
import pandas as pd
from os.path import join
import argparse

# Adds the src directory to the python path
src_dir = join(os.getcwd(), 'src')
sys.path.append(src_dir)

import rate


# ---------------------------------------------------------------------------
def make_summary_file(comboratefiles):

	# Reads all files
	dfs = [pd.read_csv(f) for f in comboratefiles]

	# Filters out empty dataframes
	dfs = [df for df in dfs	if not df.empty and not df.isna().all().all()]
	
	# If dataframes are not all empty...
	if dfs:
		# ...concatenate and return
		return pd.concat(dfs).sort_values(['Peptide', 'Fisher_pval'])
	else: # If dataframes are all empty, return empty df
		return pd.DataFrame()

# -----------------------------------------------------------------------------
def create_summary_files(outpath):

	# outputdir = join(outpath,'summary')
	cutoff_dirs = [join(outpath,cd) for cd in os.listdir(outpath) if 'reads_cutoff_' in cd]
	
	for cd in cutoff_dirs:

		# Gets the cutoff value
		cutoff = cd.split('reads_cutoff_')[1].split('.')[0]
		
		# Reads all RATE files of cutoff dir
		cutoff_files = [join(cd, cf) for cf in os.listdir(cd) if 'rate' in cf]
		
		# Concatenate all RATE files
		df_sum = make_summary_file(cutoff_files)
		
		# Saves summary file
		df_sum.to_csv(join(outpath, 'summary_reads_cutoff' + cutoff + '.csv'), index=None)
	
	return 
	

# -----------------------------------------------------------------------------
def rm_neg_rf(df):
	'''
		Remove negative correlations from 
		RATE df output
	'''
	return df[df['Relative_freq'] > 0]

# -----------------------------------------------------------------------------
def rm_non_binding_hla(df, binders:list):

	# Removes HLA- in case they still contain this prefix
	binders = [hla.replace('HLA-', '') for hla in binders]

	return df[df['HLA'].isin(binders)]

# -----------------------------------------------------------------------------
def hlas_to_group(df_out, amount:int):

	'''
		Groups the arg:amount HLAs with 
		the best p-value in RATE		
	'''

	# Discards previous groups (+ sign is a tell)
	df_output = df_out[df_out['HLA'].apply(lambda x: '+' not in x)]
	
	pvals = df_output['Fisher_pval']
	
	# Get rows with the arg:amount smallest p-values
	min_p_row = df_output[pvals.isin(pvals.nsmallest(amount))]

	# Returns list of HLAs with smallest p-values
	return min_p_row['HLA'].to_list()

# -----------------------------------------------------------------------------
def parse_prediction_file(binder_file, peptide:str=None, sep=','):

	'''
		Parse HLA prediction file (eg: netMHCpan)
		Filters the HLAs that bind to arg:peptide
		arg:sep = modify according to format (.csv, .tsv, etc)
	'''

	df = pd.read_csv(binder_file, sep=sep)
	if peptide is not None:
		df = df[df['peptide'] == peptide]

	return df

# -----------------------------------------------------------------------------
def get_predicted_hla_binders(df_pred, cutoff:float):

	'''
		Filters HLA predictions for binding 
		HLAs (ie, < arg:cutoff)
	'''

	binders = list(set(df_pred[df_pred['rank'] <= cutoff]['allele'].to_list()))
	
	return [hla.replace('HLA-', '') for hla in binders]

# -----------------------------------------------------------------------------
def make_hla_group_name(hlas_grouped:list, prefix:str=''):

	'''
		Strings all hlas contained in list
		hlas_grouped into a single '+'-separated
		name for the group.
		The final group name starts with 
		arg:prefix
	'''

	if prefix:
		prefix += '-' 
	return prefix + '+'.join(hlas_grouped)

# -----------------------------------------------------------------------------
def make_hla_bin_group(df_hla, hla_grouped:list, hla_group_name:str):

	'''
		Creates a new row of df_hla_bin matrix 
		with the new group name. The donors that express 
		the HLA	contained in the group will be 1s in 
		this new row
	'''
	
	df_hla_cp = df_hla.copy()

	# Identifies donors with 1 in binary matrix (ie, express HLAs in group)
	cols_with_ones = df_hla_cp.loc[df_hla_cp.index.isin(hla_grouped)].any(axis=0)
	
	# Make a new row with new HLA-group name. Donors will be the same as above
	df_hla_cp.loc[hla_group_name] = cols_with_ones.astype(int)
	
	# Returns only the new row
	return df_hla_cp.loc[df_hla_cp.index == hla_group_name]

# -----------------------------------------------------------------------------
def rate_execution(df_res_bin, df_hla_bin, binders_list:list):

	'''
		Runs a full rate execution using the 
		response and HLA information (binary matrix)
		Binders are filtered 
	'''

	# Compute A+R+, ..., A-R-
	df_rate = rate.compute_AR_scores(df_res_bin, df_hla_bin)

	# Compute RF, Odds, Fisher p-value
	df_rate = rate.compute_stats(df_rate)
	
	# Removes rows with negative relative frequency
	df_rate = rm_neg_rf(df_rate)

	# Removes non-binding HLAs (ie, not predicted to bind and not grouped)
	if binders_list:
		df_rate = rm_non_binding_hla(df_rate, binders_list)

	return df_rate

# -----------------------------------------------------------------------------
def write_output_files(df_rate_original, df_pred, output_d, pepnum, response_cutoff):
	
	'''
		Outputs all information to two files
		_rate.csv and _pred.csv with RATE results
		and HLA prediction, respectively
	'''

	# Add a column with the cutoff
	response_cutoff = str(response_cutoff)
	df_rate_original['Response_cutoff'] = response_cutoff

	df_rate_original = df_rate_original.sort_values('Fisher_pval')
	df_rate_original.to_csv(join(output_d, pepnum + '_rate.csv'), index=None)
	
	if not df_pred.empty:
		df_pred.to_csv(join(output_d, pepnum +  '_pred.csv'), index=None)
	
	return

# -----------------------------------------------------------------------------
def parse_cmd_line():

	cmd_parse = argparse.ArgumentParser(description='Iteratively group HLAs to lower RATE p-value')

	# Mandatory
	cmd_parse.add_argument('-a', '--allele-file',  required=True, type=str, help='Allele file name')
	cmd_parse.add_argument('-r', '--response-file', required=True, type=str, help='Response file name')
	cmd_parse.add_argument('-o', '--output', required=True, type=str, help='Output dir name')
	
	# Optional
	cmd_parse.add_argument('-p', '--prediction', required=False, type=str, help='HLA prediction file name', default=False)
	cmd_parse.add_argument('-c', '--rank-cutoff', required=False, type=float, help='Cutoff for the binding rank (default = 25)', default=25)
	cmd_parse.add_argument('-k', '--response-cutoff', required=False, type=str, help='Cutoff for the response (default = 1)', default=1)

	return cmd_parse.parse_args()
# -----------------------------------------------------------------------------

if __name__ == '__main__':

	arg = parse_cmd_line()

	# Command line arguments
	hla_file = arg.allele_file
	response = arg.response_file
	predfile = arg.prediction
	output_d = arg.output
	rank_cut = arg.rank_cutoff 
	response_cutoffs = arg.response_cutoff

	# Create output path
	os.makedirs(output_d, exist_ok=True)

	# Parse input files
	df_res_orig = pd.read_csv(response, sep='\t').dropna(how='all')
	df_hla_orig = pd.read_csv(hla_file, sep='\t')

	# Iterate over each peptide
	pepmax = int(df_res_orig['Peptide #'].max())

	response_cutoffs = [float(c) for c in response_cutoffs.split(',')]

	# Convert the allele matrix into a 1/0 matrix (allele_input)
	df_hla_bin = rate.hla2bin(df_hla_orig)	

	for response_cutoff in response_cutoffs:
		
		# Creates a separate dir for each response cutoff
		output_dir_cutoff = join(output_d, 'reads_cutoff_'+ str(response_cutoff))
		os.makedirs(output_dir_cutoff, exist_ok=True)

		for _, resrow in df_res_orig.iterrows():

			# Isolates one pipeline
			pepseq = resrow['Peptide_Seq']
			pepnum = str(int(resrow['Peptide #'])).zfill(len(str(pepmax)))
			
			## - Running the pipeline for the peptide isolated above

			# Filtering the responses for the peptide selected
			df_res = df_res_orig[df_res_orig['Peptide_Seq'] == pepseq]

			# Filter HLAs from binding prediction file
			if predfile:
				# List of HLA binders for the selected peptide
				df_pred = parse_prediction_file(predfile, peptide=pepseq)

				# Get list of HLA predicted as binders
				list_of_binders = get_predicted_hla_binders(df_pred, cutoff = rank_cut)

			else: 
				list_of_binders = []
				df_pred = pd.DataFrame()

			# Convert the response input to a 1/0 matrix (response_input)
			df_res_bin = rate.res2bin(df_res, response_cutoff)

			# Iteratively grouping the HLAs 
			iteration = 1
			hlas_grouped = list()

			# Run full rate execution: 0th iteration, original 
			df_rate_original = rate_execution(df_res_bin, df_hla_bin, list_of_binders)
			df_rate = df_rate_original
			
			# In case RATE results are empty, create empty output files
			if df_rate.empty:
				write_output_files(df_rate_original, df_pred, output_dir_cutoff, pepnum, response_cutoff)
				continue

			# Groups HLAs while p-values decrease
			while True:
				
				# Get N best p-values
				hlas_grouped = hlas_to_group(df_rate, amount=1+iteration)

				# HLA-group name for this iteration
				hla_group_name = make_hla_group_name(hlas_grouped, str(iteration))

				# Update df_hla_bin
				df_hla_bin_group = make_hla_bin_group(df_hla_bin, hlas_grouped, hla_group_name)
				
				if predfile:
					# Update list of binders with group name 
					list_of_binders.append(hla_group_name)

				# Get RATE results
				df_rate = rate_execution(df_res_bin, df_hla_bin_group, list_of_binders)

				# Stops while-loop if p-value does not decrease
				if df_rate['Fisher_pval'].min(axis=0) >= df_rate_original['Fisher_pval'].min(axis=0):
					break

				# Appends group result to original RATE results
				df_rate_original = pd.concat([df_rate_original, df_rate]).sort_values('Fisher_pval').reset_index(drop=True)

				df_rate = df_rate_original

				iteration += 1

			# Writes output files for this peptide
			write_output_files(df_rate_original, df_pred, output_dir_cutoff, pepnum, response_cutoff)


create_summary_files(output_d)

