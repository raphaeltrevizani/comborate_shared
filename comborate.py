#!/usr/bin/python

# USAGE: py comborate.py -a jessica/donor_jessica.txt -r jessica/response_jessica.txt -p jessica/iedb_netmhciipan_el_result.csv -o jessica/output

import sys
import os
import pandas as pd
from os.path import join
import argparse
import combineAB
import subprocess
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
def write_output_files(dfrate, df_pred, output_d, pepnum, response_cutoff, relfreq_cutoff, pval_cutoff):
	
	'''
		Outputs all information to two files
		_rate.csv and _pred.csv with RATE results
		and HLA prediction, respectively
	'''

	# Add a column with the cutoff
	response_cutoff = str(response_cutoff)
	dfrate['Response_cutoff'] = response_cutoff

	dfrate = dfrate.sort_values('Fisher_pval')

	dfrate = dfrate[dfrate['Relative_freq'] >= relfreq_cutoff]
	dfrate = dfrate[dfrate['Fisher_pval'] <= pval_cutoff]

	dfrate.to_csv(join(output_d, pepnum + '_rate.csv'), index=None)
	
	if not df_pred.empty:
		df_pred.to_csv(join(output_d, pepnum +  '_pred.csv'), index=None)
	
	return

# -----------------------------------------------------------------------------
def get_hla_class(hla):

	'''
		Returns the class of arg:hla
	'''

	hla = hla.replace('HLA-', '')

	if 'DP' in hla or 'DQ' in hla or 'DR' in hla:
		return 'II'
	elif 'A' or 'B' or 'C' in hla: 
		return 'I'

# -----------------------------------------------------------------------------
def group_by_hla_class(hla_list):
	
	grouped_hlas = {'I':[], 'II':[]}

	for hla in hla_list:
		grouped_hlas[get_hla_class(hla)] += [hla]

	return grouped_hlas

# -----------------------------------------------------------------------------
def parse_allowed_alleles():
	with open('alleles_name.list') as inp:
		next(inp)
		data = inp.readlines()
	data = [line.replace('\n', ',').replace('\t', ',').replace(' ', '').replace(',,', ',').split(',') for line in data]
	data = [item for sublist in data for item in sublist if item != '']
	return data

# -----------------------------------------------------------------------------
def submit_cmd_to_API_classI(peptides, hlas, sizes):

	method = 'netmhcpan_el'
	mhc_class = 'i'

	command = "curl --data \"method=" + method + "&sequence_text="+peptides+"&allele=" + hlas + "&length="+ sizes +"\" https://tools-cluster-interface.iedb.org/tools_api/mhc"+mhc_class+"/"
	result = subprocess.run(command, shell=True, capture_output=True, text=True)
	
	output = result.stdout
	errors = result.stderr

	return output

# -----------------------------------------------------------------------------
def submit_cmd_to_API_classII(peptides, hlas, sizes):

	method = 'netmhciipan_el'
	mhc_class = 'ii'

	command = "curl --data \"method=" + method + "&sequence_text="+peptides+"&allele=" + hlas + "&length="+ sizes +"\" https://tools-cluster-interface.iedb.org/tools_api/mhc"+mhc_class+"/"
	result = subprocess.run(command, shell=True, capture_output=True, text=True)

	output = result.stdout
	errors = result.stderr

	return output

# -----------------------------------------------------------------------------
def parse_valid_hla_file():
	with open('valid_hlas.txt') as inputfile:
		return inputfile.read().splitlines()
	
# -----------------------------------------------------------------------------
def run_API_class_I(peptides, hlalist, outputname):

	# Add flanking characters needed by the API query
	peptides = ''.join(['%3Epeptide' + str(num) + '%0A' + pep.rstrip() + '%0A' for num, pep in enumerate(peptides, start = 1)])

	valid_hlas = parse_valid_hla_file()

	mod_hlalist = ['HLA-' + hla for hla in list(set(hlalist))]
	mod_hlalist = [hla for hla in mod_hlalist if hla in valid_hlas]

	hlas = ','.join(mod_hlalist)
	sizes = ','.join('9' * len(mod_hlalist))

	output = submit_cmd_to_API_classI(peptides, hlas, sizes)

	with open(outputname, 'w') as out:
		out.write(output.replace('\t', ',').replace('percentile_rank', 'rank'))
	
	return

# -----------------------------------------------------------------------------
def run_API_class_II(peptides, hlalist, outputname):

	# Add flanking characters needed by the API query
	peptides = ''.join(['%3Epeptide' + str(num) + '%0A' + pep.rstrip() + '%0A' for num, pep in enumerate(peptides, start = 1)])

	# TODO group peptide by sizes and run each size separate

	sizes = '15'
	header = '\t'.join(['allele', 'seq_num', 'start', 'end', 'length', 'core_peptide', 'peptide', 'score', 'rank'])
	concat_output = header + '\n' 

	for hla in hlalist:
		output = submit_cmd_to_API_classII(peptides, hla, sizes)

		if 'Invalid' not in output:
			concat_output += '\n'.join(output.split('\n')[1:])

	with open(outputname, 'w') as out:
		out.write(concat_output.replace('\t', ','))

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
	cmd_parse.add_argument('-q', '--auto-prediction', required=False, action='store_true', help='Automatically run HLA prediction', default=False)
	cmd_parse.add_argument('-c', '--rank-cutoff', required=False, type=float, help='Cutoff for the binding rank (default = 25)', default=25)
	cmd_parse.add_argument('-k', '--response-cutoff', required=False, type=str, help='Cutoff for the response (default = 1)', default='1')
	cmd_parse.add_argument('-v', '--p-value', required=False, type=float, help='Cutoff for the p-value (default = 0.05)', default=0.05)
	cmd_parse.add_argument('-f', '--rel-freq', required=False, type=float, help='Cutoff for the Relative Frequency (default = 1)', default=1)
	cmd_parse.add_argument('-x', '--combine-AB-chains', required=False, action='store_true', help='Combine alpha and beta chains', default=False)


	return cmd_parse.parse_args()
# -----------------------------------------------------------------------------
def split_df_hla_class(df):

	# Get HLAs from df
	flat_list_hlas = [item for sublist in df.values for item in sublist if not pd.isna(item)]

	# ID the HLA class
	hla_class = set([get_hla_class(hla) for hla in flat_list_hlas])

	if len(hla_class) == 2:
		mask_abc = df.apply(lambda row: row.astype(str).str.contains(r'A\*|B\*|C\*').any(), axis=1)
		mask_d   = df.apply(lambda row: row.astype(str).str.contains(r'^D').any(), axis=1)

		df_mhci = df[mask_abc]
		df_mhcii   = df[mask_d]	
		
		return [df_mhci, df_mhcii]
	else:
		return [df]
# -----------------------------------------------------------------------------
def get_hla_id(df):
	# Get HLAs from df
	flat_list_hlas = [item for sublist in df.values for item in sublist if not pd.isna(item)]

	# ID the HLA class
	hla_class = list(set([get_hla_class(hla) for hla in flat_list_hlas]))[0]

	return hla_class
# -----------------------------------------------------------------------------
if __name__ == '__main__':

	arg = parse_cmd_line()

	# Command line arguments
	hla_file = arg.allele_file
	response = arg.response_file
	predfile = arg.prediction
	rank_cut = arg.rank_cutoff 
	rf_cutoff = arg.rel_freq 
	p_cutoff = arg.p_value
	autopred = arg.auto_prediction
	combine_chains = arg.combine_AB_chains


	# Create output path
	os.makedirs(arg.output, exist_ok=True)

	# Parse input files
	df_res_orig = pd.read_csv(response, sep='\t').dropna(how='all')
	df_hla_orig = pd.read_csv(hla_file, sep='\t', dtype=str)

	list_dfs_hla = split_df_hla_class(df_hla_orig)

	tmp_predfile = ''

	for df_hla_orig in list_dfs_hla:

		hla_class = get_hla_id(df_hla_orig)
		
		output_d = join(arg.output, 'class_' + hla_class)

		if autopred:
			## Running API prediction
			
			flat_list_hlas = [item for sublist in df_hla_orig.values for item in sublist if not pd.isna(item)]
			
			# Get peptides from df
			peptides = df_res_orig['Peptide_Seq'].drop_duplicates().to_list()
			
			if hla_class == 'II':
				tmp_predfile = 'prediction_MHC_II.tsv'
				run_API_class_II(peptides, flat_list_hlas, tmp_predfile)
				predfile = tmp_predfile
			elif hla_class == 'I':
				tmp_predfile = 'prediction_MHC_I.tsv'
				run_API_class_I(peptides, flat_list_hlas, tmp_predfile)
				predfile = tmp_predfile


		if combine_chains:
			df_hla_orig = combineAB.make_AB_combinations(df_hla_orig)

		# Iterate over each peptide
		pepmax = int(df_res_orig['Peptide #'].max())

		response_cutoffs = [float(c) for c in arg.response_cutoff.split(',')]

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
				df_rate_original = df_rate_original[df_rate_original['Relative_freq'] >= rf_cutoff]
				df_rate = df_rate_original

				# In case RATE results are empty, create empty output files
				if df_rate.empty:
					write_output_files(df_rate_original, df_pred, output_dir_cutoff, pepnum, response_cutoff, rf_cutoff, p_cutoff)
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
				write_output_files(df_rate_original, df_pred, output_dir_cutoff, pepnum, response_cutoff, rf_cutoff, p_cutoff)


		create_summary_files(output_d)

		if tmp_predfile:
			os.remove(tmp_predfile)