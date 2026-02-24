#!/usr/bin/python

# python count_donors.py jessica/output_cutoff25/ jessica/donor_jessica.txt jessica/response_jessica.txt

import pandas as pd
import sys
import os
from collections import Counter
import math

# ------------------------------------
def get_most_frq_hla(df_dnr, donors):
	hlas_count = dict(get_hlas_for_donors(df_dnr, donors))
	hlas_count = {k: v for k, v in hlas_count.items() if not (isinstance(k, float) and math.isnan(k))}
	max_hlas = {k: v for k, v in hlas_count.items() if v == max(hlas_count.values())}
	return ';'.join(max_hlas), max(hlas_count.values())

# ------------------------------------
def get_hlas_for_donors(df_dnr, donors):
	donors = list(donors)
	all_hlas = list()
	for donor in donors:
		all_hlas += list(set(df_dnr[donor].to_list()))
	return Counter(all_hlas)

# ------------------------------------
def get_hlas_from_rate_output(dfrate):

	hlas_list = dfrate['HLA'].str.split(r'\+').to_list()
	flat_list = [item for sublist in hlas_list for item in sublist]

	cleanlist = list()
	for item in flat_list:
		if '-' in item:
			item = item.split('-')[1]
		cleanlist.append(item)

	return list(set(cleanlist))

# ------------------------------------
def get_donors_w_hla(df_hla_donor, hla_list, positive_donors):
	donors_w_HLAs = list()
	
	hlas_list = set(hla_list)

	donors_w_HLAs = [c for c in df_hla_donor.columns if df_hla_donor[c].astype(str).isin(hlas_list).any()]

	donors_in_common = set(donors_w_HLAs).intersection(set(positive_donors))

	return donors_in_common

# ------------------------------------
def count_donors_w_hla(df_hla_donor, hla_list, positive_donors):
	return len(get_donors_w_hla(df_hla_donor, hla_list, positive_donors))

# ------------------------------------
def get_positive_donors(dfrate, dfresp):
	peptide = dfrate['Peptide'].head(1).values[0]
	df_pep = dfresp[dfresp['Peptide_Seq'] == peptide]
	pos_donors = df_pep.iloc[:, 3:].columns[(df_pep.iloc[:, 3:] > 0).any()]

	return set(pos_donors)

# ------------------------------------

def run(rateoutputdir, allele_inputf, responseinput):

	ratefiles = sorted([f for f in os.listdir(rateoutputdir) if 'rate' in f])

	df_alle = pd.read_csv(allele_inputf, sep='\t')
	df_resp = pd.read_csv(responseinput, sep='\t')
	
	dict_count_dnr = dict()
	for ratefile in ratefiles:

		inputfile = os.path.join(rateoutputdir, ratefile)

		df_rate = pd.read_csv(inputfile)
		
		# Skips iteration if there are no RATE results
		if df_rate.empty:
			continue
		
		pos_donors = get_positive_donors(df_rate, df_resp)
		hla_max, max_frq = get_most_frq_hla(df_alle, pos_donors)

		# Filtering Relative Frequency > 1 (positive associations)
		df_posRF = df_rate[df_rate['Relative_freq'] >= 1]

		# Filtering significant values
		df_sig = df_posRF[df_posRF['Fisher_pval'] <= 0.05]

		# Filtering individual HLAs
		df_ind = df_sig[~df_sig['HLA'].str.contains(r'\+')]

		# Filtering groups of HLAs
		df_grp = df_sig[df_sig['HLA'].str.contains(r'\+')]

		# Getting HLAs from RATE results
		hla_indvd = get_hlas_from_rate_output(df_ind)
		hla_group = get_hlas_from_rate_output(df_grp)


		hla_total = list(set(hla_indvd + hla_group))
		n_donor_ind = count_donors_w_hla(df_alle, hla_indvd, pos_donors)
		n_donor_grp = count_donors_w_hla(df_alle, hla_total, pos_donors)

		peptide_number = ratefile.split('_')[0]

		total_donors = len(pos_donors)

		hla_restrictions_gained = list(set(hla_group).difference(hla_indvd))
		
		# dict_count_dnr[peptide_number] = {
		# 	'Positives donors':f'{total_donors}',
		# 	'Individual HLAs': f'{n_donor_ind}',
		# 	'Frequency of individual HLAs':f'{n_donor_ind/total_donors}',
		# 	'RATE restrictions':f'{';'.join(hla_indvd)}',
		# 	'Grouped HLAs': f'{n_donor_grp}',
		# 	'Frequency of grouped HLAs':f'{n_donor_grp/total_donors}',
		# 	'No allele restrictions gained':f'{len(hla_restrictions_gained)}',
		# 	'ComboRATE restrictions gained':f'{';'.join(hla_restrictions_gained)}',
		# 	'Occurrences of most frequent HLA':f'{max_frq}',
		# 	'Frequency of most freq. HLA':f'{max_frq/total_donors}',
		# 	'HLA with highest frequency':f'{hla_max}'
		# 	}
		
		dict_count_dnr[peptide_number] = {
			'Sequence' : df_rate['Peptide'].iloc[0],
			'Positives donors':f'{total_donors}',
			'# restrictions':f'{len(hla_total)}',
			'ComboRATE restrictions':f'{';'.join(hla_total)}'
			}


	df_out = pd.DataFrame.from_dict(dict_count_dnr).transpose()
	df_out.index.name = 'Peptide'
	df_out.to_csv(os.path.join(rateoutputdir, 'summary_restrictions.csv'))

# ------------------------------------
def summarize_donor_count(ratedir):
	
	cutoff_dirs = [os.path.join(ratedir, sd) for sd in os.listdir(ratedir) if 'reads_cutoff_' in sd]
	
	listdf = list()
	for inputdir in cutoff_dirs:
		
		inputfile = os.path.join(inputdir,'summary_restrictions.csv')

		df = pd.read_csv(inputfile)
		df['Reads cutoff'] = inputdir.split('reads_cutoff_')[1]

		listdf.append(df)

	df_final = pd.concat(listdf)
	df_final = df_final[['Peptide', 'Sequence', 'Reads cutoff', 'Positives donors', '# restrictions', 'ComboRATE restrictions']]

	df_final.to_csv(os.path.join(ratedir, 'summary_restrictions.csv'), index=None)

# ------------------------------------

if __name__ == '__main__':

	if len(sys.argv) < 4:
		print('Error. Missing arguments.')
		print('Usage: ')
		print('python count_donors.py <rate_files_dir> <donor_alleles_file.txt> <response_file.txt>')
		exit(0)

	rateoutputdir = sys.argv[1]
	allele_inputf = sys.argv[2]
	responseinput = sys.argv[3]

	subdirs = [os.path.join(rateoutputdir, sd) for sd in os.listdir(rateoutputdir) if 'reads_cutoff_' in sd]

	# Runs analisys in reads_cutoff subdirs
	for sd in subdirs:
		run(sd, allele_inputf, responseinput)

	summarize_donor_count(rateoutputdir)
