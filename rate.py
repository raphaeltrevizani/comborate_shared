#!/usr/bin/python

import pandas as pd
import sys
from copy import deepcopy
from scipy.stats import fisher_exact
import numpy as np

# -----------------------------------------------------------------------------
def res2bin(df_res_txt, cutoff):

	# Set values below cutoff to 0
	cols = [col for col in df_res_txt.columns if 'Donor' in col]
	dx = df_res_txt[cols]
	dx = dx.astype('Int64')

	# Set values above 0 to 1
	df_res_bin = deepcopy(df_res_txt)
	df_res_bin[cols] = dx.ne(0).astype('boolean').mask(dx.isna()).astype('Int64')

	# Tidy up
	df_res_bin = df_res_bin.drop(columns=['Peptide #', 'Peptide_ID'])
	df_res_bin = df_res_bin.set_index("Peptide_Seq")
	df_res_bin.index.name = None
	df_res_bin.columns.name = None

	# Apply cutoff
	df_res_bin[df_res_bin < cutoff] = 0
	
	return df_res_bin
# -----------------------------------------------------------------------------
def hla2bin(df_hla_txt):
	
	nullvals = ['n/a', 'N/A', '-', 'NA', 'na']
	df_hla_txt = df_hla_txt.replace(nullvals, np.nan)
	df_hla_txt = df_hla_txt.replace('HLA-', '', regex=True)
	df_tmp = df_hla_txt.melt(var_name="Donor", value_name="HLA")
	df_hla_bin = (pd.crosstab(df_tmp["HLA"], df_tmp["Donor"]) > 0).astype(int)
	df_hla_bin.index.name = None
	df_hla_bin.columns.name = None

	# step 2.B: check if all loci are typed for a donor. NaN if not
	loci_required = {'DR', 'DQ', 'DP'}
	for donor in df_hla_bin.columns:
		loci_present = df_hla_txt[donor].str[:2]
		missing_loci = loci_required.difference(set(loci_present.to_list()))
		if missing_loci:
			for missing_locus in missing_loci:
				missing_hla = df_hla_bin.index.str.contains(missing_locus)
				df_hla_bin.loc[missing_hla, donor] = np.nan

	return df_hla_bin

# -----------------------------------------------------------------------------
def compute_AR_scores(df_res, df_hla):

	res = {}
	for hla, hla_row in df_hla.iterrows():
		for pep, res_row in df_res.iterrows():
			valid = hla_row.notna() & res_row.notna()
			x, y = hla_row[valid], res_row[valid]
			res[hla+'_'+pep] = {
				'A+R+': ((x>=1) & (y>=1)).sum(),
				'A-R+': ((x==0) & (y>=1)).sum(),
				'A+R-': ((x>=1) & (y==0)).sum(),
				'A-R-': ((x==0) & (y==0)).sum(),
			}

	return pd.DataFrame.from_dict(res, orient='index')

# -----------------------------------------------------------------------------
def compute_stats(df):

	df["Fisher_pval"] = df.apply(
		lambda r: fisher_exact(r.to_numpy().reshape(2,2))[1],
		axis=1)

	df['Odds_ratio'] = df['A+R+'] * df['A-R-']/(df['A-R+'] * df['A+R-'])

	df['N_donors'] = df[['A+R+', 'A-R+', 'A+R-', 'A-R-']].sum(axis=1)
	df = df.sort_values('A+R+', ascending=False)
	if 'key_0' in df.columns:
		del df['key_0']

	df['Relative_freq'] = (df['A+R+']/(df['A+R+'] + df['A+R-']))/((df['A+R+'] + df['A-R+'])/df['N_donors'])

	# Formatting output
	df[['HLA', 'Peptide']] = df.index.to_series().str.split('_', expand=True)
	df = df.reset_index(drop=True)
	df = df[['HLA', 'Peptide','A+R+', 'A-R+', 'A+R-', 'A-R-', 'N_donors', 'Relative_freq', 'Odds_ratio', 'Fisher_pval']]
	return df
# -----------------------------------------------------------------------------
def run(df_hla_orig, df_res_orig, response_cutoff):

	# step 1: convert the response input to a 1/0 matrix (response_input)
	df_res = res2bin(df_res_orig, response_cutoff)

	# step 2: convert the allele matrix into a 1/0 matrix (allele_input)
	df_hla = hla2bin(df_hla_orig)

	# step 3: compute RATE numbers
	df = compute_AR_scores(df_res, df_hla)
	df = compute_stats(df)

	return df
# -----------------------------------------------------------------------------


if __name__ == '__main__':
	hla_file = sys.argv[1] 
	response = sys.argv[2]
	output_f = sys.argv[3]

	df_res_orig = pd.read_csv(response, sep='\t')
	df_hla_orig = pd.read_csv(hla_file, sep='\t')

	df = run(df_hla_orig, df_res_orig, response_cutoff =  20)

	df.to_csv(output_f, index=None)