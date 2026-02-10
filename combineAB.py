#!/usr/bin/python

from collections import defaultdict
import pandas as pd
import sys
# ------------------------------------------------------------
def make_AB_combinations(df):

	'''
		Creates all chain combinations
		for alleles alpha/beta
	'''

	df_combined = pd.DataFrame()

	for donor in df.columns:

		alleles = [x for x in df[donor].to_list() if pd.notna(x)]

		# Group by first two letters (DP, DQ, DR)
		locus_dict = defaultdict(list)
		for hla in alleles:
			locus_dict[hla[:2]].append(hla)

		# Create A/B combinations per loci
		combined = []
		for _, items in locus_dict.items():
			A = [x for x in items if x[2] == 'A']
			B = [x for x in items if x[2] == 'B']
			if A and B:  # make all A/B pairs
				combined.extend(f"{a}/{b}" for a in A for b in B)
			else:        # keep as-is if no pairing possible
				combined.extend(items)

		# New column can be longer or shorter, so make a pd.Series
		new_col = pd.Series(combined, name=donor).drop_duplicates()

		# Concatenate donor (pd.Series) to df. Index is auto-aligned, NaN filled
		df_combined = pd.concat([df_combined, new_col], axis=1)

	return df_combined

# ------------------------------------------------------------
def parse_hla_donor_file(hlafile:str):

	'''
		arg: HLA_donor input file name
		Parses the HLA/donor file and returns a dataframe 
		formatted as

			1203     1568     1576     1577
		0  A*02:01  A*03:01  A*02:01  A*02:01
		1  A*24:02  A*32:01  A*02:01  A*03:01
		2  B*27:05  B*07:05  B*27:10  B*44:02
		3  B*39:06  B*51:01  B*27:10  B*44:02
		4  C*02:02  C*04:01  C*15:02  C*05:01
		5  C*07:02  C*15:05  C*15:02  C*07:02

	'''

	df = pd.read_csv(hlafile, sep='\t')
	df = df.replace('HLA-', '', regex=True)
	df = df.dropna(axis=1, how='all')

	return df

# ------------------------------------------------------------
def filter_combinations(df_ABcomb, filter_file):
	df_filter = pd.read_csv(filter_file, header=None)
	df_filter = df_filter.replace('HLA-', '', regex=True)
	print(df_filter)
	df_ABcomb


# ------------------------------------------------------------
if __name__ == '__main__':
	df = parse_hla_donor_file(sys.argv[1])
	# print(df)
	df_comb = make_AB_combinations(df)
	df_comb.to_csv(sys.argv[2], index=None, sep='\t')
	# print(df_comb)
	# filter_combinations(df_comb, 'linkage_hlaii.txt')