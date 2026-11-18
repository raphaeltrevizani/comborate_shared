#!/usr/bin/python

'''
	Module to standardize HLA names
'''
import re
# ------------------------------------------------------
def standardize_name(hla):

	# In case it's not a string, return as is
	if type(hla) != str:
		return hla

	# Remove non-alphanumeric
	hla = hla.replace('HLA-', '')
	hla = re.sub(r'[^A-Za-z0-9]', '', hla)

	# Splits heterodymer chains
	chain = re.findall(r'[A-Za-z]+[0-9]+', hla)

	alleles = list()

	for chain in chain:
		locus = chain[:-4]
		allele_fields = chain[-4:]
		formatted_fields = allele_fields[:2] + ':' + allele_fields[2:]
		formatted_allele = locus + '*' +  formatted_fields
		alleles.append(formatted_allele)
	
	return '/'.join(alleles)


# ------------------------------------------------------

if __name__ == '__main__':
	hlas = ['DPA101:03/DPB103:01', 'HLA-A02:01']

	for hla in hlas:
		print(standardize_name(hla))

