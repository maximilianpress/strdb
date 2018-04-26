#!/usr/bin/python
#


import sys
import numpy
from collections import Counter, defaultdict

def parse_tsv_to_table(tsv_path, ok_string=""):
    '''Parse the tsv to a dict of dicts. Not using csv.DictReader because it's garbage.

    Args:
        tsv_path (str): path to the tsv file

    Returns:
        table ({str : {str : str}): the table as a dict of dicts.

    '''
    table = {}
    file = open(tsv_path, 'r')
    header = False
    for line in file:
        fields = line.strip().split("\t")
        # check in case there is predictable garbage in a tsv file.
        if not header:
            header = fields
            continue
        if is_bad_report_line(fields, ok_string=ok_string):
            continue
        rowname = fields[0]
        # this expression handles the possibility of len(header_ != len(fields)
        table[rowname] = {header[col]: fields[col] for col in range(min(len(header), len(fields)))}
    file.close()
    return header, table
    
def read_and_clean_data(annot_path, geno_path):
	'''Read in genotype and annotations data for STRs to use.'''
	# commented out a lot of preprocessing that was better handled by 
	# just altering the source data

	annot_header, annots = parse_tsv_to_table(annot_path)
	#strain_header, strains = parse_tsv_to_table(libinfo_path)
	geno_header, genos = parse_tsv_to_table(geno_path)
	#prelim_header, genos = parse_tsv_to_table(geno_path)
	# dict comprehension mapping strain names to stupid IDs
	# strain_maps = {strains[new_name]["spikein_file"] : new_name 
# 		for new_name in strains.keys() if new_name is not None}   
# 	geno_header = [strain_maps[old] for old in header]
# 	# transform values to floats, correct column names, handle missing data
# 	for row in genos.keys():
		for col in geno_header:
			#new_col_name = strain_maps[col]
			#genos[row][new_col_name] = float(genos[row][col]) if 
			genos[row][col] = float(genos[row][col]) if
				genos[row][col] is not "NA" else None
			#del col in genos[row]
	return annots, annot_header, geno, geno_header
	

def get_high_conf_markers(geno_table, strains=None, threshold = 0.2):
	'''Just extract genotypes for markers of high quality, defined as having no more than 
	some amount of missing data.'''
	rowcounts = {}
	good_markers = []
	for row in geno_table.keys():
		row_values = geno_table[row].values() if strains is None else
			[geno_table[row][col] for col in strains]
		counted = Counter()
		rowcounts[row] = counted
		if (None is not in counted) or ((counted[None] / len(counted)) < threshold):
			good_markers.append(row)
	return good_markers
	
	
def compare_two_strains(geno_path, annot_path, strain1, strain2, num_markers = None):
	'''generate a bunch of data regarding STR differences between 2 strains'''
	
	annots, annot_header, genos, geno_header = read_and_clean_data(annot_path, geno_path)
	if any([strain1, strain2] not in geno_header):
		ValueError(
		"Data not available for at last one of {0}. Available data for: {1}".format(
		" ".join([strain1,strain2]), " ".join(geno_header)
		) 
	good_markers = get_high_conf_markers(genos, strains = [strain1, strain2])
	