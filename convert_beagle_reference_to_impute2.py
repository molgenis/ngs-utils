
def error_unknown_genotype(lc, gen, al1, al2):
	raise Exception("Unknown genotype: " + gen + " in line: " + str(lc) + ". Accepter alleles according to markers file: " + al1 + " , " + al2)

def Convert_beagle_to_impute2_reference_user_Kantale(
	input_beagle_filename = None,
	input_markers_filename = None,
	output_hap_filename = None,
	output_legend_filename = None,
	):

	input_beagle_file = open(input_beagle_filename)
	input_markers_file = open(input_markers_filename)

	output_hap_file = open(output_legend_filename, 'w')
	output_legend_file = open(output_legend_filename, 'w')
	#Write header tp legend file
	output_legend_file.write('ID pos allele0 allele1')
	lc = 0

	for beagle_line in input_beagle_file:
		lc += 1
		beagle_line_s = beagle_line.replace('\n', '').split()

		#Skip headers in beagle file
		if beagle_line_s[0] != 'M':
			continue

		#Reader markers line
		markers_line = input_markers_file.readline()
		markers_line_s = markers_line.replace('\n', '').split()

		rs_id = markers_line_s[0]
		pos = markers_line_s[1]
		allele0 = markers_line_s[2]
		allele1 = markers_line_s[3]

		beagle_rs_id = beagle_line_s[1]
		if beagle_rs_id != rs_id:
			raise Exception("Beagle and markers lines are not aligned. Markers differ in line:" + str(lc))

		legend_to_print = [rs_id, pos, allele0, allele1]

		genotypes = beagle_line_s[2:]
		hap_to_print = [0 if genotype == allele0 else 
							1 if genotype == allele1 else 
								error_unknown_genotype(lc, genotype, allele0, allele1) 
									for genotype in genotypes]

		#Print the results
		output_hap_file.write(str.join(' ', hap_to_print) + '\n')
		output_legend_file.write(str.join(' ', legend_to_print) + '\n')

	output_hap_file.close()
	output_legend_file.close()

import sys

if __name__ == '__main__':
	Convert_beagle_to_impute2_reference_user_Kantale(sys.args[1], sys.args[2], sys.args[3], sys.args[4])


