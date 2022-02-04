import argparse
from itertools import combinations, count
import pysam
import edlib
import math
import itertools
import os

def find_all_substrings(region):
	# Extract K length substrings
	# Using itertools.combinations()
	res = [region[x:y] for x, y in combinations(range(len(region) + 1), r = 2) if len(region[x:y]) >= K and len(region[x:y]) <= upper_length ]
	return res

parser = argparse.ArgumentParser()
parser.add_argument('--lower_length', help='the minimum length of STRs to check', required=False, default=3)
parser.add_argument('--upper_length', help='the minimum length of STRs to check', required=False, default=6)
parser.add_argument('--indel_file', help='the file generated from the SV detection pipeline', required=False)
parser.add_argument('--count', help='the minimum number of repeets to consider the region viable', required=False, default=3)
parser.add_argument('--bed', help='the output bed file for the coordinates of the potential STRs', required=False)
parser.add_argument('--bam', help='the bam file', required=False)
parser.add_argument('--ref', help='the Reference file for the reads', required=False, default='/.mounts/labs/simpsonlab/users/schaudhary/projects/2022.1.STRDiscovery/references/GRCh38_ref')

args = parser.parse_args()

lower_length = args.lower_length
upper_length = args.upper_length
indel_file = args.indel_file
no_of_repeats = args.count
output_file = args.bed
in_bam = args.bam
ref = args.ref

bamfile = pysam.AlignmentFile(in_bam)

# initializing K 
K = int(lower_length)

test_text = "bacbacbacabcabcabcabcabcabc"

indel_fh = open(indel_file)
#header = indel_fh.readline()           activate this when Mike adds a header to the indel_file

for line in indel_fh:
	chromosome,start,end,type_of,sv_length,sv_tmp,support = line.rstrip().split()[0:7]  #read the indel file and get data from the different columns. 
	region_start = start
	region_end = end
	repeat_of_interest = dict()  #variable to hold the highest repeating substring, which could potentially be the repeating unit of the STR.
	max_repeat_substring = ""
	max_repeat_count = 0
	read_count = 1

	if int(support) >= 2: #if 2 out of the 3 callers call the insertion 
		region_start = int(start) - 500
		region_end = int(end) + 500

		#get the alignment of the reads and get only those reads that match the called indel region denoted by region_start and region_end
		for alignment in bamfile.fetch(chromosome,region_start,region_end):
			pair_out = alignment.get_aligned_pairs(True)
	
			#get the coordinates on the read that coincide with the coordinates on the reference
			for tmp_pairs in  pair_out:
				if abs(tmp_pairs[1] - region_start) <= 10:
					region_read_start = tmp_pairs[0]
				if abs(tmp_pairs[1] - region_end) <= 10:
					region_read_end = tmp_pairs[0]

			#get the sequence from the alignment and then get the subsequence
			if region_read_start < region_read_end:
				test_seq = alignment.query_sequence[region_read_start:region_read_end]
			else:
				continue
			
			f_ref = open("references/temporary_ref.fa", "w+")
			print(">ref", file=f_ref)
			print(test_seq, file=f_ref)
			shell_out = os.system("cd references && makeblastdb -in temporary_ref.fa -dbtype nucl")


			if shell_out == 0:
				print("blast db successfully made")
			#get all the substrings of the subsequence
			#print(test_seq)
			all_substrings = find_all_substrings(test_seq)
			print("processing read", read_count)
			f = open("reads/temporary_reads.fa", "w+")

			seq_counter = 0
			for i in range(len(all_substrings)):
				seq_counter += 1
				print(">seq"+ str(seq_counter), file=f)
				print(all_substrings[i], file=f)

			#print("blastn -db references/temporary_ref.fa -query reads/temporary_reads.fa -word_size 4 -out align.out")

			shell_out = os.system("blastn -db " + ref + " -query " + os.getcwd() + "/temporary_reads.fa -word_size 4 -out align.out")

			if shell_out == 0:
				print("blast align successful")

			print("processed read ", read_count)
			read_count = read_count + 1
			break
		break

	break			

	#for key in repeat_of_interest:  #we print only the elements of the dictonary which are equal to the global maxima
		#if(repeat_of_interest[key] == max_repeat_count):
			#print(key, repeat_of_interest[key])