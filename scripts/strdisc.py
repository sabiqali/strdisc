import argparse
from itertools import combinations, count
import pysam
import parasail

def find_all_substrings(region):
	# Extract K length substrings
	# Using itertools.combinations()
	res = [region[x:y] for x, y in combinations(range(len(region) + 1), r = 2) if len(region[x:y]) >= K ]
	return res

def KMPSearch(pat, txt):
	M = len(pat)
	N = len(txt)
	index = list()
	# create lps[] that will hold the longest prefix suffix 
	# values for pattern
	lps = [0]*M
	j = 0 # index for pat[]
	
	# Preprocess the pattern (calculate lps[] array)
	computeLPSArray(pat, M, lps)
	
	i = 0 # index for txt[]
	while i < N:
		if pat[j] == txt[i]:
			i += 1
			j += 1
  
		if j == M:
			#print ("Found pattern at index " + str(i-j))
			index.append(i-j)
			j = lps[j-1]
  
		# mismatch after j matches
		elif i < N and pat[j] != txt[i]:
			# Do not match lps[0..lps[j-1]] characters,
			# they will match anyway
			if j != 0:
				j = lps[j-1]
			else:
				i += 1
	return index

  
def computeLPSArray(pat, M, lps):
	len = 0 # length of the previous longest prefix suffix

	lps[0] # lps[0] is always 0
	i = 1

	# the loop calculates lps[i] for i = 1 to M-1
	while i < M:
		if pat[i]== pat[len]:
			len += 1
			lps[i] = len
			i += 1
		else:
			# This is tricky. Consider the example.
			# AAACAAAA and i = 7. The idea is similar 
			# to search step.
			if len != 0:
				len = lps[len-1]

			# Also, note that we do not increment i here
			else:
				lps[i] = 0
				i += 1

parser = argparse.ArgumentParser()
parser.add_argument('--length', help='the minimum length of STRs to check', required=False, default=3)
parser.add_argument('--indel_file', help='the file generated from the SV detection pipeline', required=False)
parser.add_argument('--count', help='the minimum number of repeets to consider the region viable', required=False, default=3)
parser.add_argument('--bed', help='the output bed file for the coordinates of the potential STRs', required=False)
parser.add_argument('--bam', help='the bam file', required=False)

args = parser.parse_args()

length_repeat = args.length
indel_file = args.indel_file
no_of_repeats = args.count
output_file = args.bed
in_bam = args.bam

bamfile = pysam.AlignmentFile(in_bam)

# initializing K 
K = int(length_repeat)

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

	if int(support) >= 2: #if 2 out of the 3 callers call the insertion 
		region_start = int(start) - 500
		region_end = int(end) + 500

		#get the alignment of the reads and get only those reads that match the called indel region denoted by region_start and region_end
		for alignment in bamfile.fetch(chromosome,region_start,region_end):
			pair_out = alignment.get_aligned_pairs(True)
			#get the coordinates on the read that coincide with the coordinates on the reference
			for tmp_pairs in  pair_out:
				if abs(tmp_pairs[1] - region_start) <= 10:
					region_read_start = tmp_pairs[1]
				if abs(tmp_pairs[1] - region_end) <= 10:
					region_read_end = tmp_pairs[1]
			
			#get the sequence from the alignment and then get the subsequence
			test_seq = alignment.query_sequence[region_read_start:region_read_end]

			#get all the substrings of the subsequence
			all_substrings = find_all_substrings(test_text)

			#go through all the substrings of the subsequence and then check if that is a repeated substring or not, if it is document it in the repeat_of_interest dictonary(can be done using KMP for exact matches and parasail for approximate matches)
			for x in range(len(all_substrings)):
				counter = 1
				indices = KMPSearch(all_substrings[x],test_text)  #generate a list of indeces where the substring has been found
				diffs = [j-i for i, j in zip(indices[:-1], indices[1:])]  #generate a list of differences between the elements of the indeces list. 

				for diff in diffs:  #if the difference between the elements of the indeces list is equal to the length of the substring, it is a repeat, we count it.
					if diff == len(all_substrings[x]):
						counter = counter + 1
				if counter >= max_repeat_count and counter >= int(no_of_repeats):  #if local maxima is greater than global maxima and greater than the threshold, we write it to the repeat of interest. 
					if all_substrings[x] in repeat_of_interest:
						if repeat_of_interest[all_substrings[x]] < counter:
							repeat_of_interest.update({all_substrings[x]: counter})
					else:
						repeat_of_interest[all_substrings[x]] = counter
					max_repeat_count = counter  #we store the global maxima
					max_repeat_substring = all_substrings[x]

	for key in repeat_of_interest:  #we print only the elements of the dictonary which are equal to the global maxima
		if(repeat_of_interest[key] == counter):
			print(key, repeat_of_interest[key])
