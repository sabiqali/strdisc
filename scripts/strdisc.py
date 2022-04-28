import argparse
from itertools import combinations, count
from turtle import st
import pysam
import edlib
import math
import itertools

def find_all_substrings(region):
	# Extract K length substrings
	# Using itertools.combinations()
	res = [region[x:y] for x, y in combinations(range(len(region) + 1), r = 2) if len(region[x:y]) >= K and len(region[x:y]) <= int(upper_length) ]
	return res

def hamming_distance(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Strand lengths are not equal!")
    else:
        return sum(1 for (a, b) in zip(str1, str2) if a != b)

def allCharactersSame(s) :
    n = len(s)
    for i in range(1, n) :
        if s[i] != s[0] :
            return False
 
    return True

def reverse_complement(seq):  
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	rev = "".join(complement.get(base, base) for base in reversed(seq))
	return rev

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

def calc_ed(pat1, pat2):
	result = edlib.align(pat1, pat2, mode = "HW", task = "path")
	return result['editDistance']

parser = argparse.ArgumentParser()
parser.add_argument('--lower_length', help='the minimum length of STRs to check', required=False, default=3)
parser.add_argument('--upper_length', help='the minimum length of STRs to check', required=False, default=6)
parser.add_argument('--indel_file', help='the file generated from the SV detection pipeline', required=False)
parser.add_argument('--count', help='the minimum number of repeets to consider the region viable', required=False, default=3)
parser.add_argument('--bed', help='the output bed file for the coordinates of the potential STRs', required=False)
parser.add_argument('--bam', help='the bam file', required=False)

args = parser.parse_args()

lower_length = args.lower_length
upper_length = args.upper_length
indel_file = args.indel_file
no_of_repeats = args.count
output_file = args.bed
in_bam = args.bam

bamfile = pysam.AlignmentFile(in_bam)

if output_file:
	output_file_fd = open(output_file, "a+")

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
	repeat_of_interest_reverse = dict()  #variable to hold the highest repeating substring, which could potentially be the repeating unit of the STR.
	max_repeat_substring_reverse = ""
	max_repeat_count_reverse = 0

	if int(support) >= 2: #if 2 out of the 3 callers call the insertion 
		region_start = int(start) - 500
		region_end = int(end) + 500

		#get the alignment of the reads and get only those reads that match the called indel region denoted by region_start and region_end
		for alignment in bamfile.fetch(chromosome,region_start,region_end):
			pair_out = alignment.get_aligned_pairs(True)
			strand = '-' if alignment.is_reverse else '+'
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

			#get all the substrings of the subsequence
			#print(test_seq)
			all_substrings = find_all_substrings(test_seq)
			if output_file:
				print("processing read", read_count)
			
			#Using Edlib to compare the strings and calculate their edit distance to get the similarity between the substrings. 
			#if strand == '+':
			#	for i in range(len(all_substrings)):
			#		counter = 0
			#		for j in range(len(all_substrings)):
			#			if len(all_substrings[i]) == len(all_substrings[j]):
			#				ed = calc_ed(all_substrings[i],all_substrings[j])
			#				if ed/len(all_substrings[i]) <= 0.10:
			#					counter = counter + 1
			#				if counter >= max_repeat_count:  #if local maxima is greater than global maxima and greater than the threshold, we write it to the repeat of interest. 
			#					if all_substrings[i] in repeat_of_interest:
			#						if repeat_of_interest[all_substrings[i]] < counter:
			#							repeat_of_interest.update({all_substrings[i]: counter})
			#					else:
			#						repeat_of_interest[all_substrings[i]] = counter
			#					max_repeat_count = counter  #we store the global maxima
			#					max_repeat_substring = all_substrings[i]
			#					#print(max_repeat_substring, max_repeat_count)
			#			else:
			#				continue
			#else:
			#	for i in range(len(all_substrings)):
			#		counter = 0
			#		for j in range(len(all_substrings)):
			#			if len(all_substrings[i]) == len(all_substrings[j]):
			#				ed = calc_ed(all_substrings[i],all_substrings[j])
			#				if ed/len(all_substrings[i]) <= 0.10:
			#					counter = counter + 1
			#				if counter >= max_repeat_count:  #if local maxima is greater than global maxima and greater than the threshold, we write it to the repeat of interest. 
			#					if all_substrings[i] in repeat_of_interest_reverse:
			#						if repeat_of_interest_reverse[all_substrings[i]] < counter:
			#							repeat_of_interest_reverse.update({all_substrings[i]: counter})
			#					else:
			#						repeat_of_interest_reverse[all_substrings[i]] = counter
			#					max_repeat_count_reverse = counter  #we store the global maxima
			#					max_repeat_substring_reverse = all_substrings[i]
			#					#print(max_repeat_substring, max_repeat_count)
			#			else:
			#				continue

			#Using Hamming Distance to figure out the similarity between the substrings
			if strand == '+':
				for i in range(len(all_substrings)):
					counter = 0
					for j in range(len(all_substrings)):
						if len(all_substrings[i]) == len(all_substrings[j]):
							hd = hamming_distance(all_substrings[i],all_substrings[j])
							if hd <= 1:
								counter = counter + 1
							if counter >= max_repeat_count:  #if local maxima is greater than global maxima and greater than the threshold, we write it to the repeat of interest. 
								if all_substrings[i] in repeat_of_interest:
									if repeat_of_interest[all_substrings[i]] < counter:
										repeat_of_interest.update({all_substrings[i]: counter})
								else:
									repeat_of_interest[all_substrings[i]] = counter
								max_repeat_count = counter  #we store the global maxima
								max_repeat_substring = all_substrings[i]
								#print(max_repeat_substring, max_repeat_count)
						else:
							continue
			else:
				for i in range(len(all_substrings)):
					counter = 0
					for j in range(len(all_substrings)):
						if len(all_substrings[i]) == len(all_substrings[j]):
							hd = hamming_distance(all_substrings[i],all_substrings[j])
							if hd <= 1:
								counter = counter + 1
							if counter >= max_repeat_count:  #if local maxima is greater than global maxima and greater than the threshold, we write it to the repeat of interest. 
								if all_substrings[i] in repeat_of_interest_reverse:
									if repeat_of_interest_reverse[all_substrings[i]] < counter:
										repeat_of_interest_reverse.update({all_substrings[i]: counter})
								else:
									repeat_of_interest_reverse[all_substrings[i]] = counter
								max_repeat_count_reverse = counter  #we store the global maxima
								max_repeat_substring_reverse = all_substrings[i]
								#print(max_repeat_substring, max_repeat_count)
						else:
							continue

			if output_file:
				print("processed read ", read_count)
			read_count = read_count + 1
			#print(repeat_of_interest)
			#print(max_repeat_count)
			#print(repeat_of_interest_reverse)
			#print(max_repeat_count_reverse)

	#for key in repeat_of_interest:  #we print only the elements of the dictonary which are equal to the global maxima
	#	if(repeat_of_interest[key] == max_repeat_count):
	#		print(key, repeat_of_interest[key])
	
	#for key in repeat_of_interest_reverse:  #we print only the elements of the dictonary which are equal to the global maxima
	#	if(repeat_of_interest_reverse[key] == max_repeat_count_reverse):
	#		print(key, repeat_of_interest_reverse[key])

	fw_max_key = sorted(repeat_of_interest, key=repeat_of_interest.get, reverse=True)[:3]
	rw_max_key = sorted(repeat_of_interest_reverse, key=repeat_of_interest_reverse.get, reverse=True)[:3]
	#fw_output = ""
	#rw_output = ""
	for key in fw_max_key:
		if not allCharactersSame(key):
			#fw_output = fw_output + key + " " + str(repeat_of_interest[key]) + " "
			if output_file:
				print("\t".join([chromosome, start, end, key]), file=output_file_fd)
				print("Finished processing")
			else:
				print("\t".join([chromosome, start, end, key]))
	for key in rw_max_key:
		if not allCharactersSame(key):
			#rw_output = rw_output + key + " " + str(repeat_of_interest_reverse[key]) + " "
			if output_file:
				print("\t".join([chromosome, start, end, reverse_complement(key)]), file=output_file_fd)
				print("Finished processing")
			else:
				print("\t".join([chromosome, start, end, reverse_complement(key)]))
	#print("Forward strand max:")
	#print(fw_output)
	#print("Reverse strand max:")
	#print(rw_output)
	#print("%s\t%s\t%s\t%s"%(chromosome, start, end, ))
