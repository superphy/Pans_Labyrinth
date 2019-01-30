from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import hashlib

# Function to create kmers from a contig sequence
# Input: fasta file
# Input: hash string of the fasta file
# Output: series with the sequences and contig number
def create_kmer(infile, hash):
	sequence_data = []
	contig_data = []
	for seq_record in SeqIO.parse(infile, "fasta"):
		i = 0
		j = 11
		contig = seq_record.id
		my_seq_orig = seq_record.seq
		for x in my_seq_orig:
			if j == len(my_seq_orig) + 1:
				break
			my_seq = my_seq_orig[i:j]
			sequence_data.append(str(my_seq))
			contig_data.append(contig)
			j += 1
			i += 1
	df = pd.Series(sequence_data, index = contig_data)
	df.to_pickle('data/kmers/' + hash +'_kmers.pkl')

# Virtually the same function as above but creates kmers based off of the reverse compliment of the contig
def create_RC_kmer(infile, hash):
	sequence_data = []
	contig_data = []
	for seq_record in SeqIO.parse(infile, "fasta"):
		i = 0
		j = 11
		contig = seq_record.id
		my_seq_orig = seq_record.seq.reverse_complement()
		for x in my_seq_orig:
			if j == len(my_seq_orig) + 1:
				break
			my_seq = my_seq_orig[i:j]
			sequence_data.append(str(my_seq))
			contig_data.append(contig)
			j += 1
			i += 1
	df = pd.Series(sequence_data, index = contig_data)
	df.to_pickle('data/RC_kmers/' + hash +'_rc_kmers.pkl')

# Function to split a string and return the first index of the split
# Modify the index number to get a different index
# Input: string to be split and the character to split on
# Output: split string at index 0
def split_string(string, character):
	string = string.split(character)
	return string[0]
