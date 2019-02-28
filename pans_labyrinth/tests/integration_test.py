import pytest
from pans_labyrinth import main, files, dgraph, commandline, logging_functions
import os
from Bio import SeqIO


def test_query():
	"""
	Function to test that a query returns the correct result. Inserts
	a known value then queries for that value, assert that the returned value is the
	same as the inserted value.
	"""
	stub = dgraph.create_client_stub()
	client = dgraph.create_client(stub)
	dgraph.drop_all(client)
	dgraph.add_schema(client)
	dgraph.add_genome_to_schema(client, "test_genome")
	kmers = ["AAAA", "TTTT", "CCCC", "GGGG"]
	length = len(kmers)
	x = 0
	while x <= length -2:
		dgraph.add_kmer_to_graph(client, kmers[x], kmers[x+1], "test_genome")
		x += 1
	sg1 = dgraph.example_query(client, "test_genome")
	value = sg1["genome"]
	klist = []
	for counter, x in enumerate(value):
		temp = value[counter]
		klist.append(temp["kmer"])
	for i, x in enumerate(klist):
		assert x == kmers[i]

#def test_insertion(): #TODO for the command line arguments

#def test_deletion(): #TODO for the command line arguments

#def hash_verification(): #TODO verify the hash of a genome

def test_no_file():
	"""
	Function to ensure that if no file is present in the given directory then
	The program fails.
	"""
	path = os.path.abspath("pans_labyrinth/tests/data/empty_test")
	with pytest.raises(SystemExit) as se:
		stub = dgraph.create_client_stub()
		client = dgraph.create_client(stub)
		dgraph.drop_all(client)
		dgraph.add_schema(client)
		dgraph.add_genome_to_schema(client, "test_genome")
		for filepath in files.walkdir(path):
			with open(filepath, 'rb') as file:
				dgraph.create_graph(client, file, filepath)
	assert se.type == SystemExit

def test_non_fasta():
	"""
	Function to test that the program only accepts .fasta files. Could
	create issues if the program was reading non-fasta files and inserting the
	values into the graph.
	"""
	path = os.path.abspath("pans_labyrinth/tests/data/bad_fasta")
	with pytest.raises(SystemExit) as se:
		stub = dgraph.create_client_stub()
		client = dgraph.create_client(stub)
		dgraph.drop_all(client)
		dgraph.add_schema(client)
		for filepath in files.walkdir(path):
			with open(filepath, 'rb') as file:
				dgraph.create_graph(client, file, filepath)
	assert se.type == SystemExit

def test_verify_contig(): # TODO bulk load kmers
	"""
	A test to verify that the kmers
	"""

	#Create graph and fill it with kmers from a fasta file
	path = os.path.abspath("pans_labyrinth/tests/data/contig_test")
	stub = dgraph.create_client_stub()
	client = dgraph.create_client(stub)
	dgraph.drop_all(client)
	dgraph.add_schema(client)
	for filepath in files.walkdir(path):
		with open(filepath, 'rb') as file:
			x = 0
			filename = file.name
			genome = "genome_" + commandline.compute_hash(filepath)
			dgraph.add_genome_to_schema(client, genome)
			all_kmers = dgraph.get_kmers_files(filename, 11)
			kmers = all_kmers['SRR1122659.fasta|NODE_1_length_767768_cov_21.1582_ID_10270']
			length = len(kmers)
			while x <= length -2:
				dgraph.add_kmer_to_graph(client, kmers[x], kmers[x+1], genome)
				x += 1

	#Query graph and put all kmers into a list
	sg1 = dgraph.example_query(client, genome)
	value = sg1["genome"]
	klist = []
	for counter, x in enumerate(value):
		temp = value[counter]
		klist.append(temp["kmer"])

	#Gets the whole first kmer in the list and then the last character of the rest of the kmers in the list
	#Creates a string representing the contig
	first, *rest = klist
	ends = [kmer[-1] for kmer in rest]
	contig = ''.join([first] + ends)

	#Gets the actual contig from the fasta file and turns it into  string
	with open(path + "/test.fasta", "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			sequence = record.seq
	sequence_string = str(sequence)

	assert contig == sequence_string[:-1]
