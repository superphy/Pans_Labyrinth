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

	kmers = {"genome" : ["AAAA", "TTTT", "CCCC", "GGGG"]}
	dgraph.add_kmers_dgraph(client, kmers, "test_genome")

	sg1 = dgraph.example_query(client, "test_genome")

	kmer_list = []
	for x in sg1:
		if x == sg1[-1]:
			dict = sg1[-1]
			value = list(dict.values())
			kmer = value[0]
			last_kmer = kmer[0]
			kmer_list.append(last_kmer["kmer"])
		else:
			kmer = x["kmer"]
			kmer_list.append(kmer)
	#print(kmer_list)
	for i, x in enumerate(kmer_list):
		assert x == kmers["genome"][i]

#def test_insertion(): #TODO for the command line arguments

#def test_deletion(): #TODO for the command line arguments

def test_hash_verification():
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
			dgraph.add_kmers_dgraph(client, all_kmers, genome)

	assert genome == "genome_" + commandline.compute_hash(filepath)



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

def test_verify_contig():
	"""
	A test to verify that the kmers returned from the graph are the same as the
	kmers in the original contig
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
			dgraph.add_kmers_dgraph(client, all_kmers, genome)

	#Query graph and put all kmers into a list
	sg1 = dgraph.path_query(client, genome)

	kmer_list = []
	for i, x in enumerate(sg1["path"]):
		kmer = sg1["path"][i]["kmer"]
		kmer_list.append(kmer)
	#print(kmer_list)

	#Gets the whole first kmer in the list and then the last character of the rest of the kmers in the list
	#Creates a string representing the contig
	first, *rest = kmer_list
	ends = [kmer[-1] for kmer in rest]
	contig = ''.join([first] + ends)

	#Gets the actual contig from the fasta file and turns it into  string
	with open(path + "/test.fasta", "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			sequence = record.seq
	sequence_string = str(sequence)

	assert contig == sequence_string[:-1]
