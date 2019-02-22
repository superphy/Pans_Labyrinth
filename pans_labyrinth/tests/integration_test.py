import pytest
from pans_labyrinth import main, files, dgraph, commandline, loggingFunctions
import os


def test_query():
	stub = dgraph.create_client_stub()
	client = dgraph.create_client(stub)
	dgraph.drop_all(client)
	dgraph.add_schema(client)
	dgraph.add_genome_to_schema(client, "test_genome")
	dgraph.add_kmer_to_graph(client, "AAAA", "TTTT", "test_genome")
	sg1 = dgraph.example_query(client, "test_genome")
	assert str(sg1) == "{'genome': [{'kmer': 'AAAA'}]}"

#def test_insertion(): #TODO for the command line arguments

#def test_deletion(): #TODO for the command line arguments


def test_no_file():
	path = os.path.abspath("tests/data/empty_test")
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
	path = os.path.abspath("tests/data/bad_fasta")
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
