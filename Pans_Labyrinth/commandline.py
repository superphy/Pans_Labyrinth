#!/user/bin/env python
from files import *
from dgraph import *
from commandline import *

"""
Arg-parser goes here
"""

def compute_hash(filename):
	'''
	Takes a path to a fasts file and creates a hash of the file to be
	used as the genome edge name
	:param filename: path tot eh fasta file to be hashed
	'''
	try:
		BUF_SIZE = 65536
		sha1 = hashlib.sha1()
		with open(filename, 'rb') as f:
			while True:
				data = f.read(BUF_SIZE)
				if not data:
					break
				sha1.update(data)
		hash = sha1.hexdigest()
	except:
		print("Failed to create hash")
		sys.exit()
	return hash


def walkdir(folder):
	'''
	Walk through each files in a directory and yeild all the paths in the dir
	:param folder: the folder containing all the files to walk through
	'''
	for dirpath, dirs, files in os.walk(folder):
		for filename in files:
			yield os.path.abspath(os.path.join(dirpath, filename))


def arg_parser(client):
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--insert", action = 'append', help = "Insert a new genome into the graph using a fasta file",)
	parser.add_argument("-q", "--query", action = 'append', help = "Find genome path in the graph based on the fasta file hash")
	parser.add_argument("-d", "--delete", action = 'append', help = "Remove a genome grom the graph by using a the fasta file hash")

	opt = parser.parse_args()
	print(opt)
	if opt.insert:
		insert_genome(client, opt.insert)
	if opt.query:
		query_for_genome(client, opt.query)
	if opt.delete:
		delete_genome(client, opt.delete)

def insert_genome(client, genomes):
	for genome in genomes:
		filename = "data/genomes/insert/{}".format(genome)
		genome = "genome_" + compute_hash(filename)
		add_genome_to_schema(client, genome)
		all_kmers = kmer_from_file(filename, 11)
		add_all_kmers_to_graph(client, all_kmers, genome)
	print("inserted genome(s)")


def query_for_genome(client, genomes):
	for genome in genomes:
		filename = "data/genomes/insert/{}".format(genome) # TODO change pathing and figure out metadata querying
		genome = "genome_" + compute_hash(filename)
		sg1 = example_query(client, genome)
		print(sg1)


def delete_genome(client, genomes):
	for genome in genomes:
		genome = "genome_" + compute_hash(filename)
