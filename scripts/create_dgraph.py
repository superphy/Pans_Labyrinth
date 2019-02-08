#!/usr/bin/env python

"""
This is a simple example, using the layout of
https://github.com/dgraph-io/pydgraph/blob/master/examples/simple/simple.py
We have modified it for genome kmers
"""

import pydgraph
import json
import tempfile
import subprocess
import timeit
from Bio import SeqIO
from multiprocessing import Pool, Process
from functools import partial
import sys
import os
import hashlib
from tqdm import tqdm
import argparse

def create_client_stub():
	"""
	This allows us to create as many stubs as we want easily
	without having to specify the address every time.
	This allows a single source to be easily changed.
	:return: A client stub
	"""

	return pydgraph.DgraphClientStub('localhost:9080')


def create_client(client_stub):
	"""
	This allows us to create multiple clients easily,
	by passing in the stub from create_client_stub()
	:param client_stub: connection for stub
	:return: dgraph client
	"""
	return pydgraph.DgraphClient(client_stub)


def drop_all(client):
	"""
	Wipe the database for a fresh start
	:param client: the dgraph client
	:return: dgraph response
	"""

	return client.alter(pydgraph.Operation(drop_all=True))


def set_schema(client):
	"""
	We can define whatever type we like for our schema.
	For this simple example we will have only kmers and edges as genome names.
	add_genome_to_schema() allows genomes to be added programmatically.
	:param client: dgraph client
	:return: The client altered via the schema set out here
	"""
	schema = """
		kmer: string @index(exact, term) .
	"""

	return client.alter(pydgraph.Operation(schema=schema))


def add_kmer_to_graph(client, ki, kn, genome, bulk=False):
	"""
	Every kmer needs to be linked to another kmer. Single kmers not permitted.
	:param client: dgraph client
	:param ki: the initial kmer. Get the uid if it exists. Otherwise create it.
	:param kn: the next kmer, linked to ki. Get the uid if it exists. Otherwise create it.
	:param genome: the genome label between the two kmers -- this needs to be previously indexed in the schema.
	:param bulk: delay the insertion of the data, instead return a list of triples to be live-loaded
	:return: None
	"""

	# First we need to see if the two kmers exist already
	# This is an upsert. If it exists, get the uid.
	# Otherwise, create the uid.
	uid_ki = kmer_upsert(client, ki)
	uid_kn = kmer_upsert(client, kn)

	# Data to be inserted
	d = "<{0}> <{1}> <{2}> .".format(uid_ki, genome, uid_kn)

	if not bulk:
		# start the transaction
		txn = client.txn()

		# The two kmers (as nodes) are linked by a <genome name> predicate
		try:
			txn.mutate(set_nquads=d)
			txn.commit()

		finally:
			txn.discard()
	else:
		return d

def kmer_upsert(client, kmer):
	"""
	Check the existence of kmer. If it exists, return the uid.
	If it doesn't exist, create it, then return the uid.
	:param client: the dgraph client
	:param kmer: kmer sequence
	:return: uid of kmer
	"""

	uid = kmer_query(client, kmer)

	if not uid:
		# kmer does not exist in database. Create it.
		txn = client.txn()

		try:
			# The data we wish to add in triple form
			d = """
					_:blank-0 <kmer> "{0}" .
			""".format(kmer)

			m = txn.mutate(set_nquads=d)
			txn.commit()
			uid = m.uids['blank-0']

		finally:
			# Good practice according the the dgraph docs
			# No network overhead associated with it
			txn.discard()

	return uid


def kmer_query(client, kmer):
	"""
	Strictly queries the database for a kmer.
	Returns the uid if present, None if not.
	:param client: the dgraph client
	:param kmer: the kmer we would like to check
	:return: uid or None
	"""

	# Transaction for the query
	query = """
		query find_kmer($k :string){
			find_kmer(func: eq(kmer, $k)){uid}
		}
	"""
	variables = {'$k': kmer}
	res = client.query(query, variables=variables)
	json_res = json.loads(res.json)

	if json_res['find_kmer']:
	  return json_res['find_kmer'][0]['uid']
	else:
		return None


def kmer_multiple_query(client, kmer_list):
	"""
	Bulk query a list of kmers and return a dictionary of kmer:uid.
	:param client: dgraph client
	:param kmer_list: List of kmers to query dgraph for
	:return: [dict{kmer:uid}]
	"""
	query = """
		query find_all($klist: string){
			find_all(func: anyofterms(kmer, $klist))
			{
				uid
				kmer
			}
	}
	"""
	variables = {'$klist': ' '.join(kmer_list)}
	res = client.query(query, variables=variables)
	json_res = json.loads(res.json)

	if json_res['find_all']:
		return json_res['find_all']
	else:
		return None



def example_query(client, genome):
	"""
	Example of getting a subgraph from a predicate query
	:param client: dgraph client
	:param genome: genome of interest as string (needs to be indexed in schema)
	:return: sub-graph as json
	"""

	# Transaction for the query
	query = """
		{{
			genome(func: has({0})){{
				uid
				kmer
			}}
		}}
	""".format(genome)

	# This gets all but the last uid in the format {genome: [{'uid':'0x335'}, {'uid':'0x336'}]}
	res = client.query(query)
	j_res = json.loads(res.json)

	# In the graph A-p->B-p->C has(p) will return the uid for A and B, but not C
	# We want the last uid as well, so will construct a second query for it, using the last uid
	last_uid = j_res['genome'][-1]['uid']

	# Transaction for last query
	last_query = """
		{{
		  lq(func: uid({0})){{
			{1}{{
				uid
				kmer
			}}
		  }}
		}}
	""".format(last_uid, genome)

	# This gets the last query
	l_res = client.query(last_query)
	j_l_res = json.loads(l_res.json)
	return j_res['genome'] + j_l_res['lq']
	#j_l_res['lq'][genome])


def add_genome_to_schema(client, genome):
	"""
	Index the genome name as a predicate, so functions can be used on it when searching etc.
	If the genome name is not added to the schema, it will not be indexed, and functions can
	not be applied to it. We will use the upsert directive when modifying the schema, so that
	if it exists, nothing happens, but new names are added to the schema.
	:param client: dgraph client
	:param genome: genomeName
	:return: Altered dgraph
	"""
	schema = """
		{0}: uid .
	""".format(genome)

	return client.alter(pydgraph.Operation(schema=schema))


def run_subprocess(cmd):
	"""
	Run a subprocess from python
	:param cmd: the command to run
	:return: The completed process
	"""

	start_time = timeit.default_timer()
	comp_proc = subprocess.run(
		cmd,
		shell=False,
		check=False,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE
	)

	if comp_proc.returncode == 0:
		elapsed_time = timeit.default_timer() - start_time
		print("Subprocess {} finished successfully in {:0.3f} sec.".format(cmd, elapsed_time))
		return comp_proc
	else:
		print("Error in subprocess. The following command failed: {}".format(comp_proc))
		exit("subprocess failure")


def kmer_from_file(filename, kmer_size):
	"""
	 Read in a fasta file. Return a dict of lists all kmers of given size.
	:param filename: Fasta file to process
	:param kmer_size: Size of kmer
	:return: Dict of lists of all kmers dict{contig:[kmers]}
	"""
	all_kmers = {}
	with open(filename, "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			all_kmers[record.id]=[]
			for i in range(0, len(record.seq) - kmer_size - 2, kmer_size):
				all_kmers[record.id].append(str(record.seq[0+i:kmer_size+i]))
	return all_kmers


def kmer_process_file(client, filename, kmer_size):
	"""

	Using BioPython, and lots of RAM.
	:param client: dgraph client
	:param filename: Fasta file to process
	:param kmer_size: Size of kmer to process file into
	:return: A list of all kmers
	"""

	memoize_kmers = {}
	bulk_data = ""

	tfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
	with open(filename, "r") as f:
		record_count=0
		for record in SeqIO.parse(f, "fasta"):
			record_count +=1
			kmer_count = 0
			for i in range(0, len(record.seq) - kmer_size - 1, kmer_size):
				kmer_count +=1

				ks = {
					'ki': record.seq[0+i:kmer_size+i],
					'kn': record.seq[0+i+1:kmer_size+i+1]
				}

				for k in ks.keys():
					kuid = None
					if ks[k] in memoize_kmers:
						kuid = memoize_kmers[ks[k]]
					else:
						kq = kmer_query(client, ks[k])
						if kq:
							kuid = kq
							memoize_kmers[ks[k]] = kuid
						else:
							kuid = "_:b" + ks[k]
							memoize_kmers[ks[k]] = kuid
							bulk_data +=('{0} <kmer> "{1}" .{2}'.format(kuid, ks[k], "\n"))

				bulk_data +=('{0} <genomeA> {1} .{2}'.format(memoize_kmers[ks['ki']], memoize_kmers[ks['kn']], "\n"))

				if kmer_count % 1000 == 0:
					print(".", end='')
			if record_count == 1:
				break
	tfile.write(bulk_data)
	tfile.close()
	# Use the live load feature to load all the nquads
	command = [
		"dgraph", "live",
		"-r", tfile.name
	]
	bulk_return = run_subprocess(command)
	print(bulk_return)



def add_kmers_to_dict(kmer_dict, kmers):
	"""
	Updates a dictionary of kmer:uid.
	Requires the list to be in the form of [{kmer:uid}], as is returned from kmer_multiple_query()
	:param kmer_dict: {kmer:uid}
	:param kmers: [{kmer:uid}]
	:return: The updated dictionary
	"""
	if kmers:
		for ku in kmers:
			kmer_dict[ku['kmer']]=ku['uid']

	return kmer_dict


def add_all_kmers_to_graph(client, all_kmers, genome):
	"""
	Add all kmers from a given genome to the graph
	:param client: dgraph client
	:param all_kmers: dict of lists of all kmers in genome dict[contig:[kmers]]
	:param genome: name of genome to add
	:return: None
	"""

	# query a batch of kmers in bulk and get the uids if they exist
	pc = partial(process_contig_kmer, client=client, genome=genome)

	with Pool(processes=1) as pool:
		results = pool.map_async(pc, all_kmers.values())
		results.wait()



def process_contig_kmer(ckmers, client, genome):
	"""
	Process a single contig into kmers, adding nodes and edges for each
	:param ckmers: The list of kmers for the contig
	:param client: dgraph client
	:param genome: indexed edge genome name
	:return: success
	"""
	# Query for all existing kmers
	kmer_uid_dict = {}
	kmer_uid_dict = add_kmers_to_dict(kmer_uid_dict, kmer_multiple_query(client, ckmers))

	# Create list of kmers that need to be batch inserted into graph
	kmers_to_insert = []
	for kmer in ckmers:
		if kmer not in kmer_uid_dict:
			kmers_to_insert.append(kmer)

	if kmers_to_insert:
		# Bulk insert the kmers
		txn_result_dict = add_batch_kmers(client, kmers_to_insert)

		# Update the dict of kmer:uid
		kmer_uid_dict = add_kmers_to_dict(kmer_uid_dict, txn_result_dict)

	# Batch the connections between the kmers
	add_edges_to_kmers(client, ckmers, kmer_uid_dict, genome)
	return None


def add_edges_to_kmers(client, kmers, kmer_uid_dict, genome):
	"""
	Given a list of previously inserted kmers, and the corresponding dictionary of the uids
	create edges between all kmers, sequentially.
	:param client: the dgraph client
	:param kmers: list of linked kmers
	:param kmer_uid_dict: {kmer:uid}
	:param genome: the indexed edge name to connect the kmer nodes
	:return: None
	"""

	bulk_quads = []
	# We are adding sequential kmers, so we need only start at index len(kmer)-2 to
	# ensure that all kmers are linked
	for i in range(0, len(kmers)-2):
		bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(kmer_uid_dict[kmers[i]],
														 genome,
														 kmer_uid_dict[kmers[i+1]],
														 "\n"
														))
	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()


def add_batch_kmers(client, kmer_list):
	"""
	Add all the kmers in the list to the graph.
	Return the results from the transaction in the requires list of dict format.
	:param client: dgraph client
	:param kmer_list: list of kmers that need to be added to new nodes
	:return: List of {'kmer':kmer, 'uid':uid}
	"""
	# Data to be inserted
	bulk_quads = []
	for kmer in kmer_list:
		bulk_quads.append('_:{0} <kmer> "{0}" .{1}'.format(kmer, "\n"))

	# start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

		# Create the output in the form that the program is expecting
		kmer_dict_list = []
		for uid in m.uids:
			kmer_dict_list.append({'kmer':uid, 'uid':m.uids[uid]})
	finally:
		txn.discard()

	return kmer_dict_list

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

def fill_graph_progess(client):
	'''
	Fills the graph with kmers and edges based on a given folder containing genomes
	Run time is arount 1 hour
	:param client: the dgraph client
	'''
	path = "data/genomes/clean"
	filecounter = 0
	for filepath in walkdir(path):
		filecounter += 1
	for filepath in tqdm(walkdir(path), total=filecounter, unit="files"):
		with open(filepath, 'rb') as file:
			filename = file.name
			genome = "genome_" + compute_hash(filename)
			add_genome_to_schema(client, genome)
			all_kmers = kmer_from_file(filename, 11)
			add_all_kmers_to_graph(client, all_kmers, genome)


def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--insert", action = 'append', help = "Insert a new genome into the graph using a fasta file",)
	parser.add_argument("-q", "--query", action = 'append', help = "Find genome path in the graph based on the fasta file hash")
	parser.add_argument("-d", "--delete", action = 'append', help = "Remove a genome grom the graph by using a the fasta file hash")


def main():
	"""
	The program
	:return: success
	"""
	arg_parser()
	stub = create_client_stub()
	client = create_client(stub)
	drop_all(client)
	set_schema(client)
	fill_graph_progess(client)

	# manual addition to graph -- this would normally be functions
	#add_genome_to_schema(client, "genomeA")
	#add_genome_to_schema(client, "genomeB")

	#all_kmers = kmer_from_file("/tmp/pdog/pdog/data/test.fasta", 11)
	#add_all_kmers_to_graph(client, all_kmers, "genomeA")

	#query by predicate, to see the links
	#sg1 = example_query(client, "genomeA")
	#print(sg1)

	print("All done")


if __name__ == '__main__':
	main()
