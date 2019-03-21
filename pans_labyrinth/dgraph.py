#!/usr/bin/env python

"""
This file contains functions for adding to and querying dgraph as well as
creating the nodes and edges to be added.
Naming conventions:

verb-noun-object

Two main verbs for interacting with dgraph:
	1) add
	2) query

One additional verb for kmer construction:
	3) get
"""

import pydgraph
import json
from Bio import SeqIO
from functools import partial
from concurrent.futures import ThreadPoolExecutor
from pans_labyrinth import files, dgraph, commandline, logging_functions
import sys
import logging
import os
from collections import Counter

output_directory = os.path.abspath("data/logger")
LOG = logging_functions.create_logger()
fh = logging.FileHandler(os.path.join(output_directory, 'pans_labyrinth.log'), 'w', 'utf-8')


def create_client_stub():
	"""
	This allows us to create as many stubs as we want easily
	without having to specify the address every time.
	This allows a single source to be easily changed.
	:return: A client stub
	"""
	try:
		LOG.info("Creating client stub")
		return pydgraph.DgraphClientStub('localhost:9080')
	except:
		LOG.critical("Failed to create the client stub")


def create_client(client_stub):
	"""
	This allows us to create multiple clients easily,
	by passing in the stub from create_client_stub()
	:param client_stub: connection for stub
	:return: dgraph client
	"""
	try:
		LOG.info("Creating client")
		return pydgraph.DgraphClient(client_stub)
	except:
		LOG.critical("Failed to create client")


def drop_all(client):
	"""
	Wipe the database for a fresh start
	:param client: the dgraph client
	:return: dgraph response
	"""
	try:
		LOG.info("Dropping existing graph")
		return client.alter(pydgraph.Operation(drop_all=True))
	except:
		LOG.critical("Failed to drop previous graph")

def add_schema(client):
	"""
	We can define whatever type we like for our schema.
	For our database we will have kmers, and <genome name> as edges.
	add_genome_schema() allows genomes to be added programmatically.
	:param client: dgraph client
	:return: The client altered via the schema set out here
	"""
	schema = """
	kmer: string @index(exact, term) .
	"""
	try:
		LOG.info("Add schema to graph with client")
		return client.alter(pydgraph.Operation(schema=schema))
	except:
		LOG.critical("Failed to add schema to graph")


def query_kmers_dgraph(client, kmer_list):
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
	#print(res)
	j_res = json.loads(res.json)
	#print(j_res)


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

def path_query(client, genome):
	"""
	A query to return the proper path of a genome through the graph using the
	genome name to traverse the edges.
	:param client: the dgraph client
	:param genome: the genome in name in the form genome_hash
	:return: A dictionary converted from a JSON object which is the query returned by dgraph
	"""

	query = """
	{{
	  genome(func: has({0})){{
		uid
		kmer
		{0}{{
		  uid
		  kmer
		}}
	  }}
	}}
	""".format(genome)

	res = client.query(query)
	p_res = json.loads(res.json)

	path_list = p_res["genome"]
	uid_dict = {}
	uid_list = []
	start_stop_list = []

	for x, i in enumerate(path_list):
		first_uid = path_list[x]["uid"]
		second_uid = path_list[x][genome][0]["uid"]
		uid_list.append(first_uid)
		uid_list.append(second_uid)

	for key in uid_list:
		if key in uid_dict:
			uid_dict[key] += 1
		else:
			uid_dict[key] = 1
	#print(uid_dict)

	for key in uid_dict.keys():
		if uid_dict[key] == 1:
			start_stop_list.append(key)
	#print(start_stop_list)

	first = start_stop_list[0]
	second = start_stop_list[1]

	for x, i in enumerate(path_list):
		first_uid = path_list[x]["uid"]
		if first_uid == first:
			start = first
			stop = second
		if first_uid == second:
			start = second
			stop = first

	query1 = """
	{{
	 path as shortest(from: {0}, to: {1}){{
	   {2}
	 }}
	   path(func: uid(path)){{
	     kmer
	   }}
	}}
	""".format(start, stop, genome)

	res1 = client.query(query1)
	path_res1 = json.loads(res1.json)

	kmer_list1 = []
	for i, x in enumerate(path_res1["path"]):
		kmer = path_res1["path"][i]["kmer"]
		kmer_list1.append(kmer)

	return path_res1


def metadata_query(client, genome, uid):
	"""
	This function queries a uid which is known to have metadata associated with it
	:param client: the dgrpah client
	:param genome: the genome on the uid which you want the metadata for. Will contain
	information about the next and previous uids in the path
	:param uid: the uid of the node in which to look at
	"""
	query = """
	{{
	  metadata(func: uid({0})){{
	    {1}{{
	      expand(_all_)
	    }}
	  }}
	}}
	""".format(uid, genome)

	res = client.query(query)
	p_res = json.loads(res.json)

	print(p_res)


def add_metadata_to_schema(client, genome):
	"""
	This function adds the metadata for a certain genome to the schema. As well it adds
	the string telling identifying that a node is a duplicate and contains metadata
	:param client: the dgraph client
	:param genome: the name of the genome which contains a duplicate kmer
	"""
	schema = """
	{0}: uid .
	""".format(genome)

	client.alter(pydgraph.Operation(schema=schema))

	schema = """
	duplicate: string @index(exact, term) .
	"""

	client.alter(pydgraph.Operation(schema=schema))

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


def get_kmers_files(filename, kmer_size):
	"""
	 Read in a fasta file. Return a dict of lists all kmers of given size.
	:param filename: Fasta file to process
	:param kmer_size: Size of kmer
	:return: Dict of lists of all kmers dict{contig:[kmers]}
	"""
	LOG = logging_functions.create_logger()
	if not filename.endswith(".fasta"):
		LOG.critical("Non fasta file detected")
		sys.exit()

	all_kmers = {}
	with open(filename, "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			all_kmers[record.id] = []
			for i in range(0, len(record.seq) - kmer_size + 1, 1):
				all_kmers[record.id].append(str(record.seq[0 + i:kmer_size + i]))
	return all_kmers


def add_kmers_dict(kmer_dict, kmers):
	"""
	Updates a dictionary of kmer:uid.
	Requires the list to be in the form of [{kmer:uid}], as is returned from kmer_multiple_query()
	:param kmer_dict: {kmer:uid}
	:param kmers: [{kmer:uid}]
	:return: The updated dictionary
	"""
	if kmers:
		for ku in kmers:
			kmer_dict[ku['kmer']] = ku['uid']

	return kmer_dict


def add_kmers_dgraph(client, all_kmers, genome):
	"""
	Add all kmers from a given genome to the graph
	:param client: dgraph client
	:param all_kmers: dict of lists of all kmers in genome dict[contig:[kmers]]
	:param genome: name of genome to add
	:return: None
	"""

	dict_values = all_kmers.values()
	kmers = list(dict_values)
	for x, i in enumerate(kmers):
		kmer_list = kmers[x]
		get_kmers_contig(kmer_list, client, genome)

def get_kmers_contig(ckmers, client, genome):
	"""
	Process a single contig into kmers, adding nodes and edges for each
	:param ckmers: The list of kmers for the contig
	:param client: dgraph client
	:param genome: indexed edge genome name
	:return: success
	"""
	# Query for all existing kmers
	kmer_uid_dict = {}
	kmer_uid_dict = add_kmers_dict(kmer_uid_dict, query_kmers_dgraph(client, ckmers))

	# Create list of kmers that need to be batch inserted into graph
	kmers_to_insert = []
	duplicates = []
	for kmer in ckmers:
		if kmer in kmer_uid_dict:
			duplicates.append(kmer)
		else:
			kmers_to_insert.append(kmer)

	add_metadata(client, duplicates, genome)

	if kmers_to_insert:
		# Bulk insert the kmers
		txn_result_dict = add_kmers_batch_dgraph(client, kmers_to_insert)
		# Update the dict of kmer:uid
		kmer_uid_dict = add_kmers_dict(kmer_uid_dict, txn_result_dict)

	# Batch the connections between the kmers
	# Creates a list of quads that need to be added later
	print('.', end='')
	return(add_edges_kmers(client, ckmers, kmer_uid_dict, genome))

def add_metadata(client, duplicates, genome):
	"""
	Bulk insert of duplicated kmers metadata. Connect the uid of the kmer that exists in the graph
	with the uid created by the metadata node. These uids are connected along a genome name edge
	allowing for a single kmer to have multiple duplicates accross genomes.
	:param client: the dgraph client
	:param duplicates: a list of duplicated kmers in a contig
	:param genome: the genome edge name that the kmers belong to
	"""
	print("before addition to schema")
	add_metadata_to_schema(client, genome)
	print("after addition to schema")
	bulk_quads = []

	for i in range(0, len(duplicates) - 2):
		bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(kmer_uid_dict[duplicates[i]],
														  genome,
														  get_metadata_uid(client),
														  "\n"
													  ))

	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()

def get_metadata_uid(client):
	"""
	Create the initial node for the metadata path. Similar to the kmer node this
	path contains a uid, a edge named "duplicate" and a string camed duplicate.
	The uid of this node will be connected to the uid of the duplicated kmer.
	Outwards from this node will be the pervious and next uids in the path.
	:param client: the dgraph client
	"""
	duplicate = "duplicate"
	bulk_quads.append('_:{0} <duplicate> "{0}" .{1}'.format(duplicate, "\n"))

	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()
		for uid in m.uids:
			print(uid)

	finally:
		txn.discard()

def add_edges_kmers(client, kmers, kmer_uid_dict, genome):
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
	for i in range(0, len(kmers) - 2):
		bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(kmer_uid_dict[kmers[i]],
														  genome,
														  kmer_uid_dict[kmers[i + 1]],
														  "\n"
													  ))


	#print(bulk_quads)
	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()


def add_kmers_batch_dgraph(client, kmer_list):
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
	#print(bulk_quads)

	# start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()
		# Create the output in the form that the program is expecting
		kmer_dict_list = []
		for uid in m.uids:
			kmer_dict_list.append({'kmer': uid, 'uid': m.uids[uid]})
	finally:
		txn.discard()
	return kmer_dict_list

def add_kmer_to_graph(client, ki, kn, genome):
	"""
	Every kmer needs to be linked to another kmer. Single kmers not permitted.
	:param client: dgraph client
	:param ki: the initial kmer. Get the uid if it exists. Otherwise create it.
	:param kn: the next kmer, linked to ki. Get the uid if it exists. Otherwise create it.
	:param genome: the genome label between the two kmers -- this needs to be previously indexed in the schema.
	:return: None
	"""

	# First we need to see if the two kmers exist already
	# This is an upsert. If it exists, get the uid.
	# Otherwise, create the uid.
	uid_ki = kmer_upsert(client, ki)
	uid_kn = kmer_upsert(client, kn)

	#start the transaction
	txn = client.txn()

	# The two kmers (as nodes) are linked by a <genome name> predicate
	try:
		d = {
			'uid': uid_ki,
			genome: {
				'uid': uid_kn
			}
		}
		txn.mutate(set_obj=d)
		txn.commit()

	finally:
		txn.discard()

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

def execute_args(client, opt):
	'''
	Initiates which commandline argument to execute based on a flag
	:param client: The dgraph client
	:param opt: Which commandline flag was given
	:param return: none
	'''
	if opt.insert:
		dgraph.insert_genome(client, opt.insert)
	if opt.query:
		dgraph.query_for_genome(client, opt.query)
	if opt.delete:
		dgraph.delete_genome(client, opt.delete)

def insert_genome(client, genomes):
	"""
	Function which inserts genome(s) into the graph based on a commandline argument.
	A list of genomes can be given as well
	:param client: The dgraph client
	:param genomes: Either a single genome or a list of genomes to be inserted
	:param return: none
	"""
	for genome in genomes:
		filename = os.path.abspath("data/genomes/insert/{}".format(genome))
		genome = "genome_" + compute_hash(filename)
		add_genome_to_schema(client, genome)
		all_kmers = kmer_from_file(filename, 11)
		add_all_kmers_to_graph(client, all_kmers, genome)
	print("inserted genome(s)")


def query_for_genome(client, genomes):
	"""
	Function which deletes genome(s) into the graph based on a commandline argument.
	A list of genomes can be given as well
	:param client: The dgraph client
	:param genomes: Either a single genome or a list of genomes to be deleted
	:param return: none
	"""
	for genome in genomes:
		filename = os.path.abspath("data/genomes/insert/{}".format(genome)) # TODO change pathing and figure out metadata querying
		genome = "genome_" + compute_hash(filename)
		sg1 = example_query(client, genome)
		print(sg1)

def delete_genome(client, genomes):
	"""
	Function which queries for genome(s) in the graph based on a commandline argument.
	A list of genomes can be given as well
	:param client: The dgraph client
	:param genomes: Either a single genome or a list of genomes to be queried for
	:param return: none
	"""
	for genome in genomes:
		genome = "genome_" + compute_hash(filename)

def create_graph(client, file, filepath):
	"""
	This function builds the graph by being repeatedly called by the main function.
	A path to a fasta file to be inserted is given and the function breaks the file up into its
	kmers and then adds those kmers to the graph with the genome name as an edge.
	:param client: The dgraph client
	:param file: The opened fasta file
	:param filepath: The absolute path to the fasta file which is being inserted
	"""
	try:
		LOG.info("Starting to create graph")
		x = 0
		filename = file.name
		genome = "genome_" + commandline.compute_hash(filepath)
		dgraph.add_genome_to_schema(client, genome)
		all_kmers = dgraph.get_kmers_files(filename, 11)

		dgraph.add_kmers_dgraph(client, all_kmers, genome)
		LOG.info("Finished creating the graph")
		sg1 = dgraph.example_query(client, genome)
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
		sg1 = path_query(client, genome)
		#print(sg1)
	except Exception as e:
		LOG.critical("Failed to create graph at file - {}".format(filename) + str(e))
		sys.exit()
