#!/usr/bin/env python

import pydgraph
import json

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

	schema = """
	prev: uid .
	"""

	client.alter(pydgraph.Operation(schema=schema))

	schema = """
	next: uid .
	"""

	client.alter(pydgraph.Operation(schema=schema))


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


def add_metadata(client, kmer_uid_dict, duplicates, genome):
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

	#for i in range(0, len(duplicates)):
		#print("running", i)
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(kmer_uid_dict[duplicates[0]],
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
	path contains a uid, a edge named "duplicate" and a string called duplicate.
	The uid of this node will be connected to the uid of the duplicated kmer.
	Outwards from this node will be the pervious and next uids in the path.
	:param client: the dgraph client
	"""
	duplicate = "duplicate"
	bulk_quads = []
	print("getting uid")
	bulk_quads.append('_:{0} <duplicate> "{0}" .{1}'.format(duplicate, "\n"))
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()
		for uid in m.uids:
			metadata_uid = m.uids[uid]

	finally:
		txn.discard()

	return metadata_uid

def connect_prev_next(client, duplicates, genome):
	kmer = duplicates[0]
	prev, next = get_next_prev_uid(client, kmer, genome)
	bulk_quads = []
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(prev, "next", next, "\n"))

	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()


def connect_metadata_to_path(client):

	bulk_quads = []
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(kmer_uid_dict[duplicates[0]],
														  "prev",
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


def get_next_prev_uid(client, kmer, genome):
	print("get_next_prev_uid")
	query = """
	{{
		prev(func: allofterms(kmer, {0})){{
	    {1}{{
			uid
		}}
	  }}
	}}
	""".format(kmer, genome)

	res = client.query(query)
	p_res = json.loads(res.json)
	print(p_res)
