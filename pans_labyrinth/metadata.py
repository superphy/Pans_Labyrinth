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
	print("addition to schema")
	add_metadata_to_schema(client, genome)
	bulk_quads = []
	metadata_uid = get_metadata_uid(client)
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(kmer_uid_dict[duplicates[0]],
														  genome,
														  metadata_uid,
														  "\n"
													  ))

	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()

	return metadata_uid

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
		print("metadata_uid", metadata_uid)

	finally:
		txn.discard()

	return metadata_uid

def connect_prev(client, duplicates, genome, metadata_uid):
	print("getting prev uid")
	kmer = duplicates[0]
	prev, next = get_next_prev_uid(client, kmer, genome)
	print(prev, next)
	bulk_quads = []
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(metadata_uid, "prev", prev, "\n"))

	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()


def connect_next(client, duplicates, genome, metadata_uid):
	print("geting next uid")
	kmer = duplicates[0]
	prev, next = get_next_prev_uid(client, kmer, genome)
	print(prev, next)
	bulk_quads = []
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(metadata_uid, "next", next, "\n"))

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
		genome(func: allofterms(kmer, {0})){{
	    {1}{{
	      uid
	      duplicate
	    }}
	    ~{1}{{
	    	uid
	    	duplicate
	  	}}
	  }}
	}}
	""".format(kmer, genome)

	res = client.query(query)
	p_res = json.loads(res.json)

	if p_res["genome"][0][genome][1]['duplicate'] == 'duplicate':
		next = p_res["genome"][0][genome][0]['uid']
	else:
		next = p_res["genome"][0][genome][1]['uid']

	prev = p_res["genome"][0]['~'+genome][0]['uid']

	return prev, next
