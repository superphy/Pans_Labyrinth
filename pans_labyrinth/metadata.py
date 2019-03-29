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


def add_metadata(client, kmer_uid_dict, duplicates, genome, ckmers):
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
	kmer = duplicates[0]

	metadata_uid = get_metadata_uid(client, kmer_uid_dict[kmer])
	duplicate_uid = create_multi_duplicate_edge(client, metadata_uid, kmer)

	index = ckmers.index(kmer)

	connect_prev(client, duplicate_uid, kmer_uid_dict[ckmers[index -1]])
	connect_next(client, duplicate_uid, kmer_uid_dict[ckmers[index + 1]])

def get_metadata_uid(client, kmer_uid):
	"""
	Create the initial node for the metadata path. Similar to the kmer node this
	path contains a uid, a edge named "duplicate" and a string called duplicate.
	The uid of this node will be connected to the uid of the duplicated kmer.
	Outwards from this node will be the pervious and next uids in the path.
	:param client: the dgraph client
	"""
	bulk_quads = []
	print("getting uid")
	bulk_quads.append('<{0}> <duplicate> _:{1} .{2}'.format(kmer_uid, "duplicate", "\n"))
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

def connect_prev(client, metadata_uid, prev):
	print("getting prev uid", prev)
	bulk_quads = []
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(metadata_uid, "prev", prev, "\n"))

	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()


def connect_next(client, metadata_uid, next):
	print("geting next uid", next)
	bulk_quads = []
	bulk_quads.append('<{0}> <{1}> <{2}> .{3}'.format(metadata_uid, "next", next, "\n"))

	# Start the transaction
	txn = client.txn()

	try:
		m = txn.mutate(set_nquads=''.join(bulk_quads))
		txn.commit()

	finally:
		txn.discard()

def create_multi_duplicate_edge(client, metadata_uid, kmer):

	print("getting duplicate uid")
	number = duplicate_query(client, kmer)
	bulk_quads = []

	bulk_quads.append('<{0}> <{1}> _:{1} .{2}'.format(metadata_uid, number, "\n"))
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

def duplicate_query(client, kmer):
	print("duplicate query")

	query = """
	{{
	  duplicate(func: allofterms(kmer, "{0}")){{
	    duplicate{{
			number
		}}
	  }}
	}}
	""".format(kmer)

	res = client.query(query)
	path_res = json.loads(res.json)

	if "number" not in path_res:
		return "zero"
