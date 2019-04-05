#!/usr/bin/env python

# Author: Gates Kempenaar

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


def add_metadata(client, kmer_uid_dict, kmer_list, genome, ckmers):
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

	for x in range(0, len(kmer_list)):
		kmer = kmer_list[x][1]
		print(kmer_uid_dict[kmer])
		metadata_uid = get_metadata_uid(client, kmer_uid_dict[kmer])
		connect_prev(client, metadata_uid, kmer_uid_dict[kmer_list[x][0]])
		connect_next(client, metadata_uid, kmer_uid_dict[kmer_list[x][2]])

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
	"""
	Connects the duplicated kmers previous kmer to the duplicate edge. Creates
	the previous edge so that a path can be followed from the duplicated kmer.
	:param client: The dgraph client
	:param metadata_uid: The uid of the duplicate edge duplicate -> duplicate_1(metadata_uid) -> prev
	:param prev: The uid of the duplicated kmers previous kmer
	"""
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
	"""
	Connects the duplicated kmers next kmer to the duplicate edge. Creates
	the next edge so that a path can be followed from the duplicated kmer.
	:param client: The dgraph client
	:param metadata_uid: The uid of the duplicate edge duplicate -> duplicate_1(metadata_uid) -> next
	:param prev: The uid of the duplicated kmers next kmer
	"""

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

def create_duplicate_edge(client, metadata_uid, kmer):
	"""
	Creates the edge which comes off of the initial duplicate edge.
	duplicate -> duplicate_1
	The number is based off of how many times the kmer is duplicated and creats
	separate edges for all of its respective paths.
	:param client: the dgraph client
	:param metadata_uid: the uid which connects to the duplicated kmer
	:param kmer: the duplicated kmer

	"""
	print("getting duplicate uid")
	bulk_quads = []

	bulk_quads.append('<{0}> <{1}> _:{1} .{2}'.format(metadata_uid, "duplicate", "\n"))
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
