#!/user/bin/env python
from pans_labyrinth import files, dgraph, commandline
import argparse
import sys
import hashlib

"""
Arg-parser goes here
"""

def compute_hash(filename):
	'''
	Takes a path to a fasta file and creates a hash of the file to be
	used as the genome edge name
	:param filename: path to the fasta file to be hashed
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


def arg_parser(client):
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--insert", action = 'append', help = "Insert a new genome into the graph using a fasta file",)
	parser.add_argument("-q", "--query", action = 'append', help = "Find genome path in the graph based on the fasta file hash")
	parser.add_argument("-d", "--delete", action = 'append', help = "Remove a genome grom the graph by using a the fasta file hash")

	opt = parser.parse_args()
	return opt
