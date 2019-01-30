import sys
from functions import *
import hashlib

try:
	infile = open(sys.argv[1])
	# Had to create hash here for now, gave BIG issues trying to pass the file name to the function
	try:
		BUF_SIZE = 65536
		sha1 = hashlib.sha1()
		with open(sys.argv[1], 'rb') as f:
			while True:
				data = f.read(BUF_SIZE)
				if not data:
					break
				sha1.update(data)
		hash = sha1.hexdigest()
	except:
		print("Failed to create hash")

	try:
		create_kmer(infile, hash)
		try:
			create_RC_kmer(infile, hash)
		except:
			print("Failed to create reverse compliment dataframe")
	except:
		print("Failed to create dataframe")
except:
	print("Failed to open file")
