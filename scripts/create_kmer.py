import sys
from functions import *

try:
    infile = open(sys.argv[1])
    try:
        create_kmer(infile)
    except:
        print("failed to create dataframe")
except:
    print("failed to open file")
