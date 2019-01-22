from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


# Function to create kmers from a contig sequence
# Input: fasta file
# Output: series with the sequences and contig number
def create_kmer(infile):
    sequence_data = []
    contig_data = []
    for seq_record in SeqIO.parse(infile, "fasta"):
        i = 0
        j = 11
        contig = seq_record.id
        my_seq_orig = seq_record.seq
        for x in my_seq_orig:
            if j == len(my_seq_orig) + 1:
                break
            my_seq = my_seq_orig[i:j]
            rc_my_seq = my_seq.reverse_complement()
            sequence_data.append(str(rc_my_seq))
            contig_data.append(contig)
            j += 1
            i += 1
    df = pd.Series(sequence_data, index = contig_data)
    df.to_pickle('data/kmers.pkl')

def create_node():
    df = pd.read_pickle('data/kmers.pkl')
    print(df)
