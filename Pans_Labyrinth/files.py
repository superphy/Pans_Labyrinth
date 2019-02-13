#!/user/bin/env python
from files import *
from dgraph import *
from commandline import *
"""
Files and directory functions go here.
"""

def fill_graph_progess(client):
    '''
    Fills the graph with kmers and edges based on a given folder containing genomes
    Run time is arount 1 hour
    :param client: the dgraph client
    '''
    path = "data/genomes/test" # TODO change back to clean folder
    filecounter = 0
    x = 0
    for filepath in walkdir(path):
        filecounter += 1
    for filepath in tqdm(walkdir(path), total=filecounter, unit="files"):
        with open(filepath, 'rb') as file:
            filename = file.name
            genome = "genome_" + compute_hash(filename)
            print(genome)
            add_genome_to_schema(client, genome)
            all_kmers = get_kmers_files(filename, 11)
            kmers = all_kmers['SRR1122659.fasta|NODE_1_length_767768_cov_21.1582_ID_10270']
            #add_all_kmers_to_graph(client, all_kmers, genome)
            #print(kmers[0], kmers[1])
            add_kmer_to_graph(client, kmers[0], kmers[1], genome)
            sg1 = example_query(client, genome) 
            print(sg1)
