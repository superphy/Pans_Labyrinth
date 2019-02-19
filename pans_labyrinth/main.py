#!/usr/bin/env python

from pans_labyrinth import files, dgraph, commandline
import os

def main():
    """
    The program
    :return: success
    """
    path = os.path.abspath("data/genomes/test")
    print(path)

    stub = dgraph.create_client_stub()
    client = dgraph.create_client(stub)
    dgraph.drop_all(client)
    dgraph.add_schema(client)
    #options = commandline.arg_parser(client)
    #dgraph.execute_args(client, options) # TODO figure out a better spot for this
    for filepath in files.walkdir(path):
        with open(filepath, 'rb') as file:
            filename = file.name
            print(filepath, filename)
            genome = "genome_" + commandline.compute_hash(filepath)
            dgraph.add_genome_to_schema(client, genome)
            all_kmers = dgraph.get_kmers_files(filename, 11)
            kmers = all_kmers['SRR1122659.fasta|NODE_1_length_767768_cov_21.1582_ID_10270']
            #add_all_kmers_to_graph(client, all_kmers, genome)
            dgraph.add_kmer_to_graph(client, kmers[0], kmers[1], genome)
            sg1 = dgraph.example_query(client, genome)
            print(sg1)

    stub.close()
    print("All done")


if __name__ == '__main__':
    main()
