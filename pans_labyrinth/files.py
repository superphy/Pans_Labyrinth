#!/user/bin/env python

import os.path

"""
Files and directory functions go here.
"""

def get_files(file_or_directory):
    """
    Creates a list of files from either the given file, or all files within the
    directory specified (where each file name is its absolute path).
    :param file_or_directory (str): file or directory name given on commandline
    :return: a list of all the files found.
    """

    files_list = []
    if file_or_directory:
        if os.path.isdir(file_or_directory):
            # Create a list containing the file names
            for root, dirs, files in os.walk(file_or_directory):
                for filename in files:
                    files_list.append(os.path.join(root, filename))
        # check if input is concatenated file locations
        elif ',' in file_or_directory:
            for filename in file_or_directory.split(','):
                files_list.append(os.path.abspath(filename))
        else:
            files_list.append(os.path.abspath(file_or_directory))

    if not files_list:
        exit("No files were found")

    sorted_files = sorted(files_list)
    return sorted_files


# This is more like a main function.
# In files, we only want to deal with file-related stuff in a modular fashion
# The function below has everything, including querying the graph, adding kmers etc.
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
