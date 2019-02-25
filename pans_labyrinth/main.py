#!/usr/bin/env python

from pans_labyrinth import files, dgraph, commandline, logging_functions
import os
import logging

def main():
    """
    The program - The work horse
    :return: success
    """
    path = os.path.abspath("data/genomes/test/")
    output_directory = os.path.abspath("data/logger")
    print(path)
    # setup the application logging
    LOG = logging_functions.create_logger()

    fh = logging.FileHandler(os.path.join(output_directory, 'pans_labyrinth.log'), 'w', 'utf-8')
    LOG.addHandler(fh)

    #options = commandline.arg_parser(client)
    #dgraph.execute_args(client, options)

    #LOG.debug(options)
    LOG.info("Starting pans_labyrinth")
    stub = dgraph.create_client_stub()
    client = dgraph.create_client(stub)
    dgraph.drop_all(client)
    dgraph.add_schema(client)

    LOG.info("Starting to create graph")
    for filepath in files.walkdir(path):
        with open(filepath, 'rb') as file:
            dgraph.create_graph(client, file, filepath)

    stub.close()
    LOG.info("ALL DONE")


if __name__ == '__main__':
    main()
