#!/usr/bin/env python

from pans_labyrinth import files, dgraph, commandline, loggingFunctions
import os
import logging

def main():
    """
    The program
    :return: success
    """
    path = os.path.abspath("data/genomes/test")
    output_directory = os.path.abspath("data/logger")
    print(path)
    # setup the application logging
    LOG = loggingFunctions.create_logger()

    fh = logging.FileHandler(os.path.join(output_directory, 'pans_labyrinth.log'), 'w', 'utf-8')
    fh.setLevel(logging.DEBUG)
    LOG.addHandler(fh)

    #options = commandline.arg_parser(client)
    #dgraph.execute_args(client, options)

    #LOG.debug(options)
    LOG.info("Starting pans_labyrinth")

    try:
        LOG.info("Creating client stub")
        stub = dgraph.create_client_stub()
        try:
            LOG.info("Creating client")
            client = dgraph.create_client(stub)
            try:
                LOG.info("Dropping existing graph")
                dgraph.drop_all(client)
                try:
                    LOG.info("Add schema to graph with client")
                    dgraph.add_schema(client)
                    try:
                        LOG.info("Starting to create graph")
                        for filepath in files.walkdir(path):
                            with open(filepath, 'rb') as file:
                                dgraph.create_graph(client, file, filepath)
                        LOG.info("Finished creating graph")
                    except:
                        LOG.critical("Failed to create graph at file - {}".format(file.name))
                except:
                    LOG.critical("Failed to add schema to graph")
            except:
                LOG.critical("Failed to drop previous graph")
        except:
            LOG.critical("Failed to create client")
    except:
        LOG.critical("Failed to create the client stub")

    stub.close()
    LOG.info("ALL DONE")


if __name__ == '__main__':
    main()
