#!/usr/bin/env python

from pans_labyrinth import files, dgraph, commandline

def main():
    """
    The program
    :return: success
    """

    stub = dgraph.create_client_stub()
    client = dgraph.create_client(stub)
    # arg_parser(client)
    dgraph.drop_all(client)
    dgraph.add_schema(client)
    #fill_graph_progess(client)

    stub.close()
    print("All done")


if __name__ == '__main__':
    main()
