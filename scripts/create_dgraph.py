#!/usr/bin/env python
import pandas as pd
import os
from functions import *

"""
This is a simple example, using the layout of
https://github.com/dgraph-io/pydgraph/blob/master/examples/simple/simple.py
We have modified it for genome kmers
"""

import pydgraph
import json


def create_client_stub():
    """
    This allows us to create as many stubs as we want easily
    without having to specify the address every time.
    This allows a single source to be easily changed.
    :return: A client stub
    """

    return pydgraph.DgraphClientStub('localhost:9080')


def create_client(client_stub):
    """
    This allows us to create multiple clients easily,
    by passing in the stub from create_client_stub()
    :param client_stub: connection for stub
    :return: dgraph client
    """
    return pydgraph.DgraphClient(client_stub)


def drop_all(client):
    """
    Wipe the database for a fresh start
    :param client: the dgraph client
    :return: dgraph response
    """

    return client.alter(pydgraph.Operation(drop_all=True))


def set_schema(client):
    """
    We can define whatever type we like for our schema.
    For this simple example we will have only kmers and edges as genome names.
    add_genome_to_schema() allows genomes to be added programmatically.
    :param client: dgraph client
    :return: The client altered via the schema set out here
    """
    schema = """
        kmer: string @index(exact) .
    """

    return client.alter(pydgraph.Operation(schema=schema))


def add_kmer_to_graph(client, ki, kn, genome):
    """
    Every kmer needs to be linked to another kmer. Single kmers not permitted.
    :param client: dgraph client
    :param ki: the initial kmer. Get the uid if it exists. Otherwise create it.
    :param kn: the next kmer, linked to ki. Get the uid if it exists. Otherwise create it.
    :param genome: the genome label between the two kmers -- this needs to be previously indexed in the schema.
    :return: None
    """

    # First we need to see if the two kmers exist already
    # This is an upsert. If it exists, get the uid.
    # Otherwise, create the uid.
    uid_ki = kmer_upsert(client, ki)
    uid_kn = kmer_upsert(client, kn)

    #start the transaction
    txn = client.txn()

    # The two kmers (as nodes) are linked by a <genome name> predicate
    try:
        d = {
            'uid': uid_ki,
            genome: {
                'uid': uid_kn
            }
        }
        txn.mutate(set_obj=d)
        txn.commit()

    finally:
        txn.discard()


def kmer_upsert(client, kmer):
    """
    Check the existence of kmer. If it exists, return the uid.
    If it doesn't exist, create it, then return the uid.
    :param client: the dgraph client
    :param kmer: kmer sequence
    :return: uid of kmer
    """

    uid = kmer_query(client, kmer)

    if not uid:
        # kmer does not exist in database. Create it.
        txn = client.txn()

        try:
            # The data we wish to add
            d = {'kmer': kmer}
            m = txn.mutate(set_obj=d)
            txn.commit()
            uid = m.uids['blank-0']

        finally:
            # Good practice according the the dgraph docs
            # No network overhead associated with it
            txn.discard()

    return uid


def kmer_query(client, kmer):
    """
    Strictly queries the database for a kmer.
    Returns the uid if present, None if not.
    :param client: the dgraph client
    :param kmer: the kmer we would like to check
    :return: uid or None
    """

    # Transaction for the query
    query = """
        query find_kmer($k :string){
            find_kmer(func: eq(kmer, $k)){uid}
        }
    """
    variables = {'$k': kmer}
    res = client.query(query, variables=variables)
    json_res = json.loads(res.json)

    if json_res['find_kmer']:
      return json_res['find_kmer'][0]['uid']
    else:
        return None


def example_query(client, genome):
    """
    Example of getting a subgraph from a predicate query
    :param client: dgraph client
    :param genome: genome of interest as string (needs to be indexed in schema)
    :return: sub-graph as json
    """

    # Transaction for the query
    query = "{test(func: has(" + genome + ")){uid kmer}}"
    res = client.query(query)

    return json.loads(res.json)


def add_genome_to_schema(client, genome):
    """
    Index the genome name as a predicate, so functions can be used on it when searching etc.
    If the genome name is not added to the schema, it will not be indexed, and functions can
    not be applied to it. We will use the upsert directive when modifying the schema, so that
    if it exists, nothing happens, but new names are added to the schema.
    :param client: dgraph client
    :param genome: genomeName
    :return: Altered dgraph
    """
    schema = """
        {0}: uid .
    """.format(genome)

    return client.alter(pydgraph.Operation(schema=schema))


if __name__ == '__main__':
    stub = create_client_stub()
    client = create_client(stub)
    drop_all(client)
    set_schema(client)

    directory = "data/kmers"
    character = '_'
    for filename in os.listdir(directory):
        df = pd.read_pickle('data/kmers/{}'.format(filename))
        genome = split_string(filename, character)
        print(genome)
    '''
    # manual addition to graph -- this would normally be functions
    add_genome_to_schema(client, "genomeA")
    add_genome_to_schema(client, "genomeB")

    add_kmer_to_graph(client, "AAAAAAAAAAA", "TTTTTTTTCCC", "genomeA")
    add_kmer_to_graph(client, "AAAAAAAAAAA", "CGCGCGCGCCA", "genomeB")
    add_kmer_to_graph(client, "TTTTTTTTCCC", "ATGATGATGAT", "genomeA")

    #query by predicate, to see the links
    sg1 = example_query(client, "genomeA")
    sg2 = example_query(client, "genomeB")
    print(sg1)
    print(sg2)

    print("All done")
    '''
