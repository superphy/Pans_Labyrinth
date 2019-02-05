#!/usr/bin/env python

"""
This is a simple example, using the layout of
https://github.com/dgraph-io/pydgraph/blob/master/examples/simple/simple.py
We have modified it for genome kmers
"""

import pydgraph
import json
import tempfile
import subprocess
import timeit
import os
import pandas as pd
from functions import *


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


def add_kmer_to_graph(client, ki, kn, genome, bulk=False):
    """
    Every kmer needs to be linked to another kmer. Single kmers not permitted.
    :param client: dgraph client
    :param ki: the initial kmer. Get the uid if it exists. Otherwise create it.
    :param kn: the next kmer, linked to ki. Get the uid if it exists. Otherwise create it.
    :param genome: the genome label between the two kmers -- this needs to be previously indexed in the schema.
    :param bulk: delay the insertion of the data, instead return a list of triples to be live-loaded
    :return: None
    """

    # First we need to see if the two kmers exist already
    # This is an upsert. If it exists, get the uid.
    # Otherwise, create the uid.
    uid_ki = kmer_upsert(client, ki)
    uid_kn = kmer_upsert(client, kn)

    # Data to be inserted
    d = "<{0}> <{1}> <{2}> .".format(uid_ki, genome, uid_kn)

    if not bulk:
        # start the transaction
        txn = client.txn()

        # The two kmers (as nodes) are linked by a <genome name> predicate
        try:
            txn.mutate(set_nquads=d)
            txn.commit()

        finally:
            txn.discard()
    else:
        return d

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
            # The data we wish to add in triple form
            d = """
                    _:blank-0 <kmer> "{0}" .
            """.format(kmer)

            m = txn.mutate(set_nquads=d)
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
    query = """
        {{
            genome(func: has({0})){{
                uid
                kmer
            }}
        }}
    """.format(genome)

    # This gets all but the last uid in the format {genome: [{'uid':'0x335'}, {'uid':'0x336'}]}
    res = client.query(query)
    j_res = json.loads(res.json)

    # In the graph A-p->B-p->C has(p) will return the uid for A and B, but not C
    # We want the last uid as well, so will construct a second query for it, using the last uid
    last_uid = j_res['genome'][-1]['uid']

    # Transaction for last query
    last_query = """
        {{
          lq(func: uid({0})){{
            {1}{{
                uid
                kmer
            }}
          }}
        }}
    """.format(last_uid, genome)

    # This gets the last query
    l_res = client.query(last_query)
    j_l_res = json.loads(l_res.json)
    return j_res['genome'] + j_l_res['lq']
    #j_l_res['lq'][genome])


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


def run_subprocess(cmd):
    """
    Run a subprocess from python
    :param cmd: the command to run
    :return: The completed process
    """

    start_time = timeit.default_timer()
    comp_proc = subprocess.run(
        cmd,
        shell=False,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    if comp_proc.returncode == 0:
        elapsed_time = timeit.default_timer() - start_time
        print("Subprocess {} finished successfully in {:0.3f} sec.".format(cmd, elapsed_time))
        return comp_proc
    else:
        print("Error in subprocess. The following command failed: {}".format(cmd))
        exit("subprocess failure")


if __name__ == '__main__':
    stub = create_client_stub()
    client = create_client(stub)
    drop_all(client)
    set_schema(client)

    '''
    # manual addition to graph -- this would normally be functions
    add_genome_to_schema(client, "genomeA")
    add_genome_to_schema(client, "genomeB")

    # Add the following two manually
    add_kmer_to_graph(client, "AAAAAAAAAAA", "TTTTTTTTCCC", "genomeA")
    add_kmer_to_graph(client, "AAAAAAAAAAA", "TTTTTTTTCCC", "genomeB")

    # Add the following two in bulk via the live loader
    # Bulk loading is off by default, hence the fourth True argument
    d1 = add_kmer_to_graph(client, "AAAAAAAAAAA", "CGCGCGCGCCA", "genomeB", True)
    d2 = add_kmer_to_graph(client, "TTTTTTTTCCC", "ATGATGATGAT", "genomeA", True)
    '''
    directory = "data/kmers"
    character = '_'
    bulk_nquads = []
    i = 0
    j = 1
    for filename in os.listdir(directory):
        df = pd.read_pickle('data/kmers/{}'.format(filename))
        genome = "genome_" + split_string(filename, character)
        add_genome_to_schema(client, genome)
        while j <= 1000:
            d1 = add_kmer_to_graph(client, df[i], df[j], genome, True)
            bulk_nquads.append(d1)
            i += 1
            j += 1
        break
    print(bulk_nquads)
    # Write each nquad to a separate line in a temp file
    # We need to close the file for subprocess to see the contents
    # We need to set delete=False to prevent the file from being deleted when we close it
    tfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    for line in bulk_nquads:
        print(line)
        tfile.write(line + "\n")
    tfile.close()

    # Check the file
    cat_com = [
        "head", tfile.name
    ]
    my_cat = run_subprocess(cat_com)
    print(my_cat)

    # Use the live load feature to load all the nquads
    # Use the live load feature to load all the nquads
    command = [
        "dgraph", "live",
        "-r", tfile.name
    ]
    run_subprocess(command)

    #query by predicate, to see the links
    sg1 = example_query(client, genome)
    #sg2 = example_query(client, "genomeB")
    print(sg1)
    #print(sg2)

    print("All done")
