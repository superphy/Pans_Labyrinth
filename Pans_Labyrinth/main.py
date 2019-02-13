from files import *
from dgraph import *
from commandline import *

def main():
    """
    The program
    :return: success
    """

    stub = create_client_stub()
    client = create_client(stub)
    # arg_parser(client)
    drop_all(client)
    add_schema(client)
    fill_graph_progess(client)

    print("All done")


if __name__ == '__main__':
    main()
