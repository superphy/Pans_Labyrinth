import pytest
from pans_labyrinth import main, files, dgraph, commandline, logging_functions

def inc(x):
    return x + 1


def test_query():
    assert inc(3) == 4
