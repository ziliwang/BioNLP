#!/usr/bin/env python3
"""
A simple parse for the MeSH ASCII format, as downloadable from here:
on ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/d2015.bin

For a reference on MeSH headings (i.e. keys in the resulting map), see
http://www.nlm.nih.gov/mesh/elmesh99.pdf

Originally published on TechOverflow.net
"""
from collections import defaultdict

__author__ = "Uli Koehler"
__copyright__ = "Copyright 2015, Uli Koehler"
__license__ = "Apache License v2.0"
__version__ = "1.0"


def readMeSH(fin):
    """
    Given a file-like object, generates MeSH objects, i.e.
    dictionaries with a list of values for each qualifier.

    Example: {"MH": ["Acetylcysteine"]}

    """
    currentEntry = None
    for line in fin:
        line = line.strip()
        if not line:
            continue
        # Handle new record. MeSH explicitly marks this
        if line == "*NEWRECORD":
            # Yiel old entry, initialize new one
            if currentEntry:
                yield currentEntry
            currentEntry = defaultdict(list)
            continue
        # Line example: "MH = Acetylcysteine"
        key, _, value = line.partition(" = ")
        # Append to value list
        currentEntry[key].append(value)
    # If there is a non-empty entry left, yield it
    if currentEntry:
        yield currentEntry

if __name__ == "__main__":
    # Example of how to use readMeSH()
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    args = parser.parse_args()
    with open(args.file, "r") as infile:
        # readMeSH() yields MeSH objects, i.e. dictionaries
        for entry in readMeSH(infile):
            print(entry)
            # print(entry['MH'][0])
            # print('\n'.join([i.split('|')[0] for i in entry['ENTRY']]))
