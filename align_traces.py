#!/usr/bin/env python

from Bio import SeqIO
import sys
import zipfile as zf
import tempfile as tf
import os

'''
Accepts a path to a directory full of .ab1 traces
and returns a dictionary mapping name to sequence
'''
def get_traces(directory):
    items = os.listdir(directory)
    for item in items:
        with open(os.path.join(directory,item),mode='rb') as ab1:
            trace = SeqIO.parse(ab1,"abi").next()
            seq = str(trace.seq)
            print seq

if __name__ == '__main__':
    get_traces(sys.argv[1])