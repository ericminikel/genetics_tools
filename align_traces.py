#!/usr/bin/env python

from Bio import SeqIO
import sys
import zipfile as zf
import tempfile as tf

'''
Accepts a path to a zip file full of .ab1 traces
and returns a dictionary mapping name to sequence
'''
def get_traces(zip_path):
    with zf.ZipFile(zip_path,mode='rb') as z:
        zipped_items = z.infolist()
        for item in zipped_items:
            with z.open(item,mode='r') as ab1:
                trace = SeqIO.parse(ab1,"abi").next()
                seq = str(trace.seq)
                print seq

if __name__ == '__main__':
    get_traces(sys.argv[1])