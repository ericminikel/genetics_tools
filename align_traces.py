#!/usr/bin/env python

from Bio import SeqIO
from Bio import pairwise2 as pw2
import sys
import zipfile as zf
import tempfile as tf
import os
import csv
import string

'''
Accepts a path to a directory full of .ab1 traces
and returns a dictionary mapping name to sequence
'''
def get_traces(directory):
    traces = {}
    items = os.listdir(directory)
    for item in items:
        with open(os.path.join(directory,item),mode='rb') as ab1:
            trace = SeqIO.parse(ab1,"abi").next()
            traces[trace.id] = str(trace.seq)
    return traces

'''
Accepts a path to the CRISPR+primer file, infers
the position in the 96-well plate and returns
a list of spacers
'''
def get_spacers(design_file):
    spacers = []
    with open(design_file,mode='rb') as f:
        reader = csv.reader(f,delimiter='\t')
        header = reader.next()
        counter = 0
        for row in reader:
            if row is not None:
                rowdata = dict(zip(header,row))
                spacer = rowdata['spacer']
                spacers.append(spacer)
            else:
                spacers.append('')
            counter += 1
    return spacers

'''
Returns a list where 0-based indices map to plate positions.
For instance if using 96-well plates and 
filling first by row, second by column, then 
position 2 is A3 and position 12 is B1.
'''
def map_to_plate(first='row',format=96):
    mapping = []
    nrow = int((format/1.5)**0.5)
    ncol = int(nrow*1.5)
    letters = list(string.uppercase)
    if first=='row':
        for row in range(0,nrow):
            for col in range(0,ncol):
                mapping.append(letters[row]+str(col+1))
    elif first=='col':
        for col in range(0,ncol):
            for row in range(0,nrow):
                mapping.append(letters[row]+str(col+1))
    return mapping

'''
Accepts a dictionary mapping names to Sanger traces,
where the names start with plate positions such as 'A11',
and a dictionary mapping plate positions to desired sequences.
Prints pairwise alignments of traces with desired seqs.
'''
def align_traces(traces,sequences):
	for trace_name in traces.keys():
		plate_position = trace_name.split('-')[0]
		trace_seq = traces[trace_name]
		desired_seq = sequences[plate_position]
		best_alignment = pw2.align.globalms(trace_seq,desired_seq,2,-1,-10,-5)[0]
		print pw2.format_alignment(*best_alignment)
		print plate_position, trace_name, desired_seq, trace_seq
	return

if __name__ == '__main__':
    traces = get_traces(sys.argv[1])
    spacers = get_spacers(sys.argv[2])
    mapping = map_to_plate('row',96)
    mapped_spacers = dict(zip(mapping,spacers))
    align_traces(traces,mapped_spacers)