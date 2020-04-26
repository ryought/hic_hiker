#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
if len(sys.argv) != 6:
    print('not enough arguments')
    print('usage: python3 random_sample_subsequences_chrs.py input.fasta input.chrs output.fasta <K: length of each window> <N: the number of windows>')
    sys.exit()
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

def parse_chr_table(filename):
    import csv
    chrs = []
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        # header = next(reader)
        for row in reader:
            if row[0][0] == "#":
                continue
            chrs.append(row[0])
    return chrs

# extract K-mer N subsequences
input_file = sys.argv[1]
input_chrs_file = sys.argv[2]
output_file = sys.argv[3]
K = int(sys.argv[4])
N = int(sys.argv[5])

chrs = parse_chr_table(input_chrs_file)

# create `origin_sequences`
# each item is a chr that have an entry in `chrs` table
origin_sequences = []
weights = []
import sys
for i, record in enumerate(SeqIO.parse(input_file, 'fasta')):
    if record.id in chrs:
        origin_sequences.append(record)
        weights.append(len(record.seq) - K)

weights = np.array(weights)
selected = np.random.choice(range(len(origin_sequences)), N, p=weights / np.sum(weights))

fragments = []
for s_id in selected:
    start = np.random.randint(weights[s_id] + 1)
    end = start + K
    seq = origin_sequences[s_id].seq[start:end]
    record = SeqRecord(seq, origin_sequences[s_id].name + ',{},{}'.format(start, end), '', '')
    fragments.append(record)

SeqIO.write(fragments, output_file, 'fasta')


