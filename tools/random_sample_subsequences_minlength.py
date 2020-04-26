#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
if len(sys.argv) != 6:
    print('not enough arguments')
    print('usage: python3 random_sample_subsequences.py input.fasta output.fasta <K: length of each window> <N: the number of windows> <T: minimum length threshold>')
    sys.exit()
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
# extract K-mer N subsequences
input_file = sys.argv[1]
output_file = sys.argv[2]
K = int(sys.argv[3])
N = int(sys.argv[4])
T = int(sys.argv[5])

origin_sequences = []
weights = []
ignored = []
for i, record in enumerate(SeqIO.parse(input_file, 'fasta')):
    if len(record.seq) >= T:
        origin_sequences.append(record)
        weights.append(len(record.seq) - K)
    else:
        ignored.append(record.id)
print('Ignored', len(ignored), 'contigs', ignored)

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


