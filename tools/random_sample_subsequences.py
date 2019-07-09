#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
if len(sys.argv) != 5:
    print('not enough arguments')
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
# extract K-mer N subsequences
input_file = sys.argv[1]
output_file = sys.argv[2]
K = int(sys.argv[3])
N = int(sys.argv[4])

origin_sequences = []
weights = []
ids = []
for i, record in enumerate(SeqIO.parse(input_file, 'fasta')):
    origin_sequences.append(record)
    ids.append(i)
    weights.append(len(record.seq) - K)

weights = np.array(weights)
selected = np.random.choice(ids, N, p=weights / np.sum(weights))

fragments = []
for s_id in selected:
    start = np.random.randint(weights[s_id] + 1)
    end = start + K
    seq = origin_sequences[s_id].seq[start:end]
    record = SeqRecord(seq, origin_sequences[s_id].name + ',{},{}'.format(start, end), '', '')
    fragments.append(record)

SeqIO.write(fragments, output_file, 'fasta')


