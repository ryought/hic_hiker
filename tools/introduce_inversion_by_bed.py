#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import csv
from Bio import SeqIO
def main():
    if len(sys.argv) != 4:
        print('usage: python3 introduce_inversion_by_bed.py input.fasta inversions.bed output.fasta')
        sys.exit()
    input_fasta_filename  = sys.argv[1]
    input_bed_filename    = sys.argv[2]
    output_fasta_filename = sys.argv[3]

    # parse bed
    inversions = {}
    with open(input_bed_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            chr_name = row[0]
            start = int(row[1])
            end = int(row[2])
            if chr_name not in inversions:
                inversions[chr_name] = []
            inversions[chr_name].append((start, end))
    # print(inversions)

    # check if non-overlapping
    for chr_name in inversions.keys():
        prev_end = 0
        for (start, end) in inversions[chr_name]:
            if prev_end >= start:
                print('overlapping!')
                print(chr_name, start, end)

    sequences = []
    for record in SeqIO.parse(input_fasta_filename, 'fasta'):
        print(record.name)
        if record.name in inversions:
            for (start, end) in inversions[record.name]:
                record.seq = record.seq[:start] + record.seq[start:end] + record.seq[end:]
        sequences.append(record)
    SeqIO.write(sequences, output_fasta_filename, 'fasta')
    print('done!')

main()
