#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
contig class
"""
# import numpy as np
# import pysam
from Bio import SeqIO

class Contigs:
    def __init__(self, fasta_filename=None, names=None, ids=None, lengths=None, sequences=None, N=None):
        """
        Set of contig.
        each contig has its 'name'(identifier of fasta) and 'id'
        generated from fasta file, output of the assembler you use.
        """
        if fasta_filename and names and ids and lengths and sequences and N:
            self.fasta_filename = fasta_filename
            self.names = names
            self.ids = ids
            self.lengths = lengths
            self.sequences = sequences
            self.N = N
        elif fasta_filename:
            self._parse(fasta_filename)

    def __repr__(self):
        return "<Contigs (# {}) path:{} names:{} lengths:{}>".format(
                self.N, self.fasta_filename, self.names[:10], self.lengths[:10]
                )

    def _parse(self, fasta_filename):
        self.fasta_filename = fasta_filename

        # names and ids
        self.names = []
        self.ids = {}
        self.lengths = []
        self.sequences = []

        i = 0
        for record in SeqIO.parse(fasta_filename, 'fasta'):
            self.lengths.append(len(record.seq))
            self.names.append(record.id)
            self.sequences.append(record.seq)
            self.ids[record.id] = i
            i += 1

        # total number of contigs
        self.N = i
        print('loaded {0} contigs'.format(self.N))

    def get_id(self, name):
        return self.ids[name]

    def get_name(self, c_id):
        return self.names[c_id]

if __name__ == '__main__':
    c = Contigs('test/small_mock/contigs.fasta')
    print(c.names, c.ids, c.lengths)
    print(c.get_id(c.get_name(0)))
    print(c.N, 'contigs included')
    print(c)
