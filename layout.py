#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from Bio.SeqRecord import SeqRecord
from contigs import Contigs

class Scaffold:
    def __init__(self, N=None, order=None, orientation=None):
        """
        Scaffold.order = [id]
        Scaffold.id2order

        These must satisfy this relationship:
          order[x] = y      x-th contig in this scaffold is contig y
          id2order[y] = x   contig y is x-th contig in this scaffold

        Scaffold.orientation = [ori]
          ori = 0 if "<-" 1 if "->"
        """
        if N:
            self.N           = N
            self.order       = [i for i in range(self.N)]
            self.orientation = [0 for _ in range(self.N)]

        elif (order and orientation):
            self.order = order
            self.orientation = orientation
            # assertion
            assert len(self.order) == len(self.orientation)
            self.N = len(self.order)

        self._id2order = {}
        for i, x in enumerate(self.order):
            self._id2order[x] = i
    def id2order(self, contig_id):
        return self._id2order[contig_id]

    def __repr__(self):
        return "<Scaffold order:{} orientation:{} N:{}>".format(self.order, self.orientation, self.N)

class Layout:
    def __init__(self, N=None, orders=None, orientations=None):
        """
        Layout.scaffolds   list of Scaffold
        Layout.N           total number of contigs
        Layout.id2order
        """
        self.scaffolds = []
        if N:
            self.N = N
            self.scaffolds.append(Scaffold(N=N))
        elif (orders and orientations):
            self.N = 0
            for order, orientation in zip(orders, orientations):
                self.scaffolds.append(Scaffold(order=order, orientation=orientation))
                self.N += len(order)
        # construct mapper from contig_id to order
        # id2order[contig_id] = (scaffold_id, order)
        self._id2order = {}
        for i, scaf in enumerate(self.scaffolds):
            for j, contig_id in enumerate(scaf.order):
                assert contig_id not in self._id2order, \
                        'contig {} is used more than once in scaffolds'.format(contig_id)
                self._id2order[contig_id] = (i, j)
    def id2order(self, contig_id):
        """this returns scaffold id the contig belongs to, and its position in the scaffold"""
        return self._id2order[contig_id]

    def __repr__(self):
        return "<Layout {}>".format(self.scaffolds)


def get_fasta(layout, contigs, filename):
    """指定のレイアウトでcontigを並べたfastaを作る"""
    # TODO
    print('saved')

def main():
    s1 = Scaffold(N=10)
    print(s1)
    s2 = Scaffold(order=[0,3,2,1,4], orientation=[0,0,1,0,1])
    print(s2, s2.id2order(2))

    layout = Layout(orders=[[3, 4, 2], [1, 0, 5]], orientations=[[0, 1, 1], [0, 0, 0]])
    print(layout, layout.N)
    for i in range(layout.N):
        print(i, layout.id2order(i))

if __name__ == '__main__':
    main()
