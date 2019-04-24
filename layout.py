#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
# from Bio.SeqRecord import SeqRecord
from contigs import Contigs
import pysam

class Scaffold:
    def __init__(self, N=None, order=None, orientation=None, name=''):
        """
        orderとorientationを管理する。それに加えてid -> orderの逆方向の関数も用意されている
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

        if name:
            self.name = name
        else:
            self.name = ''  # name
    def id2order(self, contig_id):
        return self._id2order[contig_id]
    def __repr__(self):
        # if self.name:
            # return "<Scaffold '{}' order:{} orientation:{} N:{}>".format(self.name, self.order, self.orientation, self.N)
        # else:
        return "<Scaffold order:{} orientation:{} N:{}>".format(self.order, self.orientation, self.N)

class Layout:
    def __init__(self, scaffolds=None):
        """
        Layout.scaffolds   list of Scaffold
        基本的にはscaffoldのarrayなんだけど、
        - Layout.N(全contig数(scaffold数ではなくて)を返す)
        - Layout.id2order(contig idから、どのscaffoldの何番目かを返す関数)
        が用意されていて便利
        """
        self.scaffolds = []
        if scaffolds:
            for scaffold in scaffolds:
                if isinstance(scaffold, Scaffold):
                    self.scaffolds.append(scaffold)
                else:
                    raise Exception('class Layout only accepts instances of Scaffold')
        self.update()
    def id2order(self, contig_id):
        """this returns (<scaffold id the contig belongs to>, <its position in the scaffold>)"""
        try:
            return self._id2order[contig_id]
        except KeyError as e:
            # TODO 何かしたほうがいいかもしれない
            raise KeyError(e)
            # return None
    def update(self):
        """update self.N and self.id2order along with self.scaffolds"""
        # construct mapper from contig_id to order
        # id2order[contig_id] = (scaffold_id, order)
        self._id2order = {}
        for i, scaf in enumerate(self.scaffolds):
            for j, contig_id in enumerate(scaf.order):
                assert contig_id not in self._id2order, \
                        'contig {} is used more than once in scaffolds'.format(contig_id)
                self._id2order[contig_id] = (i, j)
        # update N
        self.N = 0
        for scaffold in self.scaffolds:
            self.N += scaffold.N

    def __repr__(self):
        return "<Layout {}>".format(self.scaffolds)

def get_reference_layout_from_sam(sam_filename: str, contigs: Contigs) -> Layout:
    """contigをreferenceにmapした結果のsam fileから、正解のlayoutを作る"""
    sam = pysam.AlignmentFile(sam_filename, 'r')
    # temp database
    chrs = {}
    for chr_name in sam.references:
        chrs[chr_name] = []
    for r in sam:
        if not (r.is_unmapped or r.is_secondary or r.is_supplementary):
            contig_id = contigs.get_id(r.query_name)
            chrs[r.reference_name].append((r.pos, contig_id, r.is_reverse))
    scaffolds = []
    for chr_name, records in chrs.items():
        sorted_chr = sorted(records, key=lambda record: record[0])  # sort by its position
        order = [r[1] for r in sorted_chr]
        orientation = [1 if r[2] else 0 for r in sorted_chr]  # if is_reverse then 1 else 0
        scaffolds.append(Scaffold(order=order, orientation=orientation, name=chr_name))
    layout = Layout(scaffolds)
    return layout

def get_fasta(layout, contigs, filename):
    """指定のレイアウトでcontigを並べたfastaを作る"""
    # TODO
    print('saved')

def main():
    s = Scaffold(N=10)
    print(s)
    s1 = Scaffold(order=[5,8,7,6], orientation=[1,0,1,1], name='hoge')
    s2 = Scaffold(order=[0,3,2,1,4], orientation=[0,0,1,0,1], name='s2')
    print(s2, s2.id2order(2))

    layout = Layout([s1, s2])
    print(layout, layout.N)
    for i in range(layout.N):
        print(i, layout.id2order(i))

if __name__ == '__main__':
    main()
