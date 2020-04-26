#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import csv
import pysam
# internal
from .contigs import Contigs

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
        return "<Scaffold name:{} order:{} orientation:{} N:{}>".format(self.name, self.order, self.orientation, self.N)

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
                # print(scaffold, isinstance(scaffold, Scaffold), type(scaffold), Scaffold, type(scaffold) is Scaffold)
                # self.scaffolds.append(scaffold)
                # if isinstance(scaffold, Scaffold):
                # if type(scaffold) is Scaffold:
                self.scaffolds.append(scaffold)
                # else:
                    # print('hoge', scaffold, type(scaffold))
                    # raise Exception('class Layout only accepts instances of Scaffold')
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

def get_layout_from_agp(agp_filename: str, contigs: Contigs) -> Layout:
    # assume that the agp file is sorted by scaffold name and its position.
    with open(agp_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        scaffolds, orders, orientations = [], [], []
        now_scaf_name, now_scaf_pos = None, None
        for row in reader:
            # initialize
            if not now_scaf_name:
                now_scaf_name = row[0]
            # when goes into the row of next scaffold
            if now_scaf_name != row[0]:
                # add the scaffold entry
                scaf = Scaffold(order=orders, orientation=orientations, name=now_scaf_name)
                scaffolds.append(scaf)
                orders, orientations = [], []
                now_scaf_name = row[0]
            component_type = row[4]
            if component_type == 'W':
                contig_name = row[5]
                contig_id = contigs.get_id(contig_name)
                contig_ori = 0 if row[8] == '+' else 1
                orders.append(contig_id)
                orientations.append(contig_ori)
                assert row[6] == '1'   # contig should be in full length
        else:
            scaf = Scaffold(order=orders, orientation=orientations, name=now_scaf_name)
            scaffolds.append(scaf)
    layout = Layout(scaffolds)
    return layout

def generate_assembly_from_layout(asm_filename, contigs, layout):
    print('from layout')
    text = ''
    # first paragraph: list all contig id and its length
    for i, (name, length) in enumerate(zip(contigs.names, contigs.lengths)):
        cid = contigs.ids[name]
        text += '>{} {} {}\n'.format(name, cid+1, length)
    # list all scaffold order and orientations
    for scaf in layout.scaffolds:
        scaf_text = ' '.join([('-' if orientation == 1 else '') + str(cid+1) for (cid, orientation) in zip(scaf.order, scaf.orientation)])
        text += scaf_text
        text += '\n'
    with open(asm_filename, mode='w') as f:
        f.write(text)

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
        if len(order) > 0 and len(orientation) > 0:
            scaffolds.append(Scaffold(order=order, orientation=orientation, name=chr_name))
    layout = Layout(scaffolds)
    return layout

def get_fasta(layout: Layout, contigs: Contigs, filename: str, add_gap=False):
    """
    generate sequences as fasta (from contigs/layout) and save them as `filename`
    If add_gap == True, N x 500 will be inserted between contigs
    """
    output_scaffolds = []
    for (scaf_id, scaf) in enumerate(layout.scaffolds):
        scaf_sequence = Seq('')  # 空seq
        for x, (i, ori) in enumerate(zip(scaf.order, scaf.orientation)):
            if add_gap and x != 0:
                # add Nx500 gap if this is not a first element
                scaf_sequence += Seq('N') * 500
            if ori == 1:
                # revcomp
                scaf_sequence += contigs.sequences[i].reverse_complement()
            elif ori == 0:
                scaf_sequence += contigs.sequences[i]
            else:
                print('ori does not match 0 or 1')
        if scaf.name:
            output_scaffolds.append(SeqRecord(scaf_sequence, scaf.name, '', ''))
        else:
            output_scaffolds.append(SeqRecord(scaf_sequence, 'hicscaffold_{:06}'.format(scaf_id), '', ''))
    SeqIO.write(output_scaffolds, filename, 'fasta')
    print('saved as', filename)

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
