#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
from contigs import Contigs
import pandas as pd
class Assembly:
    def __init__(self, asm_filename, contigs, delimiter=':::'):
        """
        class of 3d-dna output *.assembly
        Assembly.names = [contig_names]
        Assembly.lengths = [int]
        Assembly.ids = {contig_ids: contig_names} such that ids[names[x]] = x
        Assembly.scaffolds = [(contig_id, orientation \in {0, 1})]

        Assembly.new_contigs = [(new_id, new_name, original_name, start, end)]
        """
        names = []
        lengths = []
        ids = {}
        scaffolds = []

        with open(asm_filename) as f:
            for line in f:
                if line[0] == '>':
                    record = line[1:].rstrip('\n').split(' ')
                    # record = [contig_name, index, contig_length]
                    names.append(record[0])
                    lengths.append(int(record[2]))
                    ids[record[0]] = int(record[1]) - 1
                else:
                    # line = '-index index index -index index'
                    ls = line.rstrip('\n').split(' ')
                    scaffold = []
                    for contig in ls:
                        if contig[0] == '-':
                            # backward direction
                            scaffold.append((int(contig[1:]) - 1, 1))
                        else:
                            # forward direction
                            scaffold.append((int(contig) - 1, 0))
                    scaffolds.append(scaffold)
        self.names = names
        self.lengths = lengths
        self.ids = ids
        self.scaffolds = scaffolds

        new_contigs = []
        for contig_id, (contig_name, contig_length) in enumerate(zip(names, lengths)):
            ls = contig_name.split(delimiter)
            if len(ls) > 1:
                if ls[1] == 'fragment_1':
                    new_contigs.append((contig_id, contig_name, ls[0], 0, contig_length))
                    x = contig_length
                else:
                    new_contigs.append((contig_id, contig_name, ls[0], x, x + contig_length))
                    x += contig_length
            else:
                new_contigs.append((contig_id, contig_name, ls[0], 0, contig_length))
        self.new_contigs = new_contigs

        # original contig idがiのものは、新しいcontig idのどこに対応しているか？
        spliter = [[] for _ in range(contigs.N)]
        for contig_id, contig_name, old_contig_name, start, end in new_contigs:
            old_contig_id = contigs.get_id(old_contig_name)
            spliter[old_contig_id].append((contig_id, start, end))
        self.spliter = spliter

    def get_new_pos(self, contig_id, contig_pos):
        ls = self.spliter[contig_id]
        if len(ls) == 1:
            # 分断されていない
            return ls[0][0], contig_pos  # (new_contig_id, contig_pos)
        else:
            for (contig_id_new, start, end) in ls:
                if start <= contig_pos < end:
                    return contig_id_new, contig_pos - start
            return contig_id_new, end - 1

    def get_layout(self):
        from layout import Layout
        order       = [[o   for (o, ori) in scaffold] for scaffold in self.scaffolds]
        orientation = [[ori for (o, ori) in scaffold] for scaffold in self.scaffolds]
        layout = Layout(
                orders=order,
                orientations=orientation,
                )
        return layout

    def get_new_contigs(self, old_contigs: Contigs):
        sequences = [[old_contigs.sequences[contig_id][start:end] for (_, start, end) in new_contigs]
                for (contig_id, new_contigs) in enumerate(self.spliter)]
        new_contigs = Contigs(
                fasta_filename=old_contigs.fasta_filename,
                names=self.names,
                ids=self.ids,
                lengths=self.lengths,
                sequences=sequences,
                N=len(self.names))
        return new_contigs

    def get_new_df(self, df):
        U1, U2 = df['U1'].values, df['U2'].values
        X1, X2, P1, P2 = df['X1'].values, df['X2'].values, df['P1'].values, df['P2'].values
        L = len(df)
        for i in range(L):
            x1, x2, p1, p2 = X1[i], X2[i], P1[i], P2[i]
            nx1, np1 = self.get_new_pos(x1, p1)
            nx2, np2 = self.get_new_pos(x2, p2)
            if nx1 <= nx2:
                X1[i], X2[i], P1[i], P2[i], U1[i], U2[i] = nx1, nx2, np1, np2, U1[i], U2[i]
            else:
                X1[i], X2[i], P1[i], P2[i], U1[i], U2[i] = nx2, nx1, np2, np1, U2[i], U1[i]
        df_new = pd.DataFrame(
            data={ \
                  'X1':X1,'X2':X2, \
                  'P1':P1,'P2':P2, \
                  'U1':U1,'U2':U2 \
                  }
        )
        return df_new

def update_with_3ddna_output(asm_filename, contigs, df):
    """
    3d-dna cuts contigs and line them, so contigs and df are updated.

    """
    asm = Assembly(asm_filename, contigs)
    return contigs_new, df_new

if __name__ == '__main__':
    from contigs import Contigs
    contigs = Contigs('test/small_mock/contigs.fasta')
    asm = Assembly('test/small_mock/test.final.assembly', contigs)
    print(asm)

    for contig_id, contig_length in enumerate(contigs.lengths):
        for x in range(0, contig_length, 10):
            print((contig_id, x), asm.get_new_pos(contig_id, x))

    print(asm.get_layout())
    print(asm.get_new_contigs(contigs))

    contigs = Contigs('test/small_mock/contigs.fasta')
    from load import get_contacts_mnd
    df = get_contacts_mnd(contigs, 'test/small_mock/mnd.txt')
    print(df)

    print(asm.get_new_df(df))
