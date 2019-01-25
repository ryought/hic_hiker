#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mockのcontigを作るスクリプト
"""
from sklearn.neighbors.kde import KernelDensity
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def estimate_contig_length_distribution(prob_fasta, debug=False):
    from sklearn.neighbors.kde import KernelDensity
    contig_lengths = []
    import csv
    with open(prob_fasta, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in reader:
            contig_lengths.append(int(row[1]))
            
    if debug:
        contig_lengths.sort(reverse=True)
        contig_lengths = np.array(contig_lengths)
        contig_cum = np.cumsum(contig_lengths)
        x = np.arange(len(contig_lengths))
        plt.subplot(3, 1, 1)
        plt.scatter(contig_lengths, contig_cum)
    
    X = contig_lengths.reshape((-1, 1))
    bw = 200
    kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(X)
    #sample = kde.sample(1000, random_state=0).reshape(-1)
    return kde

def split_fasta(in_fasta, out_fasta, kde):
    contigs = []
    for record in SeqIO.parse(in_fasta, 'fasta'):
        print(record.id)
        lengths = np.round(kde.sample(100000).reshape(-1))
        lengths = lengths[np.where(lengths > 50)]
        L = len(record.seq)
        assert sum(lengths) >= L
        x = 0
        for l in lengths:
            length = int(l)
            seq = record.seq[x:x+length]
            name = f'{record.id}_{x}_{x+length}'
            contig = SeqRecord(seq, name, '', '')
            contigs.append(contig)
            x += length
            if x > L:
                break
    if out_fasta:
        SeqIO.write(contigs, out_fasta, 'fasta')
    return contigs

def main():
    in_fasta = 'analysis/vc2010/vc2010.draft-20180405.ref.fasta'
    out_fasta = 'analysis/mock/mock.fasta'
    prob_fasta = 'analysis/abyss/denovo-contigs.tsv'
    kde = estimate_contig_length_distribution(prob_fasta)
    _ = split_fasta(in_fasta, out_fasta, kde)

if __name__ == '__main__':
    main()
