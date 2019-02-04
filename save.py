import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
from tqdm import tqdm_notebook as tqdm
import pysam
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def supercontig_from_layout(layout, order, fasta_filename, sam_filename, contig_name, out_filename=None):
    sam = pysam.AlignmentFile(sam_filename, 'r')
    contigs = {}
    for record in SeqIO.parse(fasta_filename, 'fasta'):
        #print(record.id, sam.get_tid(record.id))
        tid = sam.get_tid(record.id)
        contigs[tid] = record.seq
    #print(contigs)
    #print(contigs[0].reverse_complement())
    #print(contigs[0] + contigs[0].reverse_complement())
    
    supercontig = None
    for o, l in zip(order, layout):
        if not supercontig:
            if l == 0:
                supercontig = contigs[o]
            else:
                supercontig = contigs[o].reverse_complement()
        else:
            if l == 0:
                supercontig += contigs[o]
            else:
                supercontig += contigs[o].reverse_complement()
                
    if not out_filename:
        # 文字列を返す
        return supercontig
    else:
        # fastaで保存する
        supercontig_rec = SeqRecord(supercontig, contig_name, '', '')
        SeqIO.write([supercontig_rec], out_filename, 'fasta')
        print('saved')