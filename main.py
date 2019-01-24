#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity

from hic_contact import HiCContact
from prob import *

def main():
    psr = argparse.ArgumentParser()
    psr.add_argument('sam_r1', help='sam file of r1.fastq')
    psr.add_argument('sam_r2', help='sam file of r2.fastq')
    psr.add_argument('fasta', help='assembler output fasta file containing contigs')
    psr.add_argument('-n', default=100000, help='number of intra-contig contacts to be sampled for distribution estimation with KDE.', type=int)
    psr.add_argument('-w', default=10, help='KDE bandwidth parameter for distribution estimation', type=int)
    psr.add_argument('--min-length', help='minimum length of contig. If not specified, all contigs will be used.', type=int)
    psr.add_argument('--debug', help='enable debug output', action='store_true')
    args = psr.parse_args()

    if args.debug:
        print('debug!')

    # main function
    process(args.sam_r1, args.sam_r2, args.fasta, args.debug, arg.n, arg.w)

def process(samr1, samr2, fasta, debug, n, w):
    # hmm infer
    ans = HiCContact(samr1, samr2, fasta)
    ans.prepare_fasta()
    # ans.build()
    kde = estimate_distribution(ans, N=n, bandwidth=b, debug=False)

if __name__ == '__main__':
    main()
