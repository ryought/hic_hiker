import numpy as np

from tqdm import tqdm as tqdm
import argparse

import pandas as pd
import csv
import pysam

"""
Get contacts from mnd(merged nodups txt) or sam(sam r1 and sam r2).
result is composed of DataFrame, pandas.
"""

def get_contacts_mnd(contigs, mnd_filename):
    """get Hi-C contacts from mnd, or merged_nodups.txt, output of juicer(Aiden lab)"""
    X1, X2, P1, P2, U1, U2 = [], [], [], [], [], []
    N = 0
    with open(mnd_filename, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        for row in tqdm(reader):
            # x1, x2: contig id
            x1, x2 = contigs.get_id(row[1]), contigs.get_id(row[5])
            # p1, p2: mapping position. position in mnd(and samfile) is 1-origin, but df uses 0-origin.
            p1, p2 = int(row[2])-1, int(row[6])-1
            # row X1 and X2 must satisfy X1 <= X2
            if x1 <= x2:
                X1.append(x1)
                X2.append(x2)
                P1.append(p1)
                P2.append(p2)
                U1.append(True)
                U2.append(True)
            elif x1 > x2:
                X1.append(x2)
                X2.append(x1)
                P1.append(p2)
                P2.append(p1)
                U1.append(True)
                U2.append(True)
            N += 1
    print('processed', N, 'lines')
    df = pd.DataFrame(
        data={ 'X1':X1,'X2':X2, \
               'P1':P1,'P2':P2, \
               'U1':U1,'U2':U2 }
    )
    return df

def get_contacts_sam(contigs, sam_r1_filename, sam_r2_filename):
    """get Hi-C contacts from sam, which you can get by running 'bwa-mem' on contig fasta and Hi-C fastq."""
    sam_r1 = pysam.AlignmentFile(sam_r1_filename, 'r')
    sam_r2 = pysam.AlignmentFile(sam_r2_filename, 'r')
    size = sam_r1.nreferences
    R, X1, X2, P1, P2, U1, U2 = [], [], [], [], [], [], []
    N = 0
    for r1, r2 in tqdm(zip(sam_r1, sam_r2)):
        while r1.is_secondary or r1.is_supplementary:
            r1 = next(sam_r1)
            if R and R[-1] == r1.query_name:
                U1[-1] = False
        while r2.is_secondary or r2.is_supplementary:
            r2 = next(sam_r2)
            if R and R[-1] == r2.query_name:
                U2[-1] = False
        if r1.query_name != r2.query_name:  # TODO マップされてない時を考えてない？
            # ordering of samfile is corrupted.
            print('assertion failed')
            print(r1.query_name, r2.query_name)
            break
        if (not r1.is_unmapped) and (not r2.is_unmapped):
            R.append(r1.query_name)
            x1, x2 = contigs.get_id(r1.reference_name), contigs.get_id(r2.reference_name)
            p1, p2 = r1.reference_start, r2.reference_start
            if x1 <= x2:
                X1.append(x1)
                X2.append(x2)
                P1.append(p1)
                P2.append(p2)
                U1.append(True)
                U2.append(True)
            elif x1 > x2:
                X1.append(x2)
                X2.append(x1)
                P1.append(p2)
                P2.append(p1)
                U1.append(True)  # unique mapping
                U2.append(True)
        N += 1
    print(N, 'lines')
    df = pd.DataFrame(
        data={'R':R, \
              'X1':X1,'X2':X2, \
              'P1':P1,'P2':P2, \
              'U1':U1,'U2':U2}
    )
    return df

def test():
    from contigs import Contigs
    contigs = Contigs('test/small_mock/contigs.fasta')
    df = get_contacts_mnd(contigs, 'test/small_mock/mnd.txt')
    print(df)

def main():
    psr = argparse.ArgumentParser()
    psr.add_argument('--sam_r1', help='r1 samfile')
    psr.add_argument('--sam_r2', help='r2 samfile')
    psr.add_argument('--mnd', help='merged_nodups.txt from Juicer')
    psr.add_argument('fasta', help='contig fasta filename')
    psr.add_argument('output', help='output DataFrame(.feather) filename')
    args = psr.parse_args()

    print('parsing fasta')
    from contigs import Contigs
    contigs = Contigs(args.fasta)

    if args.sam_r1 and args.sam_r2:
        print('parsing samfile')
        df = get_contacts_sam(contigs, args.sam_r1, args.sam_r2)
        df.to_feather(args.output)
    elif args.mnd:
        print('parsing merged nodup txt')
        df = get_contacts_mnd(contigs, args.mnd)
        df.to_feather(args.output)
    else:
        print('Error: input files are not specified. Exit.')

if __name__ == '__main__':
    # main()
    test()
