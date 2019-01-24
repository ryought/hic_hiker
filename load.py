import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

from hic_contact import HiCContact
import pysam
import argparse
import pickle

"""
samをパースしてcontactsにするやつ
"""

# main function
def get_contacts(ans):
    """
    @params ans: <class HiCContact>
    """
    contacts = [[[] for y in range(ans.size)] for x in range(ans.size)]
    for readname, c in ans.r.items():
        if len(c['r1']) == 1 and len(c['r2']) == 1:
            i = c['r1'][0][0]
            j = c['r2'][0][0]
            if i < j:
                contacts[i][j].append((c['r1'][0][1], c['r2'][0][1]))
            elif i > j:
                contacts[j][i].append((c['r2'][0][1], c['r1'][0][1]))
    return contacts

def get_contacts_numpy(ans, contacts):
    """
    @params ans: <class HiCContact>
    """
    contacts_alt = [[None for y in range(ans.size)] for x in range(ans.size)]
    for i in range(ans.size):
        for j in range(i+1, ans.size):
            contacts_alt[i][j] = (np.array([x[0] for x in contacts[i][j]]), \
                                  np.array([x[1] for x in contacts[i][j]]))
    return contacts_alt

def get_internal_contacts(r1samfile, r2samfile):
    sam_r1 = pysam.AlignmentFile(r1samfile, 'r')
    sam_r2 = pysam.AlignmentFile(r2samfile, 'r')
    r = {}
    for read in sam_r1:
        if not read.is_unmapped:
            if read.query_name not in r:
                r[read.query_name] = {'r1':[], 'r2':[]}
            r[read.query_name]['r1'].append((read.reference_id, read.reference_start, read.is_reverse))
    for read in sam_r2:
        if not read.is_unmapped:
            if read.query_name in r:
                r[read.query_name]['r2'].append((read.reference_id, read.reference_start, read.is_reverse))
    inter = [ abs(c['r1'][0][1] - c['r2'][0][1]) for key, c in r.items() \
             if len(c['r1'])==1 and len(c['r2'])==1 \
             and c['r1'][0][0] == c['r2'][0][0]]
    return inter

if __name__ == '__main__':
    psr = argparse.ArgumentParser()
    psr.add_argument('sam_r1', help='sam file of r1.fastq')
    psr.add_argument('sam_r2', help='sam file of r2.fastq')
    psr.add_argument('fasta', help='contig fasta file')
    psr.add_argument('output_filename', help='output pickle filename')
    args = psr.parse_args()
    
    print('parsing')
    from hic_contact import HiCContact
    ans = HiCContact(args.sam_r1, args.sam_r2, args.fasta)
    ans.prepare_fasta()
    ans.build()
    contacts = get_contacts_numpy(ans, get_contacts(ans))
    
    # dump it
    print('dumping')
    with open(args.output_filename, mode='wb') as f:
        pickle.dump(contacts, f, protocol=4)
    