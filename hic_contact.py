import numpy as np

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pickle

class HiCContact():
    """
    HiCのcontact情報を保管しておくclass
    欲しい情報
    - 各readがどこにマップされてるか？の情報(self.r)
    - 各contigの長さ、配列情報(self.lengths, self.fasta)
    - contig名とidの対応表(self.ids, self.names)
    - contigの本数 (self.size)
    - uniqueなmapの情報 as pandas dataframe (self.df)
    """
    def __init__(self, sam_r1_fname, sam_r2_fname, fasta_fname):
        # r = {'hicのread id': {'r1': [contig_id, contig_id, ...], 'r2': []}}
        self.sam_r1_fname = sam_r1_fname
        self.sam_r2_fname = sam_r2_fname
        self.sam_r1 = pysam.AlignmentFile(sam_r1_fname, 'r')
        self.sam_r2 = pysam.AlignmentFile(sam_r2_fname, 'r')
        self.fasta_fname = fasta_fname
    
    def prepare_fasta(self):
        print('fasta')
        fasta = []
        lengths = []
        names = []
        for record in SeqIO.parse(self.fasta_fname, 'fasta'):
            fasta.append(record)
            lengths.append(len(record.seq))
            names.append(record.id)
        self.fasta = fasta
        self.lengths = lengths
        
        # functions that maps name to id of contig
        self.names = names
        self.ids = {}
        for i, name in enumerate(names):
            self.ids[name] = i
            
        # number of contigs
        self.size = len(names)
        
    def build(self):
        print('sam')
        r = {}
        for read in self.sam_r1:
            if not read.is_unmapped:
                if read.query_name not in r:
                    r[read.query_name] = {'r1':[], 'r2':[]}
                r[read.query_name]['r1'].append((read.reference_id, read.reference_start, read.is_reverse))
        for read in self.sam_r2:
            if not read.is_unmapped:
                if read.query_name in r:
                    r[read.query_name]['r2'].append((read.reference_id, read.reference_start, read.is_reverse))
        self.r = r
        
        # pandas as unique match
        d_df = {
            'readname':[],
            'r1_id':[],
            'r1_pos':[],
            'r1_reversed':[],
            'r2_id':[],
            'r2_pos':[],
            'r2_reversed':[],
        }
        for key, row in r.items():
            if len(row['r1']) == 1 and len(row['r2']) == 1:
                r1 = row['r1'][0]
                r2 = row['r2'][0]
                # unique
                d_df['readname'].append(key)
                d_df['r1_id'].append(r1[0])
                d_df['r1_pos'].append(r1[1])
                d_df['r1_reversed'].append(r1[2])
                #d_df['r1_length'].append(self.sam_r1.lengths[r1[0]])
                d_df['r2_id'].append(r2[0])
                d_df['r2_pos'].append(r2[1])
                d_df['r2_reversed'].append(r2[2])
                #d_df['r2_length'].append(self.sam_r2.lengths[r2[0]])
        self.df = pd.DataFrame.from_dict(d_df)
        
    def dump(self, filename):
        store = {
            'r': self.r,
            'df': self.df
        }
        with open(filename, mode='wb') as f:
            pickle.dump(store, f, protocol=4)
            
    def load(self, filename):
        with open(filename, mode='rb') as f:
            store = pickle.load(f)
        self.r = store['r']
        self.df = store['df']
        