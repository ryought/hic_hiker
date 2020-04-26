import sys
import os
import pysam
import numpy as np
import matplotlib.pyplot as plt
import random
from collections import defaultdict, Counter
import argparse

def get_scaf_to_ref(sorted_sam):
    scaf_to_ref = {}
    for scaf_name in sorted_sam.keys():
        occurrence = Counter( [(ref_name, ori) for (ref_name, ref_start, scaf_name, scaf_start, mapq, ori) in sorted_sam[scaf_name]] )
        first  = occurrence.most_common()[0]
        # print(scaf_name, occurrence, first)
        scaf_to_ref[scaf_name] = first[0]
    return scaf_to_ref

def parse_chunk_sam(input_filename):
    # parse sam file
    infile = pysam.AlignmentFile(input_filename, "r")
    # sam = []
    sorted_sam = defaultdict(lambda: [])
    for r in infile:
        scaf_name, scaf_start, scaf_end = r.query_name.split(',')
        scaf_start, scaf_end = int(scaf_start), int(scaf_end)
        ref_name = r.reference_name
        ref_start = r.reference_start
        mapq = r.mapping_quality
        ori = r.is_reverse
        # sam.append((ref_name, ref_start, scaf_name, scaf_start, mapq, ori))
        if mapq >= 60:
            sorted_sam[scaf_name].append((ref_name, ref_start, scaf_name, scaf_start, mapq, ori))
    for scaf_name in sorted_sam.keys():
        sorted_sam[scaf_name].sort(key=lambda x: x[3])
    return sorted_sam

def main():
    psr = argparse.ArgumentParser()
    psr.add_argument('sam_filename', type=str, help='1k chunks mapped onto the reference')
    psr.add_argument('out_dir', type=str, help='directory where I will save the dotplot png(s)')
    args = psr.parse_args()

    sorted_sam = parse_chunk_sam(args.sam_filename)
    scaf_to_ref = get_scaf_to_ref(sorted_sam)

    for i, scaf_name in enumerate(sorted_sam.keys()):
        print(scaf_name, scaf_to_ref[scaf_name])
        #print(sorted_sam[scaf_name])
        filtered = [s for s in sorted_sam[scaf_name] if s[0] == scaf_to_ref[scaf_name][0]]
        print(len(filtered), len(sorted_sam[scaf_name]))
        x_ref = np.array([s[1] for s in filtered])
        x_scaf = np.array([s[3] for s in filtered])
        ori = np.array([1 if s[5] == True else 0 for s in filtered])
        plt.figure(dpi=200)
        plt.scatter(x_ref, x_scaf, c=ori, s=0.1)
        plt.xlabel('ref:{}'.format(scaf_to_ref[scaf_name]))
        plt.ylabel('scaf:{}'.format(scaf_name))
        plt.savefig(os.path.join(args.out_dir, '{}_{}.png'.format(scaf_to_ref[scaf_name][0], scaf_name)))
        plt.clf()

main()
