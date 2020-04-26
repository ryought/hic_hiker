import sys
import pysam
import numpy as np
import random
from collections import defaultdict

def parse_chr_table(filename):
    import csv
    scaf_to_ref = {}
    ref_to_scaf = {}
    with open(filename, 'r') as f:
        reader = csv.reader(f)
        # header = next(reader)
        for row in reader:
            if row[0][0] == "#":
                continue
            scaf_name, ref_name, is_reverse = row
            scaf_to_ref[scaf_name] = (ref_name, True if is_reverse == 'True' else False)
            ref_to_scaf[ref_name] = (scaf_name, True if is_reverse == 'True' else False)
    return scaf_to_ref, ref_to_scaf

def main():
    scaf_to_ref, ref_to_scaf = parse_chr_table(sys.argv[1])
    print(scaf_to_ref)
    infile = pysam.AlignmentFile("-", "r")
    sam = []
    sorted_sam = {k: [] for k,v in scaf_to_ref.items()}
    for r in infile:
        scaf_name, scaf_start, scaf_end = r.query_name.split(',')
        scaf_start, scaf_end = int(scaf_start), int(scaf_end)
        ref_name = r.reference_name
        ref_start = r.reference_start
        mapq = r.mapping_quality
        ori = r.is_reverse
        sam.append((ref_name, ref_start, scaf_name, scaf_start, mapq, ori))
        if scaf_name in scaf_to_ref and mapq >= 60:
            sorted_sam[scaf_name].append((ref_name, ref_start, scaf_name, scaf_start, mapq, ori))
    for scaf_name in sorted_sam.keys():
        sorted_sam[scaf_name].sort(key=lambda x: x[3])


    result_cover = calc_cover_ratio(sam, scaf_to_ref)
    result_anchor = calc_anchor(sam, scaf_to_ref)
    result_order = calc_order(sorted_sam, scaf_to_ref)
    result_local_order = calc_local_order(sorted_sam, scaf_to_ref)
    result_orientation = calc_orientation(sorted_sam, scaf_to_ref)
    print(result_cover[2], 'cover')
    print(result_anchor[2], 'anchor')
    print(result_order[2], 'order')
    print(result_local_order[2], 'local_order')
    print(result_orientation[2], 'orientation')
    print(result_cover[2])
    print(result_anchor[2])
    print(result_order[2])
    print(result_local_order[2])
    print(result_orientation[2])
    print(result_cover, result_anchor, result_order, result_local_order, result_orientation)

def calc_cover_ratio(sam, scaf_to_ref):
    """chunkがどのくらい含まれているか？(アセンブリのカバー率)"""
    refs = defaultdict(lambda: {"T": 0, "F": 0})
    for ref_name, ref_start, scaf_name, scaf_start, mapq, ori in sam:
        if mapq >= 60 and scaf_name in scaf_to_ref:
            refs[ref_name]["T"] += 1
        else:
            refs[ref_name]["F"] += 1
    T = sum([count['T'] for ref_name, count in refs.items()])
    F = sum([count['F'] for ref_name, count in refs.items()])
    return T, F, T/(T+F), refs

def calc_anchor(sam, scaf_to_ref):
    """anchoring、染色体を突き止めた度合い"""
    # scafs[scaf_name] = {T:0, F:0}
    scafs = defaultdict(lambda: {"T": 0, "F": 0})
    refs = defaultdict(lambda: {"T": 0, "F": 0})
    for ref_name, ref_start, scaf_name, scaf_start, mapq, ori in sam:
        if mapq >= 60 and scaf_name in scaf_to_ref:
            if scaf_to_ref[scaf_name][0] == ref_name:
                scafs[scaf_name]["T"] += 1
                refs[ref_name]["T"] += 1
            else:
                scafs[scaf_name]["F"] += 1
                refs[ref_name]["F"] += 1
    T = sum([count['T'] for scaf_name, count in scafs.items()])
    F = sum([count['F'] for scaf_name, count in scafs.items()])
    T2 = sum([count['T'] for ref_name, count in refs.items()])
    assert T == T2
    return T, F, T/(T+F), scafs, refs

def calc_order(sorted_sam, scaf_to_ref):
    # number of trials
    N = 100000
    # T+F = (the number of pairs in the same scaffold and mapped to the same reference)
    # T = (the numeber of pairs that ordered correctly out of them)
    scafs = defaultdict(lambda: {"T": 0, "F": 0})
    def picker():
        scaf_names = list(scaf_to_ref.keys())
        weights = np.array([len(scaf_to_ref[scaf_name]) for scaf_name in scaf_names])
        p = weights / np.sum(weights)
        return np.random.choice(scaf_names, p=p)
    # pick randomly
    # the number of picked chunks
    n = 0
    while n < N:
        # pick scaffold
        scaf_name = picker()
        ref_name_of_scaf, scaf_is_reversed = scaf_to_ref[scaf_name]
        # pick two chunks
        chunk1 = random.choice(sorted_sam[scaf_name])
        chunk2 = random.choice(sorted_sam[scaf_name])
        if chunk1 == chunk2:
            # unfortunately, chunk1 and chunk2 is same chunk
            continue
        else:
            # chunk1 is a different chunk from chunk2
            ref_name_1, ref_start_1, _, scaf_start_1, _, _ = chunk1
            ref_name_2, ref_start_2, _, scaf_start_2, _, _ = chunk2
            if ref_name_1 == ref_name_2 == ref_name_of_scaf:
                # assert ref_name_1 == ref_name_of_scaf
                # -----ref_start_1------ref_start_2-------->
                is_forward_in_ref = (ref_start_1 < ref_start_2)
                is_forward_in_scaf = (scaf_start_1 < scaf_start_2)
                # a pair of different chunks that mapped onto the same reference -> OK
                if not scaf_is_reversed and is_forward_in_ref == is_forward_in_scaf:
                    scafs[scaf_name]["T"] += 1
                elif scaf_is_reversed and is_forward_in_ref != is_forward_in_scaf:
                    scafs[scaf_name]["T"] += 1
                else:
                    scafs[scaf_name]["F"] += 1
                n += 1
    T = sum([count['T'] for scaf_name, count in scafs.items()])
    F = sum([count['F'] for scaf_name, count in scafs.items()])
    return T, F, T/(T+F), scafs

def calc_local_order(sorted_sam, scaf_to_ref):
    # T+F = (the number of adjacent pair in the same scaffold and correct reference)
    # T = (ordering is consistent with reference)
    # local metric
    scafs = defaultdict(lambda: {"T": 0, "F": 0})
    for scaf_name in sorted_sam.keys():
        ref_name_of_scaf, scaf_is_reversed = scaf_to_ref[scaf_name]
        for i in range(len(sorted_sam[scaf_name]) - 1):
            ref_name_a, ref_start_a, scaf_name_a, scaf_start_a, _, _ = sorted_sam[scaf_name][i]
            ref_name_b, ref_start_b, scaf_name_b, scaf_start_b, _, _ = sorted_sam[scaf_name][i+1]
            assert scaf_start_a <= scaf_start_b, '{} {}'.format(scaf_start_a, scaf_start_b)
            if ref_name_a == ref_name_b == ref_name_of_scaf:
                if (not scaf_is_reversed and ref_start_a < ref_start_b) or (scaf_is_reversed and ref_start_a > ref_start_b):
                    scafs[scaf_name]['T'] += 1
                else:
                    scafs[scaf_name]['F'] += 1
        # print(scaf_name, scaf_is_reversed, t, f)
    T = sum([count['T'] for scaf_name, count in scafs.items()])
    F = sum([count['F'] for scaf_name, count in scafs.items()])
    return T, F, T/(T+F), scafs

def calc_orientation(sorted_sam, scaf_to_ref):
    # T+F = (the number of chunks that is in correct scaffold)
    # T = (orientation is consistent with reference)
    # global metric
    scafs = defaultdict(lambda: {"T": 0, "F": 0})
    for scaf_name in sorted_sam.keys():
        ref_name_of_scaf, scaf_is_reversed = scaf_to_ref[scaf_name]
        for ref_name, ref_start, scaf_name, scaf_start, mapq, ori in sorted_sam[scaf_name]:
            if ref_name == ref_name_of_scaf:  # in correct scaffold?
                if (not scaf_is_reversed and not ori) or (scaf_is_reversed and ori):
                    scafs[scaf_name]['T'] += 1
                else:
                    scafs[scaf_name]['F'] += 1
    T = sum([count['T'] for scaf_name, count in scafs.items()])
    F = sum([count['F'] for scaf_name, count in scafs.items()])
    return T, F, T/(T+F), scafs

def calc_local_orientation(sorted_sam):
    # T+F = (the number of triples of chunks, all are in same scaffold and reference)
    # T = (the center one has consistent orientation with adjacent two chunks)
    # local metric
    T, F = 0, 0
    return

if __name__ == '__main__':
    main()
