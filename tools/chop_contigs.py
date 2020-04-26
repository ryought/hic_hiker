from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

def main():
    if len(sys.argv) != 4:
        print('usage: python3 chop_contigs.py hoge.assembly contigs.fasta output_chopped_contigs.fasta')
        sys.exit(1)
    assembly_filename = sys.argv[1]
    in_fasta = sys.argv[2]
    out_fasta = sys.argv[3]
    print('assembly file', assembly_filename)
    print('in_fasta', in_fasta)
    print('out_fasta', out_fasta)

    name_map = {}
    with open(assembly_filename) as f:
        for l in f:
            if l[0] == '>':
                record = l[1:].rstrip('\n').split(' ')
                name_map[int(record[1])] = (record[0], int(record[2]))

    original_contigs = {}
    for record in SeqIO.parse(in_fasta, 'fasta'):
        L = len(record.seq)
        original_contigs[record.id] = record.seq

    x = 0
    new_contigs = []
    for contig_id, (name, length) in name_map.items():
        l = name.split(':::')
        if len(l) > 1:
            target_contig = original_contigs[l[0]]
            if l[1] == 'fragment_1':
                # 最初だよ
                x = 0
                new_contigs.append((name, l[0], 0, length))
                x += length
            else:
                # 2番目以降
                new_contigs.append((name, l[0], x, x + length))
                x += length
        else:
            # 何もしなくて良い
            new_contigs.append((name, l[0], 0, length))

    contigs = []
    for contig in new_contigs:
        name = contig[0]
        orig = original_contigs[contig[1]]  # 対応する元のコンティグ
        seq = orig[contig[2]:contig[3]]
        contig = SeqRecord(seq, name, '', '')
        contigs.append(contig)

    SeqIO.write(contigs, out_fasta, 'fasta')

main()
