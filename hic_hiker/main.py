#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
import os

from memory_profiler import profile

@profile
def main():
    psr = argparse.ArgumentParser()
    psr.add_argument('fasta', help='assembler output fasta file containing contigs')
    psr.add_argument('assembly', help='assembler output fasta file containing contigs')
    psr.add_argument('mnd', help='mnd.txt file')
    psr.add_argument('workspace', help='workspace directory')
    psr.add_argument('-n', default=100000, help='number of intra-contig contacts to be sampled for distribution estimation with KDE.', type=int)
    psr.add_argument('-w', default=10, help='KDE bandwidth parameter for distribution estimation', type=int)
    psr.add_argument('-K', default=75000, help='threshold', type=int)
    psr.add_argument('--refsam', help='sam file')
    psr.add_argument('--from', help='start from intermediate files', type=int)
    args = psr.parse_args()

    # print('sam', args.sam)
    print('path =', os.path.abspath(args.workspace), args.K, type(args.K))

    from . import contigs, load_3ddna, prob, hmm, load, benchmark, layout, figures

    # workdir = os.path.dirname(os.path.abspath(__file__))
    workdir = os.path.abspath(args.workspace) + '/'

    print('[step 1] Loading files')
    print('loading', args.fasta)
    contigs = contigs.Contigs(
            fasta_filename=args.fasta
            )

    print('loading', args.mnd)
    df = load.get_contacts_mnd(
            contigs=contigs,
            mnd_filename=args.mnd
            )
    print('loading', args.assembly)
    asm = load_3ddna.Assembly(
            asm_filename=args.assembly,
            contigs=contigs
            )

    original_layout = asm.get_layout()
    contigs2, df2 = asm.update_all(contigs, df)

    print('[step 2] estimating the contact probability distribution')
    estimator, raw_estimator = prob.infer_from_longest_contig(
            df2,
            contigs2,
            maxlength=args.K
            )

    plt.figure()
    figures.fig_distribution(contigs, df, args.K)
    plt.savefig(workdir + 'fig_distribution.png', dpi=100)

    print('[step 3] calculating emission probabilities')
    probs = prob.get_prob(
            df2,
            contigs2,
            original_layout,
            estimator,
            max_k=10,
            show_progress=True
            )

    print('[step 4] running optimization algorithm')
    polished_layouts = {}
    for k in [2,3,4,5]:
        print('k =', k)
        polished_layouts[k] = hmm.optimize_layout(
                probs,
                contigs2,
                original_layout,
                k=k
                )
    print('adaptive')
    polished_layouts[1] = hmm.optimize_layout(
            probs,
            contigs2,
            original_layout,
            K=args.K
            )

    print('[step 5] Writing FASTA and assembly')
    asm.generate_assembly_with_new_layout(
            workdir + 'polished.assembly',
            contigs2,
            polished_layouts[5]
            )
    layout.get_fasta(
            polished_layouts[5],
            contigs2,
            workdir + 'polished.fasta'
            )

    if args.refsam:
        print('[step 6] benchmarking')
        ref_layout = layout.get_reference_layout_from_sam(
                args.refsam,
                contigs2
                )
        results = {}
        results[0] = benchmark.determine_correct_orientation_or_not(contigs2, original_layout, ref_layout)[0]
        for k in [1,2,3,4,5]:
            results[k] = benchmark.determine_correct_orientation_or_not(contigs2, polished_layouts[k], ref_layout)[0]

        plt.figure()
        figures.fig_errorchart(results)
        plt.savefig(workdir + 'fig_errorchart.png', dpi=100)

        plt.figure()
        figures.fig_matrix(probs, contigs2, polished_layouts[3], results[3])
        plt.savefig(workdir + 'fig_matrix.png', dpi=100)

    print('finished all process')

if __name__ == '__main__':
    main()
