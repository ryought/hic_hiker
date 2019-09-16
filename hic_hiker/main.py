#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
import os
import time

# from memory_profiler import profile

# @profile
def main():
    psr = argparse.ArgumentParser()
    psr.add_argument('fasta', help='assembler output fasta file containing contigs')
    psr.add_argument('assembly', help='assembler output fasta file containing contigs')
    psr.add_argument('mnd', help='mnd.txt file')
    psr.add_argument('workspace', help='workspace directory')
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
    print('time', time.time())
    contigs = contigs.Contigs(
            fasta_filename=args.fasta
            )

    print('loading', args.mnd)
    print('time', time.time())
    df = load.get_contacts_mnd(
            contigs=contigs,
            mnd_filename=args.mnd
            )
    print('loading', args.assembly)
    print('time', time.time())
    asm = load_3ddna.Assembly(
            asm_filename=args.assembly,
            contigs=contigs
            )

    original_layout = asm.get_layout()
    contigs2, df2 = asm.update_all(contigs, df)
    print('time', time.time())

    print('[step 2] estimating the contact probability distribution')
    estimator, raw_estimator = prob.infer_from_longest_contig2(
            df2,
            contigs2,
            K=args.K,
            K0=10**3.5
            )
    print('time', time.time())

    plt.figure()
    figures.fig_distribution(contigs, df, args.K)
    plt.savefig(workdir + 'fig_distribution.pdf', bbox_inches='tight', pad_inches=0.05)

    plt.figure()
    figures.fig_distribution(contigs, df, args.K, mode='log')
    plt.savefig(workdir + 'fig_distribution_log.pdf', bbox_inches='tight', pad_inches=0.05)

    print('[step 3] calculating emission probabilities')
    print('time', time.time())
    probs = prob.get_prob(
            df2,
            contigs2,
            original_layout,
            estimator,
            max_k=10,
            show_progress=True
            )

    print('[step 4] running optimization algorithm')
    print('time', time.time())
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
    print('time', time.time())
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
        print('time', time.time())
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
        plt.savefig(workdir + 'fig_errorchart.pdf', bbox_inches='tight', pad_inches=0.05)

        plt.figure()
        figures.fig_length_error(contigs2, original_layout, results[0], results[1])
        plt.savefig(workdir + 'fig_length_error.pdf', bbox_inches='tight', pad_inches=0.05)

        plt.figure()
        figures.fig_matrix(probs, contigs2, polished_layouts[3], results[3])
        plt.savefig(workdir + 'fig_matrix.pdf', bbox_inches='tight', pad_inches=0.05)

    print('finished all process')
    print('time', time.time())

if __name__ == '__main__':
    main()
