#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to generate figures for the paper
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from matplotlib_scalebar.scalebar import ScaleBar
from . import prob, benchmark

# Figure
def fig_distribution(contigs, df, K, K0=10**3.5, mode='default'):
    top = np.argsort(contigs.lengths)[::-1]
	# generate top three estimators
    estimator, raw_estimator_0 = prob.infer_from_contig2(df, contigs, top[0], K=K, K0=K0)
    _, raw_estimator_1         = prob.infer_from_contig2(df, contigs, top[1], K=K, K0=K0)
    _, raw_estimator_2         = prob.infer_from_contig2(df, contigs, top[2], K=K, K0=K0)

	# draw plot
    _fig_distribution(raw_estimator_0, raw_estimator_1, raw_estimator_2, estimator, K, mode)

def _fig_distribution(raw_estimator_1st, raw_estimator_2nd, raw_estimator_3rd, estimator, K, mode):
    """
    >>> fig_distribution()
    >>> plt.show()
    """
    width = K + 10000
    if mode == 'default':
        x = np.linspace(1, width, 500)
    elif mode == 'log':
        x = np.logspace(0, np.log10(width), 500)

    plt.xlabel('Separation Distance $d$ (bp)', fontsize=14)
    plt.ylabel('Contact Probability $p(d)$', fontsize=14)

    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.yscale('log')

    if mode == 'log':
        plt.xscale('log')

    plt.plot(x, np.exp(raw_estimator_1st(x)), label='1st longest contig', linewidth=1, alpha=0.9)
    plt.plot(x, np.exp(raw_estimator_2nd(x)), label='2nd longest contig', linewidth=1, alpha=0.9)
    plt.plot(x, np.exp(raw_estimator_3rd(x)), label='3rd longest contig', linewidth=1, alpha=0.9)
    plt.plot(x, np.exp(estimator(x)), label='smoothed $p(d)$')

	# grayed region
    plt.axvspan(K, width, facecolor='gray', alpha=0.2)

    plt.legend(fontsize=14)
    plt.tight_layout()





# Figure
def fig_errorchart(results):
    # for each scaffold...
    ks = [0,1,2,3,4,5]
    labels = ['3D-DNA', 'HiC-Hiker adaptive', 'HiC-Hiker k=2', 'HiC-Hiker k=3', 'HiC-Hiker k=4', 'HiC-Hiker k=5']
    T = [0 for _ in ks]
    F = [0 for _ in ks]

    for i in range(len(ks)):
        k = ks[i]
        T[i] = 0
        F[i] = 0
        for scaf_id in range(len(results[k])):
            if len(results[k][scaf_id]) > 20:
                ok, ng_ord, ng_ori = parse_result(results[k][scaf_id])
                T[i] += ok
                F[i] += ng_ori

	# calculate error rates
    T, F = np.array(T), np.array(F)
    X = F / (T+F) * 100

	# show barplot
    print(labels, X)
    plt.bar(labels, X)
    for l, x in zip(labels, X):
        plt.text(l, x, '{:.3f}'.format(x), ha='center', va='center', fontsize=13)
        # plt.text(l, x, 'hoge', ha='center', va='center', fontsize=15)
    plt.ylabel('Local Orientation Error (%)', fontsize=18)
    plt.tick_params(axis='x', labelsize=18, rotation=90)
    plt.tick_params(axis='y', labelsize=18)

    plt.tight_layout()

def parse_result(result):
    ok     = len(list(filter(lambda x:x=='ok', result)))
    ng_ord = len(list(filter(lambda x:x=='order_error', result)))
    ng_ori = len(list(filter(lambda x:x=='orientation_error', result)))
    return ok, ng_ord, ng_ori







# Figure
def fig_length_error(contigs, layout, result_3d, result_hic):
    results_with_length = []
    for scaf_id in range(len(layout.scaffolds)):
        for scaf_pos in range(len(layout.scaffolds[scaf_id].order)):
            cid = layout.scaffolds[scaf_id].order[scaf_pos]
            length = contigs.lengths[cid]
            a = result_3d[scaf_id][scaf_pos]
            b = result_hic[scaf_id][scaf_pos]
            results_with_length.append((a, b, length))
            # print((a, b, length))

    bins = np.logspace(np.log10(15000), 6, num=20, base=10)
    lens_erA = np.histogram([x[2] for x in results_with_length if x[0]=='orientation_error'], bins=bins)[0]
    lens_erB = np.histogram([x[2] for x in results_with_length if x[1]=='orientation_error'], bins=bins)[0]
    lens_all = np.histogram([x[2] for x in results_with_length if x[1]=='orientation_error' or x[1]=='ok'], bins=bins)[0]
    print('ok')
    print(lens_erA)
    print(lens_erB)
    print(lens_all)

    plt.xscale('log')
    plt.plot(bins[:-1], lens_erA / lens_all * 100, label='3D-DNA', marker='o')
    plt.plot(bins[:-1], lens_erB / lens_all * 100, label='HiC-Hiker', marker='o')
    plt.xlabel('Contig Length (bp)', fontsize=18)
    plt.ylabel('Local Orientation Error (%)', fontsize=18)
    plt.tick_params(labelsize=18)
    plt.legend(fontsize=18)
    plt.xlim(10**4, 10**6)


# Figure
def fig_matrix(probs, contigs, polished_layout, ori_layout, result):
    fig = plt.figure(figsize=(8, 8), dpi=300)

    # show selected matrixes
    plt.subplot(2, 2, 1)
    i = 520
    _inspect_with_size(probs, polished_layout, ori_layout, contigs, result, scaf_id=0, pos=i)
    plt.subplot(2, 2, 2)
    i = 806
    _inspect_with_size(probs, polished_layout, ori_layout, contigs, result, scaf_id=0, pos=i)
    plt.subplot(2, 2, 3)
    i = 1248
    _inspect_with_size(probs, polished_layout, ori_layout, contigs, result, scaf_id=0, pos=i)
    plt.subplot(2, 2, 4)
    i = 1338
    im = _inspect_with_size(probs, polished_layout, ori_layout, contigs, result, scaf_id=0, pos=i)

    # show colorbar (this parameter was set manually)
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    # show axis indicator
    plt.gcf().text(0.06, 0.98, '→ $j$')
    plt.gcf().text(0.0,  0.94, '→ $i$', rotation=270)

    # figure label
    offset = 0.02
    plt.gcf().text(0+offset, 1-offset, 'a', fontsize=15, fontweight='bold')
    plt.gcf().text(0.45+offset, 1-offset, 'b', fontsize=15, fontweight='bold')
    plt.gcf().text(0+offset, 0.51-offset, 'c', fontsize=15, fontweight='bold')
    plt.gcf().text(0.45+offset, 0.51-offset, 'd', fontsize=15, fontweight='bold')

def normalization_matrix(prob, orientations):
    X, Y = prob.shape
    M = np.zeros((X//2, Y//2))
    for i in range(X//2):
        for j in range(Y//2):
            # p: numerator
            oi = orientations[i]
            p = logsumexp([prob[2*i + oi, 2*j + 0], prob[2*i + oi, 2*j + 1]])
			# P: denom
            P = logsumexp([prob[2*i + Oi, 2*j + Oj] for Oi in [0,1] for Oj in [0,1]])
            M[i,j] = p-P
    return M

def matrix_by_range(mat, lengths, unit_length):
    """scale the matrix such that each cell of the matrix (say M_{ij}) have the height and width proportional to length[i] and length[j] respectively """
    assert mat.shape[0] == len(lengths)
    N = len(lengths)
    sizes = [l//unit_length + 1 for l in lengths]
    s = sum(sizes)
    M = np.zeros((s, s))

    for i in range(N):
        for j in range(N):
            X0 = sum(sizes[:i])
            Y0 = sum(sizes[:j])
            X1 = sum(sizes[:i+1])
            Y1 = sum(sizes[:j+1])
            M[X0:X1, Y0:Y1] = mat[i, j]
    return M

def _inspect_with_size(probs, polished_layout, ori_layout, contigs, result, scaf_id, pos):
    i = pos  # short hand
    k = 5  # how many neighbor contigs on the matrix?
    unit_length = 10000  # 1 px corresponds to unit_length bp in the plot

    target = probs[scaf_id][2*i - k*2 : 2*i + (k+1)*2, 2*i - k*2 : 2*i + (k+1)*2]
    print(target.shape, len(polished_layout.scaffolds[scaf_id].order), len(result[scaf_id]))
    orientations = [
        0 if polished_layout.scaffolds[scaf_id].orientation[x] == ori_layout.scaffolds[scaf_id].orientation[x] else 1
        for x in range(i-k, i+k+1)
    ]
    # orientations = polished_layout.scaffolds[scaf_id].orientation[i-k:i+k+1]
    print(polished_layout.scaffolds[scaf_id].orientation[i-k:i+k+1])
    print(ori_layout.scaffolds[scaf_id].orientation[i-k:i+k+1])
    print(orientations)

    M = normalization_matrix(target, orientations)
    lengths = [contigs.lengths[polished_layout.scaffolds[scaf_id].order[i+ind]] for ind in range(-k, k+1)]
    M2 = matrix_by_range(M, lengths, unit_length=unit_length)

	# show the matrix
    im = plt.imshow(np.exp(M2), cmap='bwr', interpolation=None, aspect=1.0)
    plt.clim(0, 1)

	# show the scalebar
    scalebar = ScaleBar(unit_length, label_formatter=lambda value, unit: '{} Kbp'.format(value)) # 1 pixel = 0.2 meter
    plt.gca().add_artist(scalebar)

	# show contig ids
    names = ['{} {}'.format('x' if result[scaf_id][x] == 'orientation_error' else '#' if result[scaf_id][x] == 'order_error' else '', x) for x in range(i-k, i+k+1)]
    sizes = [l//unit_length + 1 for l in lengths]
    ticks_locations = [sum(sizes[:i])-0.5 for i in range(len(sizes))]
    names_locations = [sum(sizes[:i])-0.5+(sizes[i]/2) for i in range(len(sizes))]
    plt.yticks(names_locations, names)
    plt.xticks(names_locations, names, rotation=90)

	# show the border line of each contig
    for loc in ticks_locations:
        plt.axvline(x=loc, color='black', linewidth=0.5, alpha=0.5)
        plt.axhline(y=loc, color='black', linewidth=0.5, alpha=0.5)

	# disable default ticks
    plt.tick_params(bottom=False, left=False)

    #plt.vlines(x=ticks_locations, ymin=0, ymax=M2.shape[0], color='black', linewidth=0.5)
    #plt.hlines(y=ticks_locations, xmin=0, xmax=M2.shape[0], color='black', linewidth=0.5)
    #plt.axhline(y=2.5, color='black', linewidth=1)
    #plt.title(f'scaffold {scaf_id}')
    #colorbar = plt.colorbar(orientation="vertical")
    #colorbar.set_label(r'$\sum_{\theta_i, \theta_j} P(R_{ij} \mid \theta_i, \theta_j)$', fontsize=10)

    # Minor ticks
    #ax.set_xticks(np.arange(-.5, 10, 1), minor=True);
    #ax.set_yticks(np.arange(-.5, 10, 1), minor=True);
    # Gridlines based on minor ticks
    #plt.grid(which='minor', color='w', linestyle='-', linewidth=2)
    #plt.grid(color='black', linestyle='-', linewidth=0.1)
    #plt.axes().yaxis.set_minor_locator(plt.MultipleLocator(1.0))
    #ax = plt.gca();
    #ax.set_xticks(np.arange(-.5, M2.shape[0], 10), minor=True)
    #ax.set_yticks(np.arange(-.5, M2.shape[0], 10), minor=True)
    #ax.grid(which='minor', color='gray', linestyle='dashed', linewidth=0.5)
    #plt.tick_params(which='minor', bottom=False, left=False)
    #plt.grid(color='black', linewidth=0.1)

    plt.tight_layout()
    #plt.axis('off')
    return im

