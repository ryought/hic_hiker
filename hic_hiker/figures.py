#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
script to generate figures for the paper
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from matplotlib_scalebar.scalebar import ScaleBar
from . import prob, layout, benchmark

# Figure
def fig_distribution(contigs, df, K):
    top = np.argsort(contigs.lengths)[::-1]
    K0 = 10**3.5
	# generate top three estimators
    estimator, raw_estimator_0 = prob.infer_from_contig2(df, contigs, top[0], K=K, K0=K0)
    _, raw_estimator_1         = prob.infer_from_contig2(df, contigs, top[1], K=K, K0=K0)
    _, raw_estimator_2         = prob.infer_from_contig2(df, contigs, top[2], K=K, K0=K0)

	# draw plot
    _fig_distribution(raw_estimator_0, raw_estimator_1, raw_estimator_2, estimator, K)

def _fig_distribution(raw_estimator_1st, raw_estimator_2nd, raw_estimator_3rd, estimator, K):
    """
    >>> fig_distribution()
    >>> plt.show()
    """
    width = K + 10000
    x = np.linspace(1, width, 500)

    plt.xlabel('Separation Distance $d$ (bp)', fontsize=14)
    plt.ylabel('Contact Log Probability $\log p(d)$', fontsize=14)

    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    plt.yscale('log')

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
        for scaf_id in range(20):
            ok, ng_ord, ng_ori = parse_result(results[k][scaf_id])
            T[i] += ok
            F[i] += ng_ori

	# calculate error rates
    T, F = np.array(T), np.array(F)
    X = F / (T+F) * 100

	# show barplot
    print(labels, X)
    plt.bar(labels, X)
    plt.ylabel('Local Error Rate (%)', fontsize=18)
    plt.tick_params(axis='x', labelsize=18, rotation=90)
    #plt.xticks(np.arange(k*4+2), names, rotation=90)

    plt.tight_layout()

def parse_result(result):
    ok     = len(list(filter(lambda x:x=='ok', result)))
    ng_ord = len(list(filter(lambda x:x=='order_error', result)))
    ng_ori = len(list(filter(lambda x:x=='orientation_error', result)))
    return ok, ng_ord, ng_ori







# Figure
def fig_lengtherror():
    pass









# Figure
def fig_matrix(probs, contigs, polished_layout, result):
	fig = plt.figure(figsize=(8, 8), dpi=300)

	# show selected matrixes
	plt.subplot(2, 2, 1)
	i = 200
	_inspect_with_size(probs, polished_layout, contigs, result, scaf_id=0, pos=i)

	plt.subplot(2, 2, 2)
	i = 20
	_inspect_with_size(probs, polished_layout, contigs, result, scaf_id=0, pos=i)

	plt.subplot(2, 2, 3)
	i = 138
	_inspect_with_size(probs, polished_layout, contigs, result, scaf_id=0, pos=i)

	plt.subplot(2, 2, 4)
	i = 173
	im = _inspect_with_size(probs, polished_layout, contigs, result, scaf_id=0, pos=i)

	# show colorbar (this parameter was set manually)
	fig.subplots_adjust(right=0.9)
	cbar_ax = fig.add_axes([0.93, 0.15, 0.02, 0.7])
	fig.colorbar(im, cax=cbar_ax)

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

def _inspect_with_size(probs, polished_layout, contigs, result, scaf_id, pos):
    i = pos  # short hand
    k = 5  # how many neighbor contigs on the matrix?
    unit_length = 10000  # 1 px corresponds to unit_length bp in the plot

    target = probs[scaf_id][2*i - k*2 : 2*i + (k+1)*2, 2*i - k*2 : 2*i + (k+1)*2]
    print(target.shape, len(polished_layout.scaffolds[scaf_id].order), len(result[scaf_id]))
    orientations = polished_layout.scaffolds[scaf_id].orientation[i-k:i+k+1]

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

