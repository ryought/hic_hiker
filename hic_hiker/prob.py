import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
from scipy.optimize import curve_fit
# from tqdm import tqdm_notebook as tqdm
from tqdm import tqdm as tqdm
import pysam
import argparse
import pandas as pd
import feather

"""
確率を求める関係のやつ
df, contigs
"""

from .layout import Layout
from .contigs import Contigs

# new functions
def create_db(df: pd.DataFrame, contigs: Contigs):
    # db[x][y] = (i,j), x<=y, i<=j <=> rows df[i:j] are contacts between x and y
    db: dict = {}
    df.sort_values(['X1', 'X2'], inplace=True)
    # get raw values of X1 and X2 column
    X1, X2 = df.X1.values, df.X2.values
    N = len(X1)
    # accumlator
    x1, x2 = 0, 0
    pi = 0
    for i in range(N):
        if x1 != X1[i] or x2 != X2[i]:
            if x1 not in db:  # not in keys
                db[x1] = {}
            if x2 in db[x1]:  # duplication
                raise Exception('{} and {} appeared twice. Stopped.'.format(x1, x2))
            # add entry
            db[x1][x2] = (pi, i)
            # update accumlator
            x1, x2 = X1[i], X2[i]
            pi = i
    # fill the last segment
    if x1 not in db:
        db[x1] = {}
    db[x1][x2] = (pi, i)
    return db

def get_scaffold_full_prob(scaffold, contigs):
	emission = 0
	for i in range(scaffold.N):
		for j in range(i):
			# prob between j ... i
			oi = scaffold.orientation[i]
			oj = scaffold.orientation[j]
			Li = contigs.lengths[scaffold.order[i]]
			Lj = contigs.lengths[scaffold.order[j]]
			Lgap = sum([contigs. lengths[scaffold.order[x]] for x in range(j+1, i)])
			prob = get_emission_prob(P1, P2, db, estimator, scaffold.order[j], scaffold.order[i], Lj, Li, Lgap)
			emission += prob[oi*2 + oj]
	return emission

def get_emission_prob(P1, P2, db:dict, estimator, x:int, y:int, Lx:int, Ly:int, gapLength:int):
    # assume [x] ----- [y]
    # calc four probabilities (in case of 2x2 orientations) of x-th and y-th contigs (and given gapLength) efficiently
    # target data is df[i:j]
    if x <= y:
        if x in db and y in db[x]:
            i, j = db[x][y]
            px, py = P1[i:j], P2[i:j]
        else:
            px, py = np.array([]), np.array([])
    else:
        if y in db and x in db[y]:
            i, j = db[y][x]
            px, py = P2[i:j], P1[i:j]
        else:
            px, py = np.array([]), np.array([])
    # m is np.ndarray([[dists], [dists], [dists], [dists]])
    m = np.stack([
        Lx - px + py + gapLength,       # 0,0
        px + py + gapLength,            # 1,0
        Lx - px + Ly - py + gapLength,  # 0,1
        Ly - py + px + gapLength        # 1,1
        ])
    # return value is np.ndarray([P(0,0), P(1,0), P(0,1), P(1,1)])
    # P(ori_i,ori_j) = ret[2*ori_j + ori_i], ori_i,ori_j \in {0,1}
    return np.sum(estimator(m), axis=1)

# np.add.outerの3以上引数バージョン
def many_outer(xs):
	if len(xs) == 2:
		return np.add.outer(xs[0], xs[1]).reshape(-1)
	if len(xs) == 1:
		return xs[0]
	else:
		return np.add.outer(many_outer(xs[:-1]), xs[-1]).reshape(-1)

# a=(1,2,3,4) b=(5,6,7,8) take_alternate(a,b)=(1,5,2,6,3,7,4,8)
def take_alternate(a, b):
    assert a.shape == b.shape
    res = np.empty(a.shape[0] + b.shape[0], dtype='float32')
    res[0::2] = a
    res[1::2] = b
    return res

# x=(x1,x2,x3,x4) y=(y1,y2) merge_outer(x,y)=(x1+y1, x2+y1, x3+y2, x4+y2)
def merge_outer(x, y):
    assert len(x) == len(y) * 2
    x += np.repeat(y, 2)
    return x

def get_emission_vec(t, P1, P2, db, estimator, contigs, scaffold, perms, first=False):
    # tにおけるemissionを計算する
    # tは中心座標
    # permsはpermutation
    N = scaffold.N
    W = len(perms[0])
    if first:
        # add up all contacts within this state
        lst = []
        for perm in perms:
            pos = [t + offset for offset in perm]
            if any(map(lambda x: x<0, pos)) or any(map(lambda x: x>=N, pos)):
                z = np.ones(2**W, dtype='float32') * -np.inf
                lst.append(z)
            else:
                ids = [scaffold.order[p] for p in pos]
                lengths = [contigs.lengths[i] for i in ids]
                #print('first mode', pos, ids, lengths)
                accum = np.array([0, 0], dtype='float32')
                for right in range(1, len(ids)):
                    #print(right)
                    #print([(left, right, ids[left], ids[right], sum(lengths[left+1:right])) for left in range(0, right)])
                    parts = np.stack([
                        get_emission_prob(
                            P1, P2, db, estimator, ids[left], ids[right], lengths[left], lengths[right], sum(lengths[left+1:right])
                        ) for left in range(0, right)
                    ])
                    jointed = take_alternate(
                        many_outer(parts[:,0:2]),
                        many_outer(parts[:,2:4])
                    )
                    #print('join!', jointed, accum)
                    accum = merge_outer(jointed, accum)
                    #print('joined', accum)
                lst.append(accum)
        return np.concatenate(lst)
    else:
        lst = []
        for perm in perms:
            # 配置が決まった時、一番後とそれ以外の要素とのprobを全部計算。
            # pはtからのオフセットが書いてある
            pos = [t + offset for offset in perm]
            if any(map(lambda x: x<0, pos)) or any(map(lambda x: x>=N, pos)):
                z = np.ones(2**W, dtype='float32') * -np.inf
                lst.append(z)
            else:
                ids = [scaffold.order[p] for p in pos]
                lengths = [contigs.lengths[i] for i in ids]
                gaps = [sum(lengths[i+1:-1]) for i in range(len(ids)-1)]
                parts = np.stack([
                    get_emission_prob(P1, P2, db, estimator, ids[i], ids[-1], lengths[i], lengths[-1], gaps[i]) for i in range(len(ids)-1)
                ])
                #print(perm, pos, ids, lengths, gaps)
                #print(parts)
                #print(jointed)
                # 前半同士、後半同士をくっつける
                jointed = take_alternate(
                    many_outer(parts[:,0:2]),
                    many_outer(parts[:,2:4])
                )
                lst.append(jointed)
        #print(lst)
        return np.concatenate(lst)

def run(w, contigs, scaffold, df_sorted, db, estimator):
    perms, trans = preprocess(w)
    print(perms, trans)
    M = 2**(2*w+1)
    N = len(perms)
    n = scaffold.N  # TODO number of contigs
    P1, P2 = df_sorted.P1.values, df_sorted.P2.values
    viterbi = np.zeros((n, M*N), dtype='float32')
    source  = np.zeros((n, M*N), dtype='uint16')
    print('M*N=', M*N)
    for t in range(w, n-w):  # loop over contig number
        # emission `e`
        if t == w:
            emission_vec = get_emission_vec(t, P1, P2, db, estimator, contigs, scaffold, perms, first=True)
        else:
            emission_vec = get_emission_vec(t, P1, P2, db, estimator, contigs, scaffold, perms, first=False)
        print('t=', t, 'emission=', emission_vec)
        # transition `t`
        for j in range(N):  # loop over state number
            for k in range(M):
                q = k // 2
                transition = np.concatenate([trans[j] * M + q, trans[j] * M + q + M//2])
                print(j*M + k, trans[j], transition)
                picked = viterbi[t-1][transition]
                print('picked', picked)
                idx = np.argmax(picked)
                source[t][j*M+k]  = idx
                viterbi[t][j*M+k] = picked[idx]
                print('choose', picked[idx], viterbi[t][j*M+k])
        viterbi[t] += emission_vec
    return viterbi, source

import itertools

def preprocess(w):
	start = 0
	permut = list(filter(lambda xs: len(xs) == len(set(xs)) and start in xs,
		itertools.product(*[range(i-w, i+w+1) for i in range(start-w, start+w+1)])
		))
	# 先頭がどこにあるかを記録
	prefix_dict = {}
	for i, p in enumerate(permut):
		ip = tuple(map(lambda x: x-1, p[1:]))
		if ip not in prefix_dict:
			prefix_dict[ip] = []
		prefix_dict[ip].append(i)
	# transition
	trans = [ np.array(prefix_dict[p[:-1]], dtype='uint16') for p in permut]
	return permut, trans

def get_prob2(df: pd.DataFrame, contigs: Contigs, layout: Layout, estimator, max_k=None, remove_repetitive=False, show_progress=False):
    # calc for each probabilities
    index = [(0, 0), (0, 1), (1, 0), (1, 1)]
    lengths = contigs.lengths

    # prob matrixに変換する
    # probsは各scaffoldに対応するprob matrixを入れるやつ
    probs = []
    for scaf in layout.scaffolds:
        size = scaf.N
        probs.append(np.zeros((size * 2, size * 2)))

    P1, P2 = df.P1.values, df.P2.values
    db = create_db(df, contigs)
    for i, scaf in enumerate(layout.scaffolds):
        for x1 in range(scaf.N):
            for x2 in range(x1+1, scaf.N):
                if (max_k and abs(x1 - x2) > max_k):
                    # scaffoldが違うところに乗っていたり、
                    # order上でk以上離れている確率は使わないので無視する
                    pass
                else:
                    # position list
                    I1, I2 = scaf.order[x1], scaf.order[x2]
                    L1, L2 = lengths[I1], lengths[I2]
                    # length list
                    gap = sum([ lengths[scaf.order[x]] for x in range(min(x1, x2)+1, max(x1, x2))])

                    p = get_emission_prob(P1, P2, db, estimator, I1,I2,L1,L2,gap)

                    probs[i][2*x1 + 0, 2*x2 + 0] = p[0]
                    probs[i][2*x1 + 1, 2*x2 + 0] = p[1]
                    probs[i][2*x1 + 0, 2*x2 + 1] = p[2]
                    probs[i][2*x1 + 1, 2*x2 + 1] = p[3]
    return probs

def get_prob(df: pd.DataFrame, contigs: Contigs, layout: Layout, estimator, max_k=None, remove_repetitive=False, show_progress=False):
    if remove_repetitive:
        # use only read that two "unique" flag are set
        df = df[(df['U1'] == True) & (df['U2'] == True)]

    # calc for each probabilities
    index = [(0, 0), (0, 1), (1, 0), (1, 1)]
    lengths = contigs.lengths
    def calc_probs(d):
        # contig id
        I1, I2 = d.name[0], d.name[1]
        S1, X1 = layout.id2order(I1)
        S2, X2 = layout.id2order(I2)
        # contig place
        result = [0, 0, 0, 0]
        if (S1 != S2) or (max_k and abs(X1 - X2) > max_k):
            # scaffoldが違うところに乗っていたり、
            # order上でk以上離れている確率は使わないので無視する
            pass
        else:
            # scaffold
            scaf = layout.scaffolds[S1]  # = S2
            # position list
            P1, P2 = d.P1.values, d.P2.values
            # length list
            L1, L2 = lengths[I1], lengths[I2]
            gap = sum([ lengths[scaf.order[x]] for x in range(min(X1, X2)+1, max(X1, X2))])
            for i, (o1, o2) in enumerate(index):
                distances = get_dists(P1, P2, o1, o2, L1, L2, gap)
                result[i] = estimator(distances).sum()
        return pd.Series(result, index=index)

    # use DataFrame.groupby from pandas
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
    print('calcing version3')
    """
    if show_progress:
        from tqdm import tqdm
        tqdm.pandas()
        prob_df = df.groupby(['X1', 'X2']).progress_apply(calc_probs)
    else:
        prob_df = df.groupby(['X1', 'X2']).apply(calc_probs)
    """

    # prob matrixに変換する
    # probsは各scaffoldに対応するprob matrixを入れるやつ
    probs = []
    for scaf in layout.scaffolds:
        size = scaf.N
        probs.append(np.zeros((size * 2, size * 2)))

    P1, P2 = df.P1.values, df.P2.values
    groups = df.groupby(['X1', 'X2']).groups.items()
    print('iter version3')
    for (I1, I2), idx in tqdm(groups):
        try:
            S1, X1 = layout.id2order(I1)
            S2, X2 = layout.id2order(I2)
        except KeyError as e:
            # this contig is not used in the layout.
            continue
        # contig place
        if (S1 != S2) or (max_k and abs(X1 - X2) > max_k):
            # scaffoldが違うところに乗っていたり、
            # order上でk以上離れている確率は使わないので無視する
            pass
        elif X1 == X2:
            # ignore contacts on the same contig
            pass
        else:
            # scaffold
            scaf = layout.scaffolds[S1]  # = S2
            # position list
            p1, p2 = P1[idx], P2[idx]
            # length list
            L1, L2 = lengths[I1], lengths[I2]
            gap = sum([ lengths[scaf.order[x]] for x in range(min(X1, X2)+1, max(X1, X2))])
            assert S1 == S2
            for i, (ro1, ro2) in enumerate(index):
                o1 = (1 - scaf.orientation[X1]) if ro1 == 1 else scaf.orientation[X1]
                o2 = (1 - scaf.orientation[X2]) if ro2 == 1 else scaf.orientation[X2]
                # get_distsはp1.x < p2.xを仮定しているので
                if X1 < X2:
                    distances = get_dists(p1, p2, o1, o2, L1, L2, gap)
                else:
                    distances = get_dists(p2, p1, o2, o1, L2, L1, gap)
                p = estimator(distances).sum()
                # probは対称行列にならないとおかしい
                probs[S1][2*X1 + ro1, 2*X2 + ro2] = p
                probs[S1][2*X2 + ro2, 2*X1 + ro1] = p
    print('finished!')
    return probs

def get_dists(x, y, orientation1, orientation2, L1, L2, gap):
    # assume x, y are numpy.ndarray
    if orientation1 == 0 and orientation2 == 0:
        return L1-x+y+gap
    elif orientation1 == 1 and orientation2 == 0:
        return x+y+gap
    elif orientation1 == 0 and orientation2 == 1:
        return L1-x+L2-y+gap
    elif orientation1 == 1 and orientation2 == 1:
        return L2-y+x+gap

def infer_from_contig2(df, contigs, contig_id, K=100000, K0=3000):
    # generate global KDE estimation
    C = df[(df['X1']==contig_id) & (df['X2']==contig_id)]
    inter = np.abs(C['P1'].values - C['P2'].values)
    kde = KernelDensity(kernel='gaussian', bandwidth=200).fit(inter.reshape(-1, 1))
    f = lambda x: kde.score_samples(x.reshape(-1, 1))

    # distant
    x1 = np.logspace(np.log10(K0), np.log10(K), 500)
    p = lambda x, a, b: a + b * np.log(x)
    param1, cov = curve_fit(p, x1, f(x1))

    # proximal
    degree = 30
    x0 = np.logspace(0, np.log10(K0), 500)
    param0 = np.polyfit(x0, f(x0), degree)

    P = (lambda x: np.where( \
            x < K0, \
            np.poly1d(param0)(x), \
            np.where(x < K, param1[0] + param1[1] * np.log(x), param1[0] + param1[1] * np.log(K)) \
            ))

    # P = (lambda x: np.where( \
    #         x < K0, \
    #         param1[0] + param1[1] * np.log(K0), \
    #         np.where(x < K, param1[0] + param1[1] * np.log(x), param1[0] + param1[1] * np.log(K)) \
    #         ))

    return P, f

def infer_from_longest_contig2(df, contigs, K, K0):
    lengths = np.array(contigs.lengths)
    max_contig_id = np.argmax(lengths)
    return infer_from_contig2(df, contigs, max_contig_id, K=K, K0=K0)

def infer_from_contig(df, contigs, contig_name, maxlength):
    contig_id = contigs.get_id(contig_name)
    C = df[(df['X1']==contig_id) & (df['X2']==contig_id)]
    inter = np.abs(C['P1'].values - C['P2'].values)
    print('# of contacts: {}, contig_name: {}, contig_id: {}'.format(len(inter), contig_name, contig_id))
    f, raw = get_kde_polyfit_estimator(inter, \
                                       N=30000, bandwidth=200, \
                                       maxlength=maxlength, \
                                       points=500, degree=50)
    return f, raw

def infer_from_longest_contig(df, contigs, remove_repetitive=False, maxlength=150000, image_filename=None):
    """
    最長のcontigから推定器を作る
    """
    lengths = np.array(contigs.lengths)
    max_contig_id = np.argmax(lengths)
    print('longest contig to use: id={} length={}bp name={}'.format(max_contig_id, lengths[max_contig_id], contigs.get_name(max_contig_id)))
    if remove_repetitive:
        C = df[(df['X1']==max_contig_id) & (df['X2']==max_contig_id) & (df['U1']==True) & (df['U2']==True)]
    else:
        C = df[(df['X1']==max_contig_id) & (df['X2']==max_contig_id)]
    inter = np.abs(C['P1'].values - C['P2'].values)
    print('# of contacts:', len(inter))
    f, raw = get_kde_polyfit_estimator(inter, \
                                       N=30000, bandwidth=200, \
                                       maxlength=maxlength, \
                                       points=500, degree=50)
    # estimator_benchmark(inter, raw, f, maxlength=maxlength+50000, output=image_filename)
    return f, raw

def get_kde_estimator(samples, N=10000, bandwidth=200, debug=False):
    """sample(サンプル点のリスト)をKDEでノンパラメトリックにフィッティングして、その識別器を返す"""
    np.random.seed(0)
    downsampled = np.random.choice(samples, N)
    X = downsampled.reshape(-1, 1)
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(X)
    if debug:
        print('debug')
    return (lambda x: kde.score_samples(x.reshape(-1, 1)))

def get_kde_polyfit_estimator(samples, N=100000, bandwidth=200, maxlength=150000, points=500, degree=50):
    """多項式近似したバージョンを返すやつ 一応両方かえす"""
    f = get_kde_estimator(samples, N, bandwidth)
    x = np.linspace(1, maxlength, points)
    z = np.polyfit(x, f(x), degree)
    return (lambda x: np.where(x<=maxlength, np.poly1d(z)(x), np.poly1d(z)(maxlength))), f

def estimator_benchmark(inter, estimator, f, maxlength, output=None):
    """estimator: kde, f: polyfitted"""
    x1 = np.linspace(1, maxlength, 500)
    x2 = np.linspace(1, 20000, 200)
    plt.figure(figsize=(8, 8))
    #plt.subplot(3, 1, 1)
    #plt.xlabel('basepair distance between Hi-C contact')
    #plt.ylabel('frequency')
    #_ = plt.hist(inter, bins=200)

    #plt.subplot(3, 1, 2)
    plt.xlabel('basepair distance between Hi-C contact')
    plt.ylabel('log probability')
    #plt.xscale('log')
    plt.plot(x1, estimator(x1), label='KDE estimation')
    plt.plot(x1, f(x1), label='polynomial fit')
    plt.legend()
    plt.axvspan(1, 20000, facecolor='g', alpha=0.2)

    #plt.subplot(3, 1, 3)
    #plt.xlabel('basepair distance between Hi-C contact')
    #plt.ylabel('log probability')
    #plt.plot(x2, estimator(x2))
    #plt.plot(x2, f(x2))
    plt.legend()
    plt.axvspan(1, 20000, facecolor='g', alpha=0.2)
    # いい感じに調整
    plt.tight_layout()
    if output:
        # save
        plt.savefig(output)
    else:
        plt.show()
