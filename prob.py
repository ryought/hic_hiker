import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
from tqdm import tqdm_notebook as tqdm
#from tqdm import tqdm as tqdm
import pysam
import argparse
import pandas as pd

"""
確率を求める関係のやつ
df, contigs
"""

from layout import Layout
from contigs import Contigs
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
    print('calcing version2')
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
    print('iter version2')
    for (I1, I2), idx in tqdm(groups):
        # print
        # d = df.iloc[idx]  # これ使うと遅い
        #   P1, P2 = d.P1.values, d.P2.values
        S1, X1 = layout.id2order(I1)
        S2, X2 = layout.id2order(I2)
        # contig place
        # result = [0, 0, 0, 0]
        if (S1 != S2) or (max_k and abs(X1 - X2) > max_k):
            # scaffoldが違うところに乗っていたり、
            # order上でk以上離れている確率は使わないので無視する
            pass
        else:
            # print(I1, I2)
            # scaffold
            scaf = layout.scaffolds[S1]  # = S2
            # position list
            p1, p2 = P1[idx], P2[idx]
            # length list
            L1, L2 = lengths[I1], lengths[I2]
            gap = sum([ lengths[scaf.order[x]] for x in range(min(X1, X2)+1, max(X1, X2))])
            assert S1 == S2
            for i, (o1, o2) in enumerate(index):
                distances = get_dists(p1, p2, o1, o2, L1, L2, gap)
                p = estimator(distances).sum()
                # probは対称行列にならないとおかしい
                probs[S1][2*X1 + o1, 2*X2 + o2] = p
                probs[S1][2*X2 + o2, 2*X1 + o1] = p
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
    estimator_benchmark(inter, raw, f, maxlength=maxlength+50000, output=image_filename)
    return f

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
    plt.subplot(3, 1, 1)
    plt.xlabel('basepair distance between Hi-C contact')
    plt.ylabel('frequency')
    _ = plt.hist(inter, bins=200)

    plt.subplot(3, 1, 2)
    plt.xlabel('basepair distance between Hi-C contact')
    plt.ylabel('log probability')
    #plt.xscale('log')
    plt.plot(x1, estimator(x1), label='KDE estimation')
    plt.plot(x1, f(x1), label='polynomial fit')
    plt.legend()
    plt.axvspan(1, 20000, facecolor='g', alpha=0.2)

    plt.subplot(3, 1, 3)
    plt.xlabel('basepair distance between Hi-C contact')
    plt.ylabel('log probability')
    plt.plot(x2, estimator(x2))
    plt.plot(x2, f(x2))
    plt.legend()
    plt.axvspan(1, 20000, facecolor='g', alpha=0.2)
    # いい感じに調整
    plt.tight_layout()
    if output:
        # save
        plt.savefig(output)
    else:
        plt.show()

if __name__ == '__main__':
    psr = argparse.ArgumentParser()
    psr.add_argument('--pickle', help='pickle filename of contacts')
    psr.add_argument('--pandas', help='pandas')
    psr.add_argument('--ref_r1', help='reference')
    psr.add_argument('--ref_r2', help='reference')
    psr.add_argument('--sam', help='contig sams for lengths')
    psr.add_argument('--max_length', type=int, help='max length')
    psr.add_argument('output', help='output npy filename')
    args = psr.parse_args()
    if args.ref_r1 and args.ref_r2:
        print('infer from reference contact')
        inter = get_intercontig_contacts(args.ref_r1, args.ref_r2)
        f, raw = get_kde_polyfit_estimator(inter, N=100000, bandwidth=200, maxlength=150000, points=500, degree=50)
    if args.pandas:
        import feather
        df = feather.read_dataframe(args.pandas)
        if not (args.ref_r1 and args.ref_r2):
            print('infer from longest contig')
            f = infer_from_longest_contig(df, args.sam, args.output+'.png', maxlength=args.max_length, remove_repetitive=True)
        sam = pysam.AlignmentFile(args.sam, 'r')
        lengths = sam.lengths
        prob = get_prob_pandas(df, lengths, f, remove_repetitive=True)
        np.save(args.output, prob)
