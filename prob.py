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

def infer_from_longest_contig(df, contigs, remove_repetitive=False, max_length=150000, image_filename=None):
    """
    最長のcontigから推定器を作る
    """
    lengths = np.array(contigs.lengths)
    max_contig_id = np.argmax(lengths)
    print('use longest contig', max_contig_id, lengths[max_contig_id], contigs.get_name(max_contig_id))
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
    # == (lambda x: np.poly1d(z)(x))
    #return (lambda x: np.poly1d(z)), f
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

from layout import Layout
from contigs import Contigs
def get_prob(df, contigs: Contigs, layout: Layout, estimator, max_k, remove_repetitive):
    # probsは各scaffoldに対応するprob matrixを入れるやつ
    probs = []
    for scaf in layout.scaffolds:
        size = scaf.N
        probs.append(np.zeros((size * 2, size * 2)))

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
                distances = get_dists_numpy(P1, P2, o1, o2, L1, L2, gap)
                result[i] = estimator(distances).sum()
        return pd.Series(result, index=index)

    # use DataFrame.groupby from pandas
    # https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
    prob_df = df.groupby(['X1', 'X2']).apply(calc_probs)

    for (x1, x2), p in prob_df.iterrows():
        # 

def get_prob_pandas(df, lengths, estimator, ordering=None, remove_repetitive=False, max_k=None):
    """
    pandasを読み込むよ
    contactsよりもメモリ消費が少ないかも
    今の所一番いいやつ
    """
    # マトリックスのサイズ
    size = max(df['X1'].max(), df['X2'].max()) + 1
    
    if remove_repetitive:
        print('removing repetitive hit')
        df = df[(df['U1'] == True) & (df['U2'] == True)]
    
    # orderingの変換
    # ordering[i]  = (i番目にいるcontigのid)
    # order_map[i] = (contig id iの位置)
    if ordering:
        order_map = [-1 for _ in range(len(ordering))]
        for x in range(len(ordering)):
            order_map[ordering[x]] = x
        
    # 値の切り替わりの位置を探す
    df.sort_values(by=['X1', 'X2'], inplace=True)
    x1, x2 = df['X1'].values, df['X2'].values
    start_pos = []
    L = len(x1)
    a, b = x1[0], x2[0]
    start_pos.append((x1[0], x2[0], 0))
    for i in range(1, L):
        if x1[i] == a and x2[i] == b:
            continue
        else:
            start_pos.append((x1[i], x2[i], i))
            a, b = x1[i], x2[i]
    # 確率の計算
    prob = np.zeros((size*2, size*2))
    #prob = 0  # 接触がないところは0とする np.infの方がいい？
    L = len(start_pos)
    for k in tqdm(range(L)):
        # i,jの接触
        i, j, start = start_pos[k]
        if k == L-1:
            end = L
        else:
            _, _, end = start_pos[k+1]
        if i == j:
            continue
        # limit k
        if max_k:
            if ordering:
                if abs(order_map[i]-order_map[j]) > max_k:
                    continue
            else:
                if abs(i-j) > max_k:
                    continue
        C = df[start:end]
        P1 = C['P1'].values
        P2 = C['P2'].values
        # 長さ
        L1 = lengths[i]
        L2 = lengths[j]
        if ordering:
            # 配置上の場所
            I = order_map[i]
            J = order_map[j]
            gap = sum([ lengths[ordering[k]] for k in range(min(I,J)+1, max(I,J)) ])
            
            if I > J:
                # I < Jに直す
                J, I = I, J
                L1, L2 = L2, L1
                P1, P2 = P2, P1
        else:
            # gap: i+1,..j-1のlenの合計
            gap = sum([ lengths[k] for k in range(i+1, j) ])
            
        for d1, d2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            d = get_dists_numpy(P1, P2, d1, d2, L1, L2, gap)
            p = estimator(d).sum()
            if ordering:
                prob[I*2+d1, J*2+d2] = p
                prob[J*2+d2, I*2+d1] = p
            else:
                prob[i*2+d1, j*2+d2] = p
                prob[j*2+d2, i*2+d1] = p
    return prob

def get_dists(contacts, direction1, direction2, L1, L2, gap):
    # get_dists([(1,100), ...], True, False)
    if direction1 == 0 and direction2 == 0:
        return [ L1-x+y+gap for x,y in contacts ]
    elif direction1 == 1 and direction2 == 0:
        return [ x+y+gap for x,y in contacts ]
    elif direction1 == 0 and direction2 == 1:
        return [ L1-x+L2-y+gap for x,y in contacts ]
    elif direction1 == 1 and direction2 == 1:
        return [ L2-y+x+gap for x,y in contacts ]

def get_dists_numpy(x, y, direction1, direction2, L1, L2, gap):
    # assume x, y are numpy.ndarray
    if direction1 == 0 and direction2 == 0:
        return L1-x+y+gap
    elif direction1 == 1 and direction2 == 0:
        return x+y+gap
    elif direction1 == 0 and direction2 == 1:
        return L1-x+L2-y+gap
    elif direction1 == 1 and direction2 == 1:
        return L2-y+x+gap

## ↓↓古いコード
def calc_prob_original(contacts, ans, kde):
    prob = np.zeros((ans.size*2, ans.size*2))
    for d1, d2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        print(d1, d2)
        for i in range(ans.size):
            for j in range(i+1, ans.size):
                d = get_dists(contacts[i][j], d1, d2, ans.lengths[i], ans.lengths[j])
                d = np.array(d)
                if len(d) > 0:
                    p = kde.score_samples(d.reshape(-1,1)).sum()
                    prob[i*2+d1, j*2+d2] = p
                else:
                    prob[i*2+d1, j*2+d2] = 0
    return prob

# Li = ans.lengths[i]
def calc_prob(i, j, contact, Li, Lj):
    ret = []
    for d1, d2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        d = get_dists(contact, d1, d2, Li, Lj)
        d = np.array(d)
        if len(d) > 0:
            p = kde.score_samples(d.reshape(-1,1)).sum()
        else:
            p = 0
        ret.append((i, j, d1, d2, p))
    # [(i, j, 0, 0, prob), (i, j, 0, 1, prob)]
    return ret

def calc_prob_parallel(contacts, ans, kde):
    from joblib import Parallel, delayed
    jobs = [delayed(calc_prob)(i, j, contacts[i][j], ans.lengths[i], ans.lengths[j]) \
        for i,j in itertools.combinations(range(ans.size), 2)]
    result = Parallel(n_jobs=-1, verbose=3)(jobs)
    print('finished')
    prob = np.zeros((ans.size*2, ans.size*2))
    for r in result:
        for row in r:
            i, j, d1, d2, p = row
            prob[2*i+d1, 2*j+d2] = p
    return prob
    
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
        
