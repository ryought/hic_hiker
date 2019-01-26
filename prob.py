import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
from tqdm import tqdm as tqdm
import pysam
import argparse

"""
確率を求める関係のやつ
"""

def infer_from_longest_contig(samr1_filename, samr2_filename):
    """
    実データで使う用
    一番長いコンティグ(or thresholdを設ける)のgapを使う
    ヒストグラム、fitted polylineを描画する機能をつける
    個数が不足していたら伸ばすなど
    TODO
    """
    sam_r1 = pysam.AlignmentFile(samr1_filename, 'r')
    sam_r2 = pysam.AlignmentFile(samr2_filename, 'r')
    size = sam_r1.nreferences

def get_intercontig_contacts_fast(samr1_filename, samr2_filename):
    """
    multiple hitを除いてないからどうなることやらバージョン
    """
    sam_r1 = pysam.AlignmentFile(samr1_filename, 'r')
    sam_r2 = pysam.AlignmentFile(samr2_filename, 'r')
    size = sam_r1.nreferences
    inter = []
    N = 0
    for r1, r2 in tqdm(zip(sam_r1, sam_r2)):
        # 同時に読み込む
        while r1.is_secondary or r1.is_supplementary:
            r1 = next(sam_r1)
        while r2.is_secondary or r2.is_supplementary:
            r2 = next(sam_r2)
        if r1.query_name != r2.query_name:
            print('assertion failed')
            print(r1.query_name, r2.query_name)
            break

        if (not r1.is_unmapped) and (not r2.is_unmapped):
            if r1.reference_id == r2.reference_id:
                p1, p2 = r1.reference_start, r2.reference_start
                inter.append(abs(p1-p2))
                N += 1
    print(N)
    return np.array(inter)

def get_intercontig_contacts(samr1_filename, samr2_filename):
    sam_ref_r1 = pysam.AlignmentFile(samr1_filename, 'r')
    sam_ref_r2 = pysam.AlignmentFile(samr2_filename, 'r')
    r = {}
    for read in sam_ref_r1:
        if not read.is_unmapped:
            if read.query_name not in r:
                r[read.query_name] = {'r1':[], 'r2':[]}
            r[read.query_name]['r1'].append((read.reference_id, read.reference_start, read.is_reverse))
    for read in sam_ref_r2:
        if not read.is_unmapped:
            if read.query_name in r:
                r[read.query_name]['r2'].append((read.reference_id, read.reference_start, read.is_reverse))
    # get interchromosomal contact
    inter = [ abs(c['r1'][0][1] - c['r2'][0][1]) for key, c in r.items() \
             if len(c['r1'])==1 and len(c['r2'])==1 and c['r1'][0][0] == c['r2'][0][0]]
    return np.array(inter)

def get_swapped_prob(prob, indexes=None):
    new_prob = prob.copy()
    if indexes:
        for i in indexes:
            new_prob[[2*i, 2*i+1]] = new_prob[[2*i+1, 2*i]]
            new_prob[:, [2*i, 2*i+1]] = new_prob[:, [2*i+1, 2*i]]
    return new_prob

def get_kde_estimator(samples, N=10000, bandwidth=100, debug=False):
    """sample(サンプル点のリスト)をKDEでノンパラメトリックにフィッティングして、その識別器を返す"""
    downsampled = np.random.choice(samples, N)
    X = downsampled.reshape(-1, 1)
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(X)
    if debug:
        print('debug')
    return (lambda x: kde.score_samples(x.reshape(-1, 1)))

def get_kde_polyfit_estimator(samples, N=10000, bandwidth=100, maxlength=40000, points=200, degree=20):
    """多項式近似したバージョンを返すやつ 一応両方かえす"""
    f = get_kde_estimator(samples, N, bandwidth)
    x = np.linspace(1, maxlength, points)
    z = np.polyfit(x, f(x), degree)
    return np.poly1d(z), f   # == (lambda x: np.poly1d(z)(x))

def get_prob_mixed_numpy(contacts, size, near_estimator, distant_estimator):
    """近いところと遠いところで違う関数を使うやつ"""
    prob = np.zeros((size*2, size*2))
    for d1, d2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        print(d1, d2)
        for i in tqdm(range(size)):
            for j in range(i+1, size):
                if j - i == 1:
                    gap = 0
                    c = contacts[i][j]
                    d = get_dists_numpy(c[0], c[1], d1, d2, 20000, 20000, gap)
                    if len(d) > 0:
                        p = near_estimator(d).sum()
                        prob[i*2+d1, j*2+d2] = p
                else:
                    gap = 20000*(j-i-1)
                    c = contacts[i][j]
                    d = get_dists_numpy(c[0], c[1], d1, d2, 20000, 20000, gap)
                    if len(d) > 0:
                        p = distant_estimator(d).sum()
                        prob[i*2+d1, j*2+d2] = p
    return prob

def get_prob_pandas(df, size, lengths, estimator):
    """
    pandasを読み込むよ
    contactsよりもメモリ消費が少ないかも
    今の所一番いいやつ
    """
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
        C = df[start:end]
        P1 = C['P1'].values
        P2 = C['P2'].values
        # 長さ
        L1 = lengths[i]
        L2 = lengths[j]
        # gap: i+1,..j-1のlenの合計
        gap = sum([lengths[k] for k in range(i+1, j)])
        for d1, d2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
            d = get_dists_numpy(P1, P2, d1, d2, L1, L2, gap)
            p = estimator(d).sum()
            prob[i*2+d1, j*2+d2] = p
    return prob

def get_prob_numpy(contacts, size, estimator):
    """普通のやつ"""
    prob = np.zeros((size*2, size*2))
    for d1, d2 in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        print(d1, d2)
        for i in tqdm(range(size)):
            for j in range(i+1, size):
                if j - i == 1:
                    gap = 0
                else:
                    gap = 20000*(j-i-1)
                c = contacts[i][j]
                d = get_dists_numpy(c[0], c[1], d1, d2, 20000, 20000, gap)
                if len(d) > 0:
                    p = estimator(d).sum()
                    prob[i*2+d1, j*2+d2] = p
    return prob

def estimate_distribution(ans, N=100000, bandwidth=10, debug=False):
    # estimate Hi-C distance distribution
    inter = [ abs(c['r1'][0][1] - c['r2'][0][1]) \
            for key, c in ans.r.items() if len(c['r1'])==1 and len(c['r2'])==1 \
            and c['r1'][0][0] == c['r2'][0][0]]
    inter_downsampled = np.random.choice(inter, N)

    X = inter_downsampled.reshape(-1, 1)
    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(X)
    if debug:
        X_test = np.linspace(0, 2000, 50)[:, np.newaxis]
        plt.subplot(2,1,1)
        plt.plot(X_test, np.exp(kde.score_samples(X_test)))  # fitted
        plt.subplot(2,1,2)
        plt.hist(X, bins=50)  # samples
        plt.show()

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
    psr.add_argument('output', help='output npy filename')
    args = psr.parse_args()
    
    if args.ref_r1 and args.ref_r2:
        inter = get_intercontig_contacts(args.ref_r1, args.ref_r2)
        f, raw = get_kde_polyfit_estimator(inter, N=100000, bandwidth=200, maxlength=150000, points=500, degree=50)
        
    if args.pandas:
        import feather
        df = feather.read_dataframe(args.pandas)
        size = max(df['X1'].max(), df['X2'].max()) + 1
        sam = pysam.AlignmentFile(args.sam, 'r')
        lengths = sam.lengths
        #print(lengths, size)
        prob = get_prob_pandas(df, size, lengths, f)
        np.save(args.output, prob)
        