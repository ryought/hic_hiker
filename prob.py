import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors.kde import KernelDensity
from tqdm import tqdm_notebook as tqdm
import pysam

"""
確率を求める関係のやつ
"""

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
    psr.add_argument('contacts_pkl', help='pickle filename of contacts')
    psr.add_argument('n', type=int, help='number of contigs. O(n)')
    psr.add_argument('--prob', help='npy file of probabilities. You can create this file by running load.py')
    psr.add_argument('--debug', help='run with test data', action='store_true')
    args = psr.parse_args()
    
    if args.debug:
        _ = run_hmm(None, k=args.k, n=args.n, debug=True)
    elif args.prob:
        print(args.prob, args.k, args.n)
        prob = np.load(args.prob)
        _ = run_hmm(prob, k=args.k, n=args.n)
    else:
        print('you should select either --prob prob.npy or --debug option to run hmm.')