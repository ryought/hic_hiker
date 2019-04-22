import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

np.set_printoptions(formatter={'float': '{: 0.10e}'.format}, linewidth=200)
#np.set_printoptions(suppress=True)

"""
HMM module
- viterbi algorithm for optimal path
- traceback
- forward and backward algorithm for probability of each state
"""

def get_layout(k, j):
    # generate layout of contigs in the state, given k(number of contigs in the state) and j(state id)
    # for example, j=2 in k=3 -> (1, 1, 0)
    assert 0 <= j < 2**k, 'state index out of range'
    layout = [0 for _ in range(k)]
    x = j
    for i in range(k-1, -1, -1):
        layout[i] = x % 2
        x = x // 2
    return layout


def emission_prob(i, j, prob, ks):
    # DONE
    """時刻iのstate jでの出力確率を返す"""
    p = 0.0
    k = ks[i]
    layout = get_layout(k, j)
    if i == 0:
        # 最初は全部を出力する必要がある。
        # contig0, 1, ..., k-1の間のprobを全部合計する
        for x in range(k):
            for y in range(x+1, k):
                p += prob[ x*2 + layout[x], y*2 + layout[y] ]
    else:
        # i != 0の場合
        # そのstateで考えているcontigたちは[end-(k-1), end-(k-2), ..., end]のk個。
        # 新たに追加される確率は(end-(k-1), end), (end-(k-2), end), ... ,(end-1, end)
        end = (ks[0] - 1) + i
        for x in range(1, k):
            p += prob[ (end-x)*2 + layout[k-1-x], end*2 + layout[k-1] ]
    return p

def transition_prob(j0, j1, k):
    # DONE
    # t-1のステートj0から、tのステートj1(k個を同時に考えている、つまりステート数は2**k個)
    # j1の1つシフトしたもの(k-1ビットある) == j0の下位k-1ビットなら良い
    # jの下nビットを取り出す: j & (2**n - 1)
    if (j1 >> 1) == (j0 & (2**(k-1) - 1)):
        return -np.log(2)
    else:
        return -np.inf

def run_hmm_adaptive(prob, ks):
    """
    prob -> ks -> path
    ks: 各時刻のstateで考えるcontigの上下のリスト
    ks :: [(contig_id, length)]
      ここでcontig_idはそのstateの最初のcontigのid、lengthは同時に考えるcontigの個数
    """
    assert ((ks[0] - 1) + len(ks)) * 2 == prob.shape[0] == prob.shape[1]
    N = len(ks)  # number of states
    states = [ np.zeros(2**ks[i])            for i in range(N) ]
    origin = [ [-1 for _ in range(2**ks[i])] for i in range(N) ]  # origin[i][j]  state[i][j]はどこからきたか？

    # 初期化 i = 0の時
    # 最初は確率1のinitial stateから、等確率で分配される。それぞれ\frac{1}{|S0|}の確率
    p0 = - np.log(len(states[0]))
    for j in range(len(states[0])):
        # 初期確率と、その時のemissionを掛けたものがviterbi確率
        states[0][j] = p0 + emission_prob(i=0, j=j, prob=prob, ks=ks)

    from tqdm import tqdm_notebook as tqdm
    for i in tqdm(range(1, N)):
        # time iに移行
        p = np.zeros(len(states[i-1]))
        for j in range(len(states[i])):
            #print(i, j)
            # t=iの時に state jになる確率を計算
            # その状態に入る最高値を求める
            for h in range(len(states[i-1])):
                # p[h] = (t=i-1の時のhのステートから、t=iの時のjのステートにくる確率)
                p[h] = states[i-1][h] + transition_prob(h, j, ks[i])
            Mp, Mpi = np.max(p), np.argmax(p)
            states[i][j] = Mp + emission_prob(i=i, j=j, prob=prob, ks=ks)
            origin[i][j] = Mpi

    # backtrace
    path = [ np.argmax(states[-1]) ]
    for i in range(len(origin) - 1, 0, -1):
        path.insert(0, origin[i][ path[0] ])
    assert len(path) == N

    # stateをほぐして、各contigの向きを出力
    # orientation
    Ncontigs = (ks[0] - 1) + N  # contig数
    orientation = [-1 for _ in range(Ncontigs)]  # contig数個の長さがあるはず
    for i in reversed(range(N)):
        # path[i]の一番右の桁を取ってくる (& 1 する)
        # state iは
        orientation[(ks[0] - 1) + i] = path[i] & 1
    for i in range(ks[0] - 1):
        # path[0]にはks[0]個分の情報が入っているので、全部の桁を取ってくる
        orientation[i] = (path[0] >> (ks[0]-i-1)) & 1

    return orientation, path, states, origin

def get_ks_constant_k(Ncontigs, k):
    return [ k for _ in range(Ncontigs - (k-1)) ]

def get_ks_adaptive(lengths, K):
    return ()

def show_benchmark(layout, pro, k, head=10):
    """debug用のbenchmark関数"""
    #P = np.exp(pro)
    P = pro
    N = 0
    locus_prob = []
    for i in range(len(layout)-k):
        N += 1
        if layout[i] == 0:
            #continue
            pass
        
        # layout 表示部分
        U = P[i+1][0:2**(k-1)]
        D = P[i+1][2**(k-1):2**k]
        sU = log_sum(U)
        sD = log_sum(D)
        sUD = log_add(sU, sD)
        
        locus_prob.append((i, layout[i], np.exp(sU-sUD), np.exp(sD-sUD)))
        
        if N < head:
            print('\x1b[31m', i, layout[i], '\x1b[0m', np.exp(sU-sUD), np.exp(sD-sUD))
        
    for step in range(k, 0, -1):
        N += 1
        i = len(layout)-step
        LP = P[len(layout)-k+1]
        # step-1個右にシフトして、1とandとった奴が0の時上、1の時下 (最下位ビットが0ならU)
        U, D = [], []
        for j in range(len(LP)-1):
            if (j >> (step-1) & 1) == 0:
                #rint(j)
                U.append(LP[j])
            else:
                D.append(LP[j])
        sU = log_sum(U)
        sD = log_sum(D)
        sUD = log_add(sU, sD)
        locus_prob.append((i, layout[i], np.exp(sU-sUD), np.exp(sD-sUD)))
        if N < head:
            print('\x1b[31m', i, layout[i], '\x1b[0m', np.exp(sU-sUD), np.exp(sD-sUD))
            
    return locus_prob

def get_ranges(lengths, K=20):
    """
    adaptiveなやつ用
    K: [(2, 2), (3, 2), (4, 1), (5, 1)]  各状態の開始index(0-origin)と、その隣を何個考えるか
    """
    k = []
    for i in range(1, len(lengths)):
        # iで終わる区間について
        #print(i)
        #print(lengths[:i])
        lst = np.cumsum(list(reversed(lengths[:i])))
        #print(lst)
        idxs = np.arange(i)[lst > K]
        #print(idxs)
        if len(idxs) == 0:
            continue
        else:
            k.append((i, idxs[0]+1))
    k[0] = (k[0][0], k[0][0])
    return k

if __name__ == '__main__':
    psr = argparse.ArgumentParser()
    psr.add_argument('k', type=int, help='number of neighbor contigs that consist of 1 state. O(2^k)')
    psr.add_argument('n', type=int, help='number of contigs. O(n)')
    psr.add_argument('--prob', help='npy file of probabilities. You can create this file by running load.py')
    psr.add_argument('--debug', help='run with test data', action='store_true')
    args = psr.parse_args()
    
    if args.debug:
        layout, P = run_hmm(None, k=args.k, n=args.n, debug=True)
        show_benchmark(layout, P, k=args.k)
        
    elif args.prob:
        print(args.prob, args.k, args.n)
        prob = np.load(args.prob)
        layout, P = run_hmm(prob, k=args.k, n=args.n)
        show_benchmark(layout, P, k=args.k)
    else:
        print('you should select either --prob prob.npy or --debug option to run hmm.')
