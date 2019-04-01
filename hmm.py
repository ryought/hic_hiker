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


"""
TODO
- orderingの情報を入れられるようにする
多分emission probabilityのところに組み込めば良いだけだと思う
"""

def create_transition_matrix(k=2):
    n = 2**k
    A = np.zeros((n+1, n+1))
    for i in range(n):
        x = i % (n//2)
        A[i, 2*x] = 0.5
        A[i, 2*x+1] = 0.5
    A[n, :n] = 1/n
    return A

def get_layout(k, j):
    # ステートidからレイアウトを出力
    # for example, 2 in k=3 -> (1, 1, 0)
    layout = [0 for _ in range(k)]
    x = j
    for i in range(k-1, -1, -1):
        layout[i] = x % 2
        x = x // 2
    return layout

def get_emission_probability(i, j, prob, k=2):
    # 時間iのステートjでリードが出てくる確率
    # i,i+1,i+2のリードのペアのlog確率の和を出力すれば良い
    if i == 0:
        return -np.inf
    if j == 2**k:  # initial state will emit nothing
        return -np.inf
    if i == 1:
        # 最初は全部を出力する必要がある
        p = 0
        layout = get_layout(k, j)
        # {1,...,k}の全部列挙
        for x in range(k):
            for y in range(x+1, k):
                p += prob[ (x)*2 + layout[x], (y)*2 + layout[y] ]
        return p
    else:
        # 2以上だったら、 (i-1, i+) ()
        p = 0
        layout = get_layout(k, j)
        for x in range(k-1):
            # 最初のcontigをlとすると l, ..., l+k-1
            # 足し合わせたい組は(l, l+k-1), (l+1, l+k-1)...
            # 時間iの時先頭のコンティグはid i-1
            p += prob[ (i-1+x)*2 + layout[x], (i-1+k-1)*2 + layout[k-1] ]
        return p
    
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

def emission_prob_test(i, j): # iとi-1, i-2, ..., i-kの接触
    # k[i] に入っているコンティグの間の接触
    # 向きの情報はjに従う
    #return np.log(0.5)
    P = [
        [0.04, 0.02, 0.8, 0.01, 0.01, 0.01, 0.01, 0.1], # 8
        [0.01, 0.01, 0.01, 0.01, 0.93, 0.01, 0.01, 0.01], # 8
        [0.1, 0.7, 0.1, 0.1], # 4
        [0.1, 0.1, 0.1, 0.7] 
    ]
    return np.log(P[i][j])

def emission_prob(i, j, prob, k_list):
    if i == 0:
        # 最初は全部を出力する必要がある
        p = 0
        k = k_list[i][1]+1
        layout = get_layout(k, j)
        # {1,...,k}の全部列挙
        for x in range(k):
            for y in range(x+1, k):
                #print('contact', x, y, layout[x], layout[y])
                p += prob[ (x)*2 + layout[x], (y)*2 + layout[y] ]
        return p
    else:
        p = 0
        k = k_list[i][1]+1
        layout = get_layout(k, j)
        end = k_list[i][0]
        for x in range(1, k):
            #print('contact', (end-x), (end),layout[k-1-x], layout[k-1])
            p += prob[ (end-x)*2 + layout[k-1-x], (end)*2 + layout[k-1] ]
        return p

def transition_prob(i, j, kj):
    # t-1のステートiから、tのステートj(kj個を同時に考えている、つまりステート数は2**kj個)
    # jの1つシフトしたもの(kn-1桁ある) == iの下位kn-1桁
    # なら良い
    #print(bin(j>>1), bin(i), bin(2**(kj-1) -1), bin(i & (2**kj - 1)))
    if (j >> 1) == (i & (2**(kj-1) - 1)):
        return -np.log(2)
    else:
        return -np.inf
    
def run_hmm_adaptive_all(prob, K, n, lengths):
    k_list = get_ranges(lengths, K)
    return run_hmm_adaptive(prob, k_list)
    
def run_hmm_adaptive(prob, k):
    """
    prob:
    -> k: 各状態
    """
    #print(k)
    states = [np.zeros(2**(k[i][1]+1)) for i in range(len(k))]  # states[i][j]  i番目のステートでの向きjのvterbi確率を表す
    origin = [[-1 for _ in range(2**(k[i][1]+1))] for i in range(len(k))]  # origin[i][j]  state[i][j]はどこからきたか？
    Lstates = [k[i][1]+1 for i in range(len(k))]
    print(Lstates)
    
    # 初期化
    # 最初は確率1のinitial stateから、等確率で分配される
    p0 = - np.log(len(states[0]))
    print(np.exp(p0), p0)
    for j in range(len(states[0])):
        states[0][j] = p0 + emission_prob(0, j, prob, k)
    print('init', np.exp(states[0]))
    from tqdm import tqdm_notebook as tqdm
    for i in tqdm(range(1, len(states))):
        # time iに移行
        p = np.zeros(len(states[i-1]))
        for j in range(len(states[i])):
            #print(i, j)
            # t=iの時に state jになる確率を計算
            # その状態に入る最高値を求める
            for h in range(len(states[i-1])):
                # t=i-1の時のhのステートから、t=iの時のjのステートにくる確率
                # transition probability
                p[h] = states[i-1][h] + transition_prob(h, j, Lstates[i])
            Mp, Mpi = np.max(p), np.argmax(p)
            states[i][j] = Mp + emission_prob(i, j, prob, k)
            origin[i][j] = Mpi
    #print(states, origin)
    #print(states)
    path = [ np.argmax(states[-1]) ]
    for i in range(len(origin)-1, 0, -1):
        path.insert(0, origin[i][ path[0] ])
    print(path)
    
    # stateをほぐして、各contigの向きを出力
    # orientation
    orientation = [-1 for _ in range(k[-1][0] + 1)]  # contig数個の長さがあるはず
    for i in reversed(range(len(k))):
        # path[i]の一番右の桁を取ってくる (& 1 する)
        orientation[k[i][0]] = path[i] & 1
    for i in range(k[0][1]):
        # path[0]にはk[0][1]個分の情報が入っているので、全部の桁を取ってくる
        orientation[i] = path[0] & (2**(k[0][1] - i))
    
    return orientation, None, states, path, origin

def run_hmm(prob, k=4, n=5000, debug=False, short=False):
    """
    k: number of contig in one state
    n: number of contig
    ns: number of state in each time  
    """
    ns = 2**k + 1
    state = np.zeros((n-k+2, ns))  # state[i]: t=iにおける状態
    pointer = np.zeros((n-k+2, ns))
    state[0] = -np.inf
    state[0][-1] = 0 # 初期状態、一番後ろのやつが始状態
    
    # A[i][j] = ステートiからjに移動する確率(log)
    A = np.log(create_transition_matrix(k))
    
    # デバッグ用のデータ
    if debug:
        n = 4
        prob = np.zeros((2*n, 2*n))
        prob[0 + 0][2 + 0] = 0.05
        prob[0 + 0][2 + 1] = 0.75
        prob[0 + 1][2 + 0] = 0.1
        prob[0 + 1][2 + 1] = 0.1
        prob[2 + 0][4 + 0] = 0.1
        prob[2 + 0][4 + 1] = 0.1
        prob[2 + 1][4 + 0] = 0.7
        prob[2 + 1][4 + 1] = 0.1
        prob[4 + 0][6 + 0] = 0.8
        prob[4 + 0][6 + 1] = 0.05
        prob[4 + 1][6 + 0] = 0.05
        prob[4 + 1][6 + 1] = 0.1
        prob = np.log(prob)

    for i in range(1, n-k+2):
        # transition to time i
        # jというステートになる確率を計算する
        for j in range(ns):
            # 今までのステートでの最高値を求める
            M = state[i-1] + A[:,j]
            p = np.argmax(M)
            m = np.max(M)
            state[i][j] = get_emission_probability(i, j, prob, k) + m
            pointer[i][j] = p  # come from state p in t-1

    # 最適列の復元
    path = [ np.argmax(state[-1]) ]
    for i in range(n-k+2-1, 0, -1):
        path.insert(0, int(pointer[i][ path[0] ]))
    #print('path', path)

    # stateをほぐして、各contigの向きを出力
    # orientationの長さはnになるはず
    orientation = []
    L = len(path)
    for i in range(1, L):
        orientation.append(path[i] >> (k-1))  # 2進数にした時の一番左の位が一番若いstate
    for j in range(k-2, -1, -1): # 最後のステートにしかない情報は別に取り出す
        #print(j)
        orientation.append(path[i] >> (k-1) & 1)
    #print(k, n, sum(orientation))
    
    if short:
        # 各状態の計算はせずに、向きだけ必要な時
        return orientation, None, state, path
    
    
    # 各状態での確率の算出
    # 前向きアルゴリズム
    F = np.zeros((n-k+2, ns))  # F[i][j] = F_j(i) に対応
    F[0] = -np.inf
    F[0][-1] = 0
    #print(A)
    for i in range(1, n-k+2):
        for j in range(ns):
            tmp = F[i-1] + A[:,j]
            #print(i, j, np.exp(F[i-1] + A[:,j]), np.exp(F[i-1]), np.exp(A[:,j]))
            tmp2 = tmp[np.where(tmp != -np.inf)]
            if len(tmp2) > 0:
                S = tmp2[0]
                for l in range(1, len(tmp2)):
                    S = log_add(S, tmp2[l])
            else:
                S = -np.inf
            #print(i, j, np.exp(S))
            F[i][j] = get_emission_probability(i, j, prob, k) + S
            #print(np.exp(F[i][j]))
    
    # 後ろ向きアルゴリズム
    B = np.zeros((n-k+2, ns)) 
    B[-1] = 0
    for i in range(n-k, -1, -1):
        #print('b', i)
        for j in range(ns):
            S = -np.inf
            trans = A[j]
            tmp = B[i+1] + trans
            #trans[-1] = -np.inf
            #print('A', np.exp(trans))
            #print(i, j, np.exp(tmp))
            for l in range(ns):
                e = get_emission_probability(i+1, l, prob, k)
                #print(l, np.exp(e))
                #print(i, j, l, np.exp(tmp[l] + get_emission_probability(i+1, l, prob, k)))
                S = log_add(S, tmp[l] + e)
            B[i][j] = S
    
    P = F+B  # P[t][j] = 時刻tにjのステートにいる確率
    
    if debug:
        print('state', np.exp(state))
        print('for', np.exp(F))
        print('back', np.exp(B))
        print(np.exp(A))
        print(np.exp(P))
    
    return orientation, P, state, path

def log_sum(L):
    S = L[0]
    for i in range(1, len(L)):
        S = log_add(S, L[i])
    return S

def log_add(lA, lB):
    """
    suppose lA = log(A), lB = log(B) then return log(A+B)
    """
    if lA == -np.inf and lB == -np.inf:
        return -np.inf
    elif lA == -np.inf:
        return lB
    elif lB == -np.inf:
        return lA
    else:
        return max(lA, lB) + np.log(np.exp(lA - max(lA, lB)) + np.exp(lB - max(lA, lB)))

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