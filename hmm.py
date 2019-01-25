import numpy as np
import argparse

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '.'))

#np.set_printoptions(formatter={'float': '{: 0.3e}'.format}, linewidth=200)
np.set_printoptions(suppress=True)

"""
HMM module
- viterbi algorithm for optimal path
- traceback
- forward and backward algorithm for probability of each state
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
    else:
        p = 0
        layout = get_layout(k, j)
        for x in range(k):
            for y in range(x+1, k):
                # contig i-1+x, contig i-1+y 間の、
                # それぞれ layout[x] layout[y] の向きとした時の確率を出す
                # x < yになってる
                p += prob[ (i-1+x)*2 + layout[x], (i-1+y)*2 + layout[y] ]
        return p
    
def show_benchmark(layout, pro, k, head=10):
    """debug用のbenchmark関数"""
    #P = np.exp(pro)
    P = pro
    N = 0
    print(layout)
    for i in range(len(layout)):
        N += 1
        if N == head:
            continue
        if layout[i] == 0:
            #continue
            pass
        
        # layout 表示部分
        print(i, layout[i], end=' ')
        
        for j in range(len(P[i])):
            if j == np.argmax(P[i]):
                print('\x1b[31m{:.7e}\x1b[0m'.format(P[i][j]), end=' ')
            else:
                print('{:.7e}'.format(P[i][j]), end=' ')
        
        print(end='\n')
        print(np.exp(P[i+1][0:2**(k-1)]).sum(), np.exp(P[i+1][2**(k-1):2**k]).sum())

def run_hmm(prob, k=4, n=5000, debug=False):
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
        print(j)
        orientation.append(path[i] >> (k-1) & 1)
    print(k, n, sum(orientation))
    
    
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
        print(np.exp(state))
        print(np.exp(F))
        print(np.exp(B))
        print(np.exp(A))
        print(np.exp(P))
    
    return orientation, P


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