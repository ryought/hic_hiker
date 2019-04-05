import numpy as np
import matplotlib.pyplot as plt
import numpy as np

# hmmの実装

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
    if j == 2**k:
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

def main(k=10, n=5000):
    ns = 2**k + 1
    state = np.zeros((n-k+2, ns))  # state[i] t=iにおける状態
    pointer = np.zeros((n-k+2, ns))
    state[0] = -np.inf
    state[0][-1] = 0 # 初期状態、一番後ろのやつが始状態
    print(state[0])
    
    # A[i][j] = ステートiからjに移動する確率
    A = create_transition_matrix(k)
    print(A)
    A = np.log(A)
    
    prob = np.load('prob_70x.npy')
    print(prob.shape)
    
    for i in range(1, n-k+2):
        for j in range(ns):
            # 今までのステートで最高ちを求める
            M = state[i-1] + A[:,j]
            p = np.argmax(M)
            m = np.max(M)
            state[i][j] = get_emission_probability(i, j, prob, k) + m
            pointer[i][j] = p

    # 最適列の復元
    path = [ np.argmax(state[-1]) ]
    for i in range(n-k+2-1, 0, -1):
        path.insert(0, int(pointer[i][ path[0] ]))
    print(state[-1])
    print('(result) path:', path, 'log prob:', np.max(state[-1]))

    ordering = []
    L = len(path)
    for i in range(1, L):
        ordering.append(path[i] >> (k-1))
        # if i == L-1:
            # for j in range(k-2, -1, -1):
                # ordering.append(path[i] >> j)
    print(path)
    print(ordering)
    print(sum(ordering))

main()
