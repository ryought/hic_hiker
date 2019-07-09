import numpy as np
import matplotlib.pyplot as plt

# matrixを作る系の関数
# contig同士の接触をみるやつ
def get_contact_matrix(hic):
    mat = np.zeros((hic.size, hic.size))
    for read_name, maps in hic.r.items():
        m1 = [row[0] for row in maps['r1']]
        m2 = [row[0] for row in maps['r2']]
        if len(m1) > 0 and len(m2) > 0:
            # unique map
            SR1 = set(m1)
            SR2 = set(m2)
            if len(SR1) == 1 and len(SR2) == 1:
                # unique map, so add up it
                # Case: repetitive region on same contig 1:[0] 2:[1,1,1,1,1,1] is almost surely (0,1)
                x, y = m1[0], m2[0]
                mat[min(x,y), max(x,y)] += 1
            # multiple hit, but origin was determined
            elif (len(SR1)==2 and len(SR2)==1 and list(SR2)[0] in SR1) or (len(SR1)==1 and len(SR2)==2 and list(SR1)[0] in SR2):
                # and case that joining point is inside of R1 or R2 read. (1:[0] 2:[0,1]) is maybe contact between 0 and 1
                x, y = list(SR1 | SR2)
                mat[min(x,y), max(x,y)] += 1
    for x in range(hic.size):
        for y in range(x+1, hic.size):
            mat[y,x] = mat[x,y]
    return mat

# contig内の接触について、あるbinsizeでcontact matrixを作るやつ
def get_all_contact_map_in_contig(hic, bin_size=1000):
    mat = []
    for i in range(len(hic.lengths)):
        length = hic.lengths[i]
        mat.append(np.zeros((length//bin_size + 1, length//bin_size + 1)))
    for read_name, maps in hic.r.items():
        if len(maps['r1']) == 1 and len(maps['r2']) == 1:
            if maps['r1'][0][0] == maps['r2'][0][0]:
                chrid = maps['r1'][0][0]
                pos_x, pos_y = maps['r1'][0][1], maps['r2'][0][1]
                i, j = pos_x // bin_size, pos_y // bin_size
                mat[chrid][i,j] += 1
                if i != j:
                    mat[chrid][j,i] += 1
    return mat

# orient matrixを作る用の関数
# これは単純に半分の長さで割っているので、怪しい。
def get_orient_matrix(hic):
    mat = np.zeros((hic.size*2, hic.size*2))
    for read_name, maps in hic.r.items():
        if len(maps['r1']) == 1 and len(maps['r2']) == 1:
            x, pos_x, y, pos_y = maps['r1'][0][0], maps['r1'][0][1], maps['r2'][0][0], maps['r2'][0][1]
            len_x, len_y = hic.lengths[x], hic.lengths[y]
            i = x*2 if pos_x < len_x/2 else x*2+1
            j = y*2 if pos_y < len_y/2 else y*2+1
            mat[i,j] += 1
            if i!=j:
                mat[j,i] += 1
    return mat

# weighted-center matrix用の関数
# 下のrandomじゃなくしたバージョン
# pos_list = get_contact_pos_list()
# weighted_centers = get_weighted_center_new(pos_list)
# mat = get_orient_matrix_weighted_center(r, size, lengths, weighted_centers)
# 接触した場所のリスト。後でソートする
def get_contact_pos_list(hic):
    pos_list = [[] for i in range(len(hic.ids))]
    for readname, record in hic.r.items():
        if len(record['r1']) == 1 and len(record['r2']) == 1:
            r1, r2 = record['r1'][0], record['r2'][0]
            c1, pos1 = r1[0], r1[1]
            c2, pos2 = r2[0], r2[1]
            pos_list[c1].append(pos1)
            pos_list[c2].append(pos2)
    return pos_list
# サンプルをもらって、その中心を返す関数
def calc_center(array):
    array.sort()
    N = len(array)
    if N == 0:
        return 0
    else:
        #print(N, array[N//2-2:N//2+2])
        #print(len([x for x in array if x == array[N//2]]))
        return array[N//2]
def get_weighted_center_new(pos_list):
    centers = []
    for i in range(len(pos_list)):
        pos_list[i].sort()
        centers.append(calc_center(pos_list[i]))
        if i >= 10:
            break
    return centers
def get_orient_matrix_weighted_center(hic, centers):
    mat = np.zeros((hic.size*2, hic.size*2))
    for read_name, maps in hic.r.items():
        if len(maps['r1']) == 1 and len(maps['r2']) == 1:
            x, pos_x, y, pos_y = maps['r1'][0][0], maps['r1'][0][1], maps['r2'][0][0], maps['r2'][0][1]
            len_x, len_y = hic.lengths[x], hic.lengths[y]
            i = x*2 if pos_x < centers[x] else x*2+1
            j = y*2 if pos_y < centers[y] else y*2+1
            mat[i,j] += 1
            if i!=j:
                mat[j,i] += 1
    return mat

# random weighted-center matrix用の関数
# うまくいっているやつ
def get_orient_matrix_random_weighted_center(hic):
    """あるcontigから出ているコンタクトreadを、そのposでソートする。中央値を決めてそのreadを半分(H(head) or T(tail))に分けたい。
    まずcontig iに落ちるreadを全部取ってきて、各readに対してどちらに入るかを決める。
    その後に、contig同士のコンタクトを集計し、その際にHTどちらかに入れるかは、先に決めた割り振りに従うようにする。
    """
    import random
    contacts = [[] for _ in range(hic.size)]
    for read_name, maps in hic.r.items():
        if len(maps['r1']) == 1 and len(maps['r2']) == 1:
            x, pos_x, y, pos_y = maps['r1'][0][0], maps['r1'][0][1], maps['r2'][0][0], maps['r2'][0][1]
            # 開始場所でソート、同じ順位のところからはランダムに拾ってくる
            contacts[x].append((read_name, 'r1', pos_x))
            contacts[y].append((read_name, 'r2', pos_y))
    is_head = {}
    for i in range(hic.size):
        random.shuffle(contacts[i])
        contacts[i].sort(key=lambda x:x[2]) # 開始位置でソート
        L = len(contacts[i])
        for j, c in enumerate(contacts[i]):
            is_head[c[0] + '_' + c[1]] = True if j < L//2 else False
    # matrixを作る
    mat = np.zeros((hic.size*2, hic.size*2))
    for read_name, maps in hic.r.items():
        if len(maps['r1']) == 1 and len(maps['r2']) == 1:
            x, pos_x, y, pos_y = maps['r1'][0][0], maps['r1'][0][1], maps['r2'][0][0], maps['r2'][0][1]
            i = x*2 if is_head[read_name + '_' + 'r1'] else x*2+1
            j = y*2 if is_head[read_name + '_' + 'r2'] else y*2+1
            mat[i,j] += 1
            if i!=j:
                mat[j,i] += 1
    return contacts, is_head, mat

# 向きを決めるのに必要な系の関数
def vote(array, verbose=False):
    if verbose:
        print('true:', len(np.where(array == True)[0]), 'false:', len(np.where(array == False)[0]))
    return len(np.where(array == True)[0]) >= len(np.where(array == False)[0])
def is_repetitive(i):
    return len(np.where(rep_regions[i]==True)[0])/rep_regions[i].shape[0] > 0.1

def benchmark_by_order(m, size):
    no, yes = 0, 0
    wrong = []
    for i in range(10, size-10):
        PH = m[2*i-1][2*i]
        PT = m[2*i-1][2*i+1]
        NH = m[2*i+2][2*i]
        NT = m[2*i+2][2*i+1]
        
        # 頭と前全部 これはなぜかめっちゃ正答率が悪い。多分巨大なcoverageのあるところに吸い込まれてしまってるんだと思う(小さいやつを足しても仕方ない、という感じ)
        PAH = np.sum(m[2*i][:2*i])  # 頭と前全部
        PAT = np.sum(m[2*i+1][:2*i])  # ケツと前全部
        NAH = np.sum(m[2*i][2*i+2:]) # 頭と後全部
        NAT = np.sum(m[2*i+1][2*i+2:])  # ケツと後全部
        
        # 頭と、contig全体ずつ（H+T)の値で多数決
        PCH = m[2*i][:2*i].reshape([-1,2]).sum(axis=1)
        PCT = m[2*i+1][:2*i].reshape([-1,2]).sum(axis=1)
        NCH = m[2*i][2*i+2:].reshape([-1,2]).sum(axis=1)
        NCT = m[2*i+1][2*i+2:].reshape([-1,2]).sum(axis=1)
        #print(m[2*i][2*i+2:].shape, NCH.shape, NAH == NCH.sum())
        
        # segmentごと(HT分離して考えてる)
        PSH = m[2*i][:2*i]
        PST = m[2*i+1][:2*i]
        NSH = m[2*i][2*i+2:]
        NST = m[2*i+1][2*i+2:]
        # feature vector(bool)
        FV = np.hstack(((PSH >= PST), (NSH <= NST)))
        #print(PSH, PST, NSH, NST)
        
        # PSH>PSTとなるのが正解
        OK_PS = len(np.where(PSH > PST)[0])
        NG_PS = len(np.where(PSH < PST)[0])
        # NSH<NSTとなるのが正解
        OK_NS = len(np.where(NSH < NST)[0])
        NG_NS = len(np.where(NSH > NST)[0])
        
        FV2 = np.hstack((PSH[-10:] >= PST[-10:], NSH[:10] <= NST[:10]))
        # 前後10でカウント
        #if (not PH >= PT) or (not NT >= NH):
        #if not(vote( >= ) == True and vote( <= ) == True):
        #if (not PAH >= PAT) and (not NAH <= NAT):  # これ成績悪い
        #if not(vote(PCH >= PCT) and vote(NCH <= NCT)):
        #if not(vote(PSH >= PST) and vote(NSH <= NST)): 最初と最後だけ間違えているようだ
        #if not(vote(PSH >= PST) and vote(NSH <= NST)):
        #if not(vote(FV)):
        #if (OK_PS+OK_NS) > (NG_PS+NG_NS):
        #if len(np.where(FV2 == True)[0]) > 10:
        if PH>=PT and NT>=NH:
            yes += 1
        else:
            print(i, PH>=PT, NT>=NH, FV2, FV2.shape)
            #print(i, FV2, len(np.where(FV2 == True)[0]))    
            no += 1
            wrong.append(i)
            #print(i, PAH, PAT, NAH, NAT, vote(PSH >= PST), vote(NSH <= NST))
            #print(PSH)
            #print(PST)
            #print(NSH)
            #print(NST)
            #print(i,PH,PT,NT,NH, PH>=PT, NT>=NH)
            #print('x', vote(m[2*i][:2*i] >= m[2*i+1][:2*i]))
            #print(m[2*i][:2*i])
            #print(m[2*i+1][:2*i])
            #print('y', vote(m[2*i][2*i+2:] <= m[2*i+1][2*i+2:]))
            #print(m[2*i][2*i+2:])
            #print(m[2*i+1][2*i+2:])
    print(no, yes)
    
    plt.hist(np.array(wrong), bins=50)
    plt.show()

# contig内のcoverageを計算する関数
def get_coverage():
    bin_size = 100
    # local_contacts[i][j]  contig iのj番目のbinに対する接触
    coverage = [[0 for _ in range(ans70.lengths[i] // bin_size + 1)] for i in range(len(ans70.ids))]
    for readname, record in ans70.r.items():
        if len(record['r1']) == 1 and len(record['r2']) == 1:
            r1, r2 = record['r1'][0], record['r2'][0]
            c1, pos1 = r1[0], r1[1] // bin_size
            c2, pos2 = r2[0], r2[1] // bin_size
            coverage[c1][pos1] += 1
            coverage[c2][pos2] += 1
    return coverage


# 描画系
def draw_matrix(matrix):
    plt.figure(figsize=(10,10))
    plt.imshow(matrix)
    plt.clim((0,2))
    plt.colorbar()
    plt.show()
    
def plot_contact_matrix_with_region(mat, region):
    plt.figure(figsize=(5,5))
    plt.subplot(1,2,1)
    plt.imshow(mat, aspect='auto')
    plt.clim((0,5))
    plt.subplot(1,2,2)
    plt.imshow(region.reshape((-1, 1)), aspect='auto')
    plt.show()

# matrix manipulation
# 上三角行列を対称行列にする
def symmetrize_matrix(matrix):
    size = matrix.shape[0]
    for x in range(size):
        for y in range(x+1, size):
            matrix[y,x] = matrix[x,y]
    return matrix

# contig iだけ向きを入れ替えた行列を作る これでbenchmarkを走らせて、iが返ってきたら成功
def get_reverted_matrix(matrix, i):
    M = matrix.copy()
    n = matrix.shape[0]
    indices = np.arange(n)
    indices[2*i] = 2*i+1
    indices[2*i+1] = 2*i
    return M[indices][:,indices]