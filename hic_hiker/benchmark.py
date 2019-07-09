#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ベンチマークする
"""
from .layout import Scaffold, Layout
from .contigs import Contigs

# 評価結果
OK = 'ok'
ORI_ERR = 'orientation_error'
ORD_ERR = 'order_error'
NONE = 'none'

def determine_correct_orientation_or_not(contigs: Contigs, layout: Layout, ref_layout: Layout):
    """
    2つのlayoutを比較する関数
    特に片方は正解だとしている
    各contigが正解しているか？のデータに変換する
    """
    from collections import Counter
    result = []
    joineds = []

    for i, scaf in enumerate(layout.scaffolds):
        # まずこのscafがref_scafのどの辺に対応するかを調べる
        # ref_scaffoldの中のchr_id番目のscaffoldの、chr_direction(0,1)の向きに対応していると考える
        # chr_counter = Counter([ref_layout.id2order(contig_id)[0] for contig_id in scaf.order])
        # chr_id, _ = chr_counter.most_common()[0]

        # 下のやつだとKeyErrorが出てしまうので
        # joined = [ (contig_id, ori, ref_scaf, ref_ord, ref_ori, contig_id2) for contig_id, ori, (ref_scaf, ref_ord, contig_id2, ref_ori)
                # in zip(scaf.order, scaf.orientation, [(x, y, ref_layout.scaffolds[x].order[y], ref_layout.scaffolds[x].orientation[y]) for x,y in [ref_layout.id2order(i) for i in scaf.order]])]
        # joined作る
        joined = []
        for i in range(scaf.N):
            contig_id = scaf.order[i]
            ori = scaf.orientation[i]
            try:
                ref_scaf_id, ref_order = ref_layout.id2order(contig_id)
                ref_ori = ref_layout.scaffolds[ref_scaf_id].orientation[ref_order]
            except KeyError:
                ref_scaf_id, ref_order = -1, 0
                ref_ori = 0
            joined.append((contig_id, ori, ref_scaf_id, ref_order, ref_ori))
        assert len(joined) == len(scaf.order), 'joined length assertion'
        # 保存
        joineds.append(joined)

        # scaf_result作る
        scaf_result = []
        if scaf.N > 2:
            # 最初と最後は不明なのでnone
            scaf_result.append('none')
            for i in range(1, len(joined) - 1):

                # 評価
                is_on_same_scaffold = (joined[i-1][2] == joined[i][2] == joined[i+1][2])
                is_ascending        = (joined[i-1][3] <  joined[i][3] <  joined[i+1][3])
                is_descending       = (joined[i-1][3] >  joined[i][3] >  joined[i+1][3])

                # 向き
                ori = joined[i][1]
                ref_ori = joined[i][4]

                if is_on_same_scaffold and is_descending:
                    # 反対向きのものが正解
                    if ori != ref_ori:
                        scaf_result.append(OK)
                    else:
                        scaf_result.append(ORI_ERR)
                elif is_on_same_scaffold and is_ascending:
                    # 同じ順序のものが正解
                    if ori == ref_ori:
                        scaf_result.append(OK)
                    else:
                        scaf_result.append(ORI_ERR)
                else:
                    scaf_result.append(ORD_ERR)
            # 最初と最後は不明なのでnone
            scaf_result.append(NONE)
        else:
            for i in range(len(joined)):
                scaf_result.append(NONE)
        assert len(scaf_result) == len(scaf.order), 'scaf test, {}, {}'.format(len(scaf_result), len(scaf.order))
        # 保存
        result.append(scaf_result)
    return result, joineds

if __name__ == '__main__':
    run_benchmark()


