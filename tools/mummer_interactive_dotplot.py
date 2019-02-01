#!/usr/bin/env python
# -*- coding: utf-8 -*-
import plotly.offline as offline
import plotly.graph_objs as go
import csv
import sys

def parse_multi_mums(filename):
    ret = []
    query_store = None
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for l in reader:
            if l[0] == '>':
                if l[-1] == 'Reverse':
                    # 逆向きの始まり
                    direction = '-'
                else:
                    # 新しいqueryの始まり
                    if query_store:
                        ret.append(query_store)
                    query_store = {
                        'query_name': l[1],
                        '+': [],
                        '-': []
                    }
                    direction = '+'
            else:
                # (ref_id, refでの開始位置, queryでの開始位置, match length)
                query_store[direction].append((l[0], int(l[1]), int(l[2]), int(l[3])))
        ret.append(query_store)
    return ret

# parse output of mummer4
def parse_mums(filename):
    ret = {
        '+': [],
        '-': []
    }
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for l in reader:
            if l[0] == '>':
                if l[-1] == 'Reverse':
                    direction = '-'
                else:
                    direction = '+'
            else:
                ret[direction].append((int(l[0]), int(l[1]), int(l[2])))
    return ret
# output: result['+'] = [(referenceでの位置, queryでの位置, match length  )]

def draw_map_plotly(mums, html=None):
    data = []
    for align in mums['+']:
        x = align[0]
        y = align[1]
        length = align[2]
        data.append(go.Scattergl(x=[x,x+length], y=[y,y+length], marker=dict(color='red')))
    for align in mums['-']:
        x = align[0]
        y = align[1]
        length = align[2]
        data.append(go.Scattergl(x=[x,x+length], y=[y,y-length], marker=dict(color='blue')))

    # output
    layout = go.Layout(
        showlegend=False, barmode='stack',
        xaxis=dict(title='reference', domain=[0, 1], showspikes=True, spikemode='across'),
        yaxis=dict(title='query', domain=[0, 1], showspikes=True, spikemode='across')
    )
    fig = go.Figure(data=data, layout=layout)
    if html:
        offline.plot(fig, filename=html, auto_open=False)
    else:
        offline.iplot(fig)

def draw_map_plotly_multi(mums, refs, querys, html=None):
    data = []
    refs_axes = ['x'+str(i+1) if i != 0 else 'x' for i in range(len(refs))]
    querys_axes = ['y'+str(i+1) if i != 0 else 'y' for i in range(len(querys))]
    for i in range(len(refs)):
        for j in range(len(querys)):
            for query_store in mums:
                if query_store['query_name'] == querys[j]:
                    for align in query_store['+']:
                        if align[0] == refs[i]:
                            x = align[1]
                            y = align[2]
                            length = align[3]
                            data.append(go.Scattergl(x=[x,x+length], y=[y,y+length], 
                                                     xaxis=refs_axes[i],
                                                     yaxis=querys_axes[j],
                                                     marker=dict(color='red')))
                    for align in query_store['-']:
                        if align[0] == refs[i]:
                            x = align[1]
                            y = align[2]
                            length = align[3]
                            data.append(go.Scattergl(x=[x,x+length], y=[y,y-length],
                                                     xaxis=refs_axes[i],
                                                     yaxis=querys_axes[j],
                                                     marker=dict(color='blue')))
    # output
    option = dict(showlegend=False, barmode='stack')
    # references
    pitch = 1 / len(refs)
    start, end = 0, pitch
    for i in range(len(refs)):
        option['xaxis' + str(i+1) if i != 0 else 'xaxis'] = dict(
            title=refs[i], domain=[start, end], showspikes=True, spikemode='across'
        )
        start += pitch
        end += pitch
    # querys
    pitch = 1 / len(querys)
    start, end = 0, pitch
    for j in range(len(querys)):
        option['yaxis' + str(j+1) if j != 0 else 'yaxis'] = dict(
            title=querys[j], domain=[start, end], showspikes=True, spikemode='across'
        )
        start += pitch
        end += pitch
    layout = go.Layout(option)
    fig = go.Figure(data=data, layout=layout)
    if html:
        offline.plot(fig, filename=html, auto_open=False)
    else:
        offline.iplot(fig)

def main_for_3ddna():
    mum = parse_multi_mums('analysis/mock/mock_chrII.final.mum')
    draw_map_plotly_multi(mum,
        # refs=['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'],
        refs=['chrII'],
        # querys=['supercontig', 'supercontig_original'],
        querys=['HiC_scaffold_1', 'HiC_scaffold_2'],
        # querys=['HiC_scaffold_2', 'HiC_scaffold_5', 'HiC_scaffold_3', 'HiC_scaffold_6', 'HiC_scaffold_4', 'HiC_scaffold_1'],
        html='analysis/html/mock_chrII_3ddna.html')

def main_super():
    mum = parse_multi_mums('analysis/mock/mock_chrII.3ddna.mum')
    draw_map_plotly_multi(mum,
        # refs=['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'],
        refs=['chrII'],
        querys=['supercontig'],
        # querys=['HiC_scaffold_1', 'HiC_scaffold_2'],
        # querys=['HiC_scaffold_2', 'HiC_scaffold_5', 'HiC_scaffold_3', 'HiC_scaffold_6', 'HiC_scaffold_4', 'HiC_scaffold_1'],
        html='analysis/html/mock_chrII_0.html')

def main():
    mum = parse_multi_mums('analysis/mock/mock_chrI_rev_tmp.mum')
    draw_map_plotly_multi(mum,
        # refs=['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'],
        refs=['chrI'],
        querys=['supercontig', 'supercontig_original'],
        # querys=['HiC_scaffold_1', 'HiC_scaffold_2'],
        # querys=['HiC_scaffold_2', 'HiC_scaffold_5', 'HiC_scaffold_3', 'HiC_scaffold_6', 'HiC_scaffold_4', 'HiC_scaffold_1'],
        html='analysis/html/mock_chrI_rev_combined.html')

def mainA():
    mum = parse_multi_mums('analysis/answer/answer70gap.mum')
    for r, q in zip(['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrX'], ['HiC_scaffold_1', 'HiC_scaffold_2', 'HiC_scaffold_3', 'HiC_scaffold_4', 'HiC_scaffold_5', 'HiC_scaffold_6']):
        draw_map_plotly_multi(mum, refs=[r], querys=[q],
            html='analysis/html/answer70_{}.html'.format(r))

def main2():
    if len(sys.argv) != 3:
        print('usage: python mummer_interactive_dotplot.py hoge.mum hoge.html')
        return -1
    filename = sys.argv[1]
    output_filename = sys.argv[2]
    mum = parse_mums(filename)
    draw_map_plotly(mum, html=output_filename)

if __name__ == '__main__':
    # main()
    # main_for_3ddna()
    main_super()
