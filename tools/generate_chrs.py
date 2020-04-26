import csv
import sys
def parse_mums(filename):
    ret = []
    query_store = None
    single_ref = False
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for l in reader:
            if l[0] == '>':  # start flag
                length = int(' '.join(l[1:]).split('\t')[-1].split(' ')[-1])
                name_with_flag = ''.join(' '.join(l[1:]).split('\t')[:-1])
                L = name_with_flag.split(' ')
                if L[-1] != 'Reverse':  # start line of match
                    if query_store:
                        ret.append(query_store)
                    query_store = {
                        'query_name': ' '.join(L),
                        'length': length,
                        '+': [],
                        '-': []
                    }
                    direction = '+'
                else:  # reverse match is reported after forward match
                    direction = '-'
            else:
                if len(l) > 3:
                    # (ref_id, refでの開始位置, queryでの開始位置, match length)
                    query_store[direction].append((l[0], int(l[1]), int(l[2]), int(l[3])))
                    single_ref = False
                else:
                    # (ref_id = 'query', refでの開始位置, queryでの開始位置, match length)
                    query_store[direction].append(('reference', int(l[0]), int(l[1]), int(l[2])))
                    single_ref = True
        ret.append(query_store)
    return ret, single_ref

def main():
    if len(sys.argv) != 2:
        print('usage: python3 generate_chrs.py hoge.mum > hoge.chrs')
        sys.exit()
    mum_filename = sys.argv[1]
    print('loading', mum_filename, file=sys.stderr)
    mum, _ = parse_mums(mum_filename)

    # minimum length
    # longer than 1M
    L = 1000000
    mum_long = [m for m in mum if m['length'] >= L]
    print(len(mum_long), file=sys.stderr)
    from collections import defaultdict
    for long_scaf in mum_long:
        # print(long_scaf['query_name'])
        d = defaultdict(lambda: 0)
        for segment in long_scaf['+']:
            chr_name, _, _, match_length = segment
            d[(chr_name, '+')] += match_length
        for segment in long_scaf['-']:
            chr_name, _, _, match_length = segment
            d[(chr_name, '-')] += match_length
        match_lengths = sorted(d.items(), key=lambda x: x[1])[::-1]
        if len(match_lengths) == 0:
            print('no alignment! {}'.format(long_scaf['query_name'], file=sys.stderr))
            print('{},{},{}'.format(long_scaf['query_name'], '?', '?'))
        elif len(match_lengths) > 0:
            print(long_scaf['query_name'], match_lengths[0], file=sys.stderr)
            (chr_name, ori), match_length = match_lengths[0]
            print('{},{},{}'.format(long_scaf['query_name'], chr_name, False if ori == '+' else True))
            if len(match_lengths) > 1:
                print(long_scaf['query_name'], match_lengths[0], match_lengths[0][1] / match_lengths[1][1], file=sys.stderr)
                if match_lengths[0][1] / match_lengths[1][1] < 3:
                    print('error rate seems to be high', file=sys.stderr)

        # from collections import Counter
        # cp = Counter([segment[0] for segment in long_scaf['+']])
        # cn = Counter([segment[0] for segment in long_scaf['-']])
        # # cp_top = cp.most_common()[0]
        # cn_top = cn.most_common()[0]
        # print('+', cp.items())
        # print('-', cn.items())

main()
