#!/usr/bin/env python
# -*- coding: utf-8 -*-
from joblib import Parallel, delayed
import itertools

def hoge(array):
    return [x*5 for x in array]

def process(i, j, array):
    return (i, sum(hoge(array)))

def main():
    print('start')
    array = [[1,5,10,215],[15,412,424],[524,553,4]]
    jobs = [delayed(process)(i, j, array[j]) \
        for i,j in itertools.combinations([0,1,2], 2)]
    r = Parallel(n_jobs=-1, verbose=3)(jobs)
    print(r)
    print('finished')

if __name__ == '__main__':
    main()
