#!/bin/bash
set -Ceuo pipefail
cd $(dirname $0)
K=100
N=10000
python random_sample_subsequences.py hoge.fasta fuga.fasta $K $N

bwa mem reference.fasta fuga.fasta > reference.sam
bwa mem assembly.fasta fuga.fasta > assembly.sam


# get uniquely mapped fragments


# compare their orientations
