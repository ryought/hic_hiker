
# cat - << EOS > mock.kv
# hoge1	chr1
# hoge2	chr4
# EOS

echo -e "hoge\tchr1\nhoge2\tchr4" > mock.kv

samtools faidx mock.fasta
samtools faidx mock.fasta hoge hoge2 > mock.picked.fasta
seqkit replace -p '(.+)$' -r '{kv}' -k mock.kv  mock.picked.fasta
