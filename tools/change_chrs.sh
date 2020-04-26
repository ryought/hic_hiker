set -Ceuo pipefail
[ $# -ne 1 ] && echo 'not enough arguments' && exit 1
chrs=$1

cat $chrs \
  | sed -e "s/HiC_scaffold_1/hicscaffold_000000/g" \
  | sed -e "s/HiC_scaffold_2/hicscaffold_000001/g" \
  | sed -e "s/HiC_scaffold_3/hicscaffold_000002/g" \
  | sed -e "s/HiC_scaffold_4/hicscaffold_000003/g" \
  | sed -e "s/HiC_scaffold_5/hicscaffold_000004/g" \
  | sed -e "s/HiC_scaffold_6/hicscaffold_000005/g"
