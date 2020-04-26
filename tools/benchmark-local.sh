assembly=$1
polished_assembly=$2
contig_fasta=$3
reference='/work/ryought/hic/celegans/reference/vc2010.draft-20180405.pilon.fasta'

# chop contig fasta by assembly
python /work/ryought/hi-c-assembly/hic_hiker/tools/chop_contigs.py $assembly $contig_fasta $assembly.chopped.fasta

# align to reference
bwa mem -t 32 $reference $assembly.chopped.fasta > $assembly.chopped.sam

# run benchmark for aseembly
python /work/ryought/hi-c-assembly/hic_hiker/tools/calc_accuracy_local.py $contig_fasta $assembly.chopped.sam $assembly $polished_assembly > $assembly.chopped.bench
