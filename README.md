# hic_hiker v1.0.0

3D-DNA scaffolds refinement


## Installation
```
$ pip install -e .
```

This will install the python package (and its cli command) `hic_hiker`.

(In the future, `$ pip install hic_hiker`)

### Requirements

To run HiC-Hiker, you will need to install the following:

- `3D-DNA, juicer` and their requirements
- `python3, >=3.5` with
    - `numpy, matplotlib, sklearn, scipy` (basic matrix operation and visualization)
    - `pandas, feather` (handling of large datasets)
    - `BioPython` (fasta parsing)
    - `pysam` (sam parsing)
    - `tqdm` (progress bar)
    - `matplotlib-scalebar` (to show scalebar on the benchmark matrix plot)
    (These will automatically satisfied while installation using `pip`)

### System Requirements

- RAM >100GB for human dataset

## Usage

### Required Files

You have to prepare

- Input of 3D-DNA:
    - contigs: `.fasta`
    - Hi-C contacts: `.mnd.txt` (or `.R1.sam, .R2.sam`)
- Output of 3D-DNA:
    - scaffold layout: `.final.assembly`
        You are recommended to use not `.FINAL.assembly` but `.final.assembly`
    - (chopped contigs: `.final.fasta`)
- A directory for workspace (to store intermediate files or results)

Additionally, to run benchmarks with the reference sequence (and to plot the figures on our paper), you will need
- alignment of contigs to the reference: `.sam`

by running
```
$ bwa mem -t 32 ../hg38/hg38.fa GSE95797_Hs1.final.fasta > GSE95797_Hs1.final.fasta.sam
```
where `GSE95797_Hs1.final.fasta` is one of the outputs of 3D-DNA. (chopped contigs)

### Pipeline
```
$ hic_hiker <contigs.fasta> <scaffold_layout.assembly> <contacts.mnd.txt> <workspace directory> <K>
```

After the process finishes, you will see in workspace directory:

- `polished.assembly` assembly file with refinement of orientations
- `polished.fasta` polished chromosome-length scaffold sequences (with no gap added)
- figures
    - `fig_distribution.png`
    - `fig_errorchart.png`
    - `fig_matrix.png`

## Uninstallation
```
$ pip uninstall hic_hiker
```

## Citation

HiC-Hiker: A probabilistic model to determine contig orientation in chromosome-length scaffolds by Hi-C (not published)
