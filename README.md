# hic-hiker
orientation

## modules
- main.py
    CUI app module
- load.py
    from sam file to contacts.pkl (parser)
- prob.py
    from contacts.pkl to prob.npy (kde and polyfit)
- hmm.py
    from prob.npy to most probable path (hmm)

## descritption of intermidiate files
- contacts.pkl
    `contacts[i][j], (i<j)` is a list of contacts between contig i and j. `contacts[i][j][0]` is ndarray of position of each contact in contig i.
- prob.npy
    `prob[2*i+(0or1), 2*j+(0or1)]` represents probability of layout of contig i and j. (0or1) is orientation.

## how to run
### whole pipeline
1. generate contigs with an assembler you prefer
    We assume you get `contigs.fasta`
2. mapping
```
bwa mem -t 16 contigs.fasta R1.fastq > contigs.R1.sam
bwa mem -t 16 contigs.fasta R2.fastq > contigs.R2.sam
```
3. run
`python main.py contigs.R1.sam contigs.R2.sam`

### each part
`python [script] [input] [output]`


## parameter
- debug
