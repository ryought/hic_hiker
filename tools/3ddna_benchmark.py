# filename='../../analysis/3D_DNA/answer70gap/vc2010_25000.final.assembly'
# filename='../../analysis/3D_DNA/mock/mock.final.assembly'
filename='../../analysis/3D_DNA/mock_chrI/mock_chrI.final.assembly'
print(filename)
C = 0
T = 0
with open(filename) as f:
    for l in f:
        if l[0] != '>':
            ls = l.split(' ')
            A, B = 0, 0
            if len(ls) > 1:
                for contig in ls:
                    if contig[0] != '-':
                        A += 1
                    else:
                        B += 1
                print(A, B)
                T += min(A, B)
            else:
                C += 1
print(C, T)
