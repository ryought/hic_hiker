import hic_hiker
import sys

def main():
    if len(sys.argv) != 5:
        print('usage: python3 calc_accuracy_local.py contigs.fasta chopped.sam original.assembly polished.assembly')
        sys.exit(1)
    contigs_fasta_filename = sys.argv[1]
    chopped_sam_filename = sys.argv[2]
    assembly_filename = sys.argv[3]
    polished_assembly_filename = sys.argv[4]

    contigs = hic_hiker.contigs.Contigs(fasta_filename=contigs_fasta_filename)
    asm = hic_hiker.load_3ddna.Assembly(
            asm_filename=assembly_filename,
            contigs=contigs,
            )
    layout = asm.get_layout()

    chopped_contigs = asm.get_new_contigs(contigs)
    ref_layout = hic_hiker.layout.get_reference_layout_from_sam(
            chopped_sam_filename,
            chopped_contigs)
    result = hic_hiker.benchmark.determine_correct_orientation_or_not(chopped_contigs, layout, ref_layout)
    F, T = hic_hiker.benchmark.parse_error_rate_from_result(result)
    print(F, T, F / (T + F), T / (T + F))

    asm_polished = hic_hiker.load_3ddna.Assembly(
            asm_filename=polished_assembly_filename,
            contigs=contigs,
            )
    layout_polished = asm_polished.get_layout()
    result_polished = hic_hiker.benchmark.determine_correct_orientation_or_not(chopped_contigs, layout_polished, ref_layout)
    F, T = hic_hiker.benchmark.parse_error_rate_from_result(result_polished)
    print(F, T, F / (T + F), T / (T + F))

    TP, TN, FP, FN = 0, 0, 0, 0
    for i in range(len(result[0])):
        for x in zip(result[0][i], result_polished[0][i]):
            if x[0] == hic_hiker.benchmark.ORI_ERR and x[1] == hic_hiker.benchmark.OK:
                TP += 1
            elif x[0] == hic_hiker.benchmark.ORI_ERR and x[1] == hic_hiker.benchmark.ORI_ERR:
                FN += 1
            elif x[0] == hic_hiker.benchmark.OK and x[1] == hic_hiker.benchmark.OK:
                TN += 1
            elif x[0] == hic_hiker.benchmark.OK and x[1] == hic_hiker.benchmark.ORI_ERR:
                FP += 1
    print('TP, TN, FP, FN')
    print(TP, TN, FP, FN)

main()
