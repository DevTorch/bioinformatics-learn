from Bio import Align
from Bio.Align import substitution_matrices

aligner = Align.PairwiseAligner()
aligner.mode = 'global'  # или 'local' для локального выравнивания
aligner.match_score = 1.0
aligner.mismatch_score = 0.0
aligner.open_gap_score = 0.0
aligner.extend_gap_score = 0.0

alignments = aligner.align("ACCGGT", "ACGT")
for alignment in alignments:
    print(alignment)
    print(f"Score: {alignment.score}")
    print()

matrix = substitution_matrices.load("BLOSUM62")
aligner.substitution_matrix = matrix
alignments_new = aligner.align("KEVLA", "EVL")
for alignment in alignments_new:
    print(alignment)
    print(f"Score: {alignment.score}")
    print()


def print_alignment(aligner, seq1, seq2, title=""):
    """Утилита для печати выравниваний"""
    print(f"\n{'=' * 50}")
    print(f"{title}")
    print(f"{'=' * 50}")

    alignments = aligner.align(seq1, seq2)
    for i, alignment in enumerate(alignments[:3]):  # первые 3 выравнивания
        print(f"Alignment #{i + 1}:")
        print(alignment)
        print(f"Score: {alignment.score}")
        print("-" * 30)


# 1. Global (похоже на globalxx)
aligner_global = Align.PairwiseAligner()
aligner_global.mode = 'global'
aligner_global.match_score = 1.0
aligner_global.mismatch_score = 0.0
aligner_global.open_gap_score = 0.0
aligner_global.extend_gap_score = 0.0

# print_alignment(aligner_global, "ACCGGT", "ACGT", "Global alignment (like globalxx)")

# 2. Local (похоже на localxx)
aligner_local = Align.PairwiseAligner()
aligner_local.mode = 'local'
aligner_local.match_score = 1.0
aligner_local.mismatch_score = 0.0
aligner_local.open_gap_score = 0.0
aligner_local.extend_gap_score = 0.0

# print_alignment(aligner_local, "ACCGGT", "ACGT", "Local alignment (like localxx)")

# 3. С штрафом за гэпы (похоже на globalms)
aligner_gaps = Align.PairwiseAligner()
aligner_gaps.mode = 'global'
aligner_gaps.match_score = 2.0
aligner_gaps.mismatch_score = -1.0
aligner_gaps.open_gap_score = -0.5
aligner_gaps.extend_gap_score = -0.1

# print_alignment(aligner_gaps, "ACCGGT", "ACGT", "Global with gap penalties")