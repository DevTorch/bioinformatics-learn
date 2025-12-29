from pathlib import Path
from Bio import motifs, SeqIO
from Bio.Seq import Seq

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
fasta_file = DATA_DIR / 'ls_orchid.fasta'

instances = []
with open(fasta_file) as handle:
    records = list(SeqIO.parse(handle, "fasta"))
    for record in records:
        instances.append(record.seq)
        # print(record.id)
        # print(record.seq)
# Это работает только на последовательностях одинаковой длины
# motif = motifs.create(instances)
# print(motif)

hardcoded_instances = [
    Seq("ATCGTACGTA"),
    Seq("ATCGTACGTA"),
    Seq("ATCGTACGTA"),
    Seq("ATCGTACGTA"),
    Seq("ATCGTACGTA")
]
motif = motifs.create(hardcoded_instances)
print(motif)
print(motif.counts)
print(motif.consensus)
print(motif.anticonsensus)
print(motif.degenerate_consensus)
print(motif.pssm)

