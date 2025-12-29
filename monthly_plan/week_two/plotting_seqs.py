from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
fasta_file = DATA_DIR / 'ls_orchid.fasta'

sequences = []
sizes = []
for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(record.seq)
    sizes.append(len(record.seq))

# Histogram of sequence lengths
plt.hist(sizes, bins=20)
plt.title(f'{len(sizes)} orchid sequences\n {min(sizes)} to {max(sizes)}')
plt.xlabel('Sequence length (bp)')
plt.ylabel('Count')
plt.show()

# GC content
gc_content = sorted((gc_fraction(seq))  for seq in sequences)
plt.plot(gc_content)
plt.title(f'{len(gc_content)} orchid sequences\nGC Content(%) {min(gc_content)} to {max(gc_content)}')
plt.xlabel('Genes')
plt.ylabel('GC content')
plt.show()