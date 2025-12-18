from collections import defaultdict
import gzip
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

fastq_path = DATA_DIR / "SRR003265.filt.fastq.gz"

recs = SeqIO.parse(gzip.open(fastq_path, 'rt', encoding='utf-8'), 'fastq')

quality_pos = defaultdict(list)

for rec in recs:
    for i, quality in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25 or quality == 40:
            continue
        pos = i + 1
        quality_pos[pos].append(quality)

poses = sorted(quality_pos.keys())
vps = [quality_pos[pos] for pos in poses]

fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=vps, ax=ax)

ax.set_xticks(range(len(poses)))
ax.set_xticklabels([str(pos) for pos in poses])
plt.show()
