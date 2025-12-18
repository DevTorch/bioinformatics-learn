from collections import defaultdict
import gzip
from pathlib import Path

# import seaborn as sns
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

fasta_path = DATA_DIR / "SRR003265.filt.fastq.gz"

from Bio import SeqIO

recs = SeqIO.parse(gzip.open(fasta_path, 'rt', encoding='utf-8'),'fastq')

# for record in recs:
#     print(record.id)
#     print(record.id, record.description, record.seq)
#     print(record.letter_annotations['phred_quality'])

rec = next(recs)
print(rec.id)
print(rec.id, rec.description, rec.seq)
print(rec.letter_annotations['phred_quality'])

# counter = defaultdict(int)
# for rec in recs:
#     for letter in rec.seq:
#         counter[letter] += 1
# total = sum(counter.values())
# for letter, count in counter.items():
#     print(f'{letter}: {(count / total) * 100:.2f}%')

# n_counter = defaultdict(int)
# for rec in recs:
#     for i, letter in enumerate(rec.seq):
#         pos = i + 1
#         if letter == 'N':
#             n_counter[pos] += 1
#
# seq_len = max(n_counter.keys())
# positions = range(1, seq_len + 1)
# fig, ax = plt.subplots(figsize=(10, 6))
# ax.plot(positions, [n_counter[x] for x in positions])
# ax.set_xlim(1, seq_len)
# plt.show()

# Quality Reads: N (странный фактор-примесь) не встречается до 25-й позиции (предыдущий график)
# потому отбрасываем позиции < 25
quality_counter = defaultdict(int)
for rec in recs:
    for i, quality in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25:
            continue
        quality_counter[quality] += 1
total = sum(quality_counter.values())
for quality, count in quality_counter.items():
    print('%d: %.2f %d' % (quality, (count / total) * 100, count))

