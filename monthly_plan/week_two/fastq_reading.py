from Bio import SeqIO
from pathlib import Path
import time

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
fastq_file = DATA_DIR / 'SRR494102.fastq'
orchid_fasta_file = DATA_DIR / 'ls_orchid.fasta'

# Индексное чтение
start = time.perf_counter()
fastq_dict = SeqIO.index(fastq_file, "fastq")
end = time.perf_counter()

elapsed_time = end - start
num_records = len(fastq_dict)

print(f"Проиндексировано записей: {num_records:}")
print(f"Время индексации: {elapsed_time:.4f} секунд")
print(f"Скорость: {num_records / elapsed_time:.1f} записей/сек")

# Последовательность ключей
keys_ = list(fastq_dict.keys())[:5]
print(f'Ключи (5): {keys_}')

fastq_dict_srr____seq = fastq_dict['SRR494102.20000']
print(f'Доступ по индексу: {fastq_dict_srr____seq}')

