import os
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import molecular_weight
from dotenv import load_dotenv

load_dotenv()

email = os.getenv("EMAIL")

if not email:
    raise ValueError("EMAIL not found in .env")

Entrez.email = email
api_key = os.getenv("NCBI_API_KEY")

if api_key:
    Entrez.api_key = api_key

handle = Entrez.efetch(db="nucleotide", id="MN908947", rettype="gb", retmode="text")
records = list(SeqIO.parse(handle, "gb"))
handle.close()
# print(records)

covid_dna = records[0].seq
print(f'The genome of COVID-19 consists of {len(covid_dna)} nucleotides')
print(f'The molecular weight of COVID-19 is {molecular_weight(covid_dna, "DNA")}')
print(f'The GC content of COVID-19 is {gc_fraction(covid_dna) * 100}')  # Чем выше GC-content, тем выше стабильность DNA

nucleotides_counter = {
    'A': covid_dna.count('A'),
    'T': covid_dna.count('T'),
    'C': covid_dna.count('C'),
    'G': covid_dna.count('G')
}
print(f'The number of nucleotides in COVID-19 is {nucleotides_counter}')

width = 0.8
plt.bar(nucleotides_counter.keys(), nucleotides_counter.values(), width, color=['blue', 'green', 'red', 'cyan'])
plt.xlabel('Nucleotides')
plt.ylabel('Frequency')
plt.title('Nucleotides frequency in COVID-19')
# plt.show()

# Транскрипция
covid_dna_full_length = len(covid_dna) - (len(covid_dna) % 3)
trimmed_covid_dna = covid_dna[:covid_dna_full_length]
m_rna = trimmed_covid_dna.transcribe()
print(m_rna)

# Трансляция
protein = m_rna.translate()
print(protein)

# Распределение аминокислот
amino_acids_counter = Counter(protein)
print(amino_acids_counter.most_common())

del amino_acids_counter['*']
width = 0.5
plt.bar(amino_acids_counter.keys(), amino_acids_counter.values(), width,
        color=['blue', 'green', 'red', 'cyan', 'yellow', 'black'])
plt.xlabel('Amino acid')
plt.ylabel('Frequency')
plt.title('Proteins frequency in COVID-19')
# plt.show()

proteins = protein.split('*')
print(proteins)
print(f'Number of proteins: {len(proteins)}')

clean_protein = []
for protein in proteins:
    if len(protein) > 20:  # Только последовательности, содержащие более 20 аминокислот, являются белками
        clean_protein.append(protein)

print(clean_protein)
print(f'Number of proteins now: {len(clean_protein)}')
print(sorted(clean_protein, key=len, reverse=True))

sorted_proteins = sorted(clean_protein, key=len)
print(len(sorted_proteins[-1]))

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
fasta_file = DATA_DIR / 'longest_orf.fasta'

with open(fasta_file, 'w') as f:
    f.write(f'>COVID PROTEIN\n{sorted_proteins[-1]}')
