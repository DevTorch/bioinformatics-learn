
from monthly_plan.week_one.fasta_functions import read_fasta_file, gc_content, get_fasta_len

base_path = '/mnt/c/Users/Torchez/PycharmProjects/bioinformatics-learn/'

file_path = base_path + 'monthly_plan/week_one/data/ls_orchid.fasta'
print(file_path)
fasta_dict = read_fasta_file(file_path)
for key, value in fasta_dict.items():
    print(key, ':',  value)

print('\nGC content:')
for key, value in fasta_dict.items():
    print(f'{key}: {gc_content(value)}%')

print('\nLength:')
print(f'Length of FASTA: {get_fasta_len(fasta_dict)} chromosomes')
