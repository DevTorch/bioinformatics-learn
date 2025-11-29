import fasta_functions as ff

base_path = '/mnt/c/Users/Torchez/PycharmProjects/bioinformatics-learn/'

file_path = base_path + 'monthly_plan/week_one/data/ls_orchid.fasta'

fasta_dict = ff.read_fasta_file(file_path)
# for key, value in fasta_dict.items():
#     print(key, ':',  value)

fasta_get = fasta_dict.get('gi|2765603|emb|Z78478.1|PVZ78478 P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA')
print(f'\nSubsequence: {ff.get_subsequence(fasta_get, 10, 50)}')
print('\nGet info:')
print(ff.get_info(fasta_dict.get('gi|2765603|emb|Z78478.1|PVZ78478 P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA')))