import fasta_functions as ff
import bioseq_functions as bf

base_path = '/mnt/c/Users/Torchez/PycharmProjects/bioinformatics-learn/'

fasta_path = base_path + 'monthly_plan/week_one/data/ls_orchid.fasta'
ls_orchid_gbk = base_path + 'monthly_plan/week_one/data/ls_orchid.gbk'

fasta_dict = ff.read_fasta_file(fasta_path)
# for key, value in fasta_dict.items():
#     print(key, ':',  value)

fasta_get = fasta_dict.get('gi|2765603|emb|Z78478.1|PVZ78478 P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA')
print(f'\nSubsequence: {ff.get_subsequence(fasta_get, 10, 50)}')
print('\nGet info:')
print(ff.get_info(fasta_dict.get('gi|2765603|emb|Z78478.1|PVZ78478 P.victoria 5.8S rRNA gene and ITS1 and ITS2 DNA')))

print(bf.gc_content_max(fasta_path, 'fasta'))
print(bf.gc_content_max(ls_orchid_gbk, 'genbank'))

# print(ff.fasta_to_pandas(fasta_dict))

print('k-mer        : ' + bf.kmer_count(fasta_get, 'AT').__str__())
print('k-mer overlap: ' + bf.kmer_count_overlap(fasta_get, 'ATG').__str__())
print('GC content   : ' + bf.gc_content(fasta_get).__str__())

# print(bf.fasta_to_pandas(fasta_path, 'fasta'))
print(bf.reverse_complement(fasta_get))
print(bf.transcribe(fasta_get))
print(bf.translate(fasta_get))
