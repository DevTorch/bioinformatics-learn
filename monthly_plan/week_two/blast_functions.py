from Bio.Blast import NCBIWWW, NCBIXML

# help(NCBIWWW)

result_handle = NCBIWWW.qblast("blastn", "nt", """ggtaagtcctctagtacaaacacccccaatattgtgatataattaaa
attatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc""")
blast_records = NCBIXML.parse(result_handle)

E_VALUE_THRESH = 0.00000000001
count = 0
for record in blast_records:
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                count += 1
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
                print()

print(f"There are {count} similar sequences in Blast output")
