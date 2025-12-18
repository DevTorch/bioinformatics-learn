import os
from dotenv import load_dotenv
from Bio import Entrez
from Bio import SeqIO

load_dotenv()

email = os.getenv("EMAIL")

if not email:
    raise ValueError("EMAIL not found in .env")

Entrez.email = email
api_key = os.getenv("NCBI_API_KEY")

if api_key:
    Entrez.api_key = api_key

handle = Entrez.esearch(db="nucleotide", term="COX1[Gene]", retmax=5)
result = Entrez.read(handle)
print(result["IdList"])


## DataBase List
handle = Entrez.einfo()
rec = Entrez.read(handle)
handle.close()
print(rec)
print(rec['DbList'])

handle = Entrez.esearch(db="nucleotide", term='CRT[Gene Name] AND "Plasmodium falciparum"[Organism]', retmax=5)
result = Entrez.read(handle)
handle.close()
print(result['Count'])
print(result['IdList'])

handle = Entrez.efetch(db="nucleotide", id=result['IdList'], rettype="gb")
records = list(SeqIO.parse(handle, "gb"))
handle.close()
for record in records:
    print(record.id)
    print(record.seq)
