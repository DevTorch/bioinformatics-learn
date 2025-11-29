import os

def read_fasta_file(file_path: str) -> dict:
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found")
    with open(file_path, "r") as f:
        fasta_file = [line.strip() for line in f.readlines()]

    fasta_dict = {}
    fasta_label = ""

    for line in fasta_file:
        if line.startswith(">"):
            fasta_label = line[1:]
            fasta_dict[fasta_label] = ""
        else:
            fasta_dict[fasta_label] += line
    return fasta_dict

def gc_content(seq):
    """Содержание GC в последовательности ДНК/РНК"""
    return round((seq.count('G') + seq.count('C')) / len(seq) * 100)

def get_fasta_len(fasta_dict: dict):
    return len(fasta_dict)
