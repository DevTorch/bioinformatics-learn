import os
from typing import Any


def read_fasta_file(file_path: str) -> dict:
    """
    Read a FASTA file and return a dictionary containing the sequences.

    Parameters:
    file_path (str): The path to the FASTA file to read.

    Returns:
    dict: A dictionary containing the sequences from the FASTA file.

    Raises:
    FileNotFoundError: If the file at the given path does not exist.
    """
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
    """
    Return the GC content of a given sequence.

    The GC content is the percentage of G and C nucleotides in a given sequence.

    Parameters:
    seq (str): The sequence to calculate the GC content for.

    Returns:
    str: The GC content of the given sequence as a string with two decimal places.
    """
    return format((seq.count('G') + seq.count('C')) / len(seq) * 100, '.2f')

def gc_content_dict(fasta_dict: dict) -> dict:
    """
    Return a dictionary with GC content for each sequence in a given FASTA file.
    The returned dictionary is sorted by GC content in descending order.

    Parameters:
    fasta_dict (dict): A dictionary containing sequences from a FASTA file.

    Returns:
    dict: A dictionary with GC content for each sequence in the given FASTA file.

    """
    gc_content_dictionary = {label: gc_content(seq) for label, seq in fasta_dict.items()}
    return dict(sorted(gc_content_dictionary.items(), key=lambda x: x[1], reverse=True))

def get_fasta_len(fasta_dict: dict):
    return len(fasta_dict)

def genomes_summary_len(fasta_dict: dict):
    return sum(len(seq) for seq in fasta_dict.values())

def genomes_len(fasta_dict: dict):
    return {label: len(seq) for label, seq in fasta_dict.items()}

def genome_len(seq: str):
    return len(seq)

def reverse_complement(seq: str):
    # return seq[::-1].replace('A', 'T').replace('T', 'A').replace('C', 'G').replace('G', 'C')
    mapping = str.maketrans("ATCG", "TAGC")
    return seq.translate(mapping)[::-1]

def get_subsequence(seq: str, start: int = 0, end: int = -1) -> str:
    """
    Return a subsequence of a given sequence.

    Parameters:
    seq (str): The sequence to get the subsequence from.
    start (int): The start index of the subsequence. Defaults to 0.
    end (int): The end index of the subsequence. Defaults to -1, which means the end of the sequence.

    Returns:
    str: The subsequence of the given sequence.
    """
    if end == -1:
        end = len(seq)
    return seq[start:end] if start < end < len(seq) else "Invalid indices"

def get_info(seq: str):
    return (
        f'Length            : {len(seq)}\n'
        f'GC content        : {gc_content(seq)}%\n'
        f'Sequence          : {seq}\n'
        f'Reverse complement: {reverse_complement(seq)}\n'
    )