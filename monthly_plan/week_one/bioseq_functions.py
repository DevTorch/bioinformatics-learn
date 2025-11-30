from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def gc_content(seq: str) -> float:
    sequence = Seq(seq) if isinstance(seq, str) else seq
    return float(gc_fraction(sequence) * 100)

def gc_content_max(file_path: Path, file_format: str) -> dict[str, float]:

    """
    Return a dictionary containing the maximum GC content and the corresponding sequence ID from a given FASTA file.

    Parameters:
    file_path (str): The path to the FASTA file to parse.
    file_format (str): The format of the FASTA file.

    Returns:
    dict: A dictionary containing the maximum GC content and the corresponding sequence ID from the given FASTA file.
    """
    max_gc_content = .0
    max_seq_id = ""

    with open(file_path) as handle:
        for seq_record in SeqIO.parse(handle, file_format):
            seq = str(seq_record.seq)
            gc_content_fl = gc_fraction(seq)
            if gc_content_fl > max_gc_content:
                max_gc_content = gc_content_fl
                max_seq_id = seq_record.id

    return {max_seq_id: max_gc_content}

def kmer_count(seq: str, kmer: str) -> int:
    sequence = Seq(seq)
    return sequence.count(kmer)

def kmer_count_overlap(seq: str, kmer: str) -> int:
    sequence = Seq(seq)
    return sequence.count_overlap(kmer)

def fasta_to_pandas(file_path: Path, file_format: str) -> pd.DataFrame:
    with open(file_path) as handle:
        fasta_dict = {seq_record.id: [len(seq_record.seq), gc_fraction(seq_record.seq) * 100] for seq_record in SeqIO.parse(handle, file_format)}
    return pd.DataFrame.from_dict(fasta_dict, orient='index', columns=['Length', 'GC content'])

def plot_fasta_gc_content(file_path: Path, file_format: str) -> None:
    df = fasta_to_pandas(file_path, file_format)
    plt.figure(figsize=(18, 6))
    plt.title('GC content by sequence')
    plt.bar(df.index, df['GC content'])
    plt.ylabel('GC %')
    plt.xticks(rotation=90)
    plt.tight_layout()
    # plt.show()

def reverse_complement(seq: str) -> Seq:
    sequence = Seq(seq)
    return sequence.reverse_complement()

def transcribe(seq: str) -> Seq:
    m_rna = Seq(seq)
    return m_rna.transcribe()

def translate(seq: str) -> Seq:
    protein = Seq(seq)
    return protein.translate()