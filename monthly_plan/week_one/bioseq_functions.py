from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt

def gc_content(seq: str):
    sequence = Seq(seq)
    return gc_fraction(sequence) * 100

def gc_content_max(file_path: str, file_format: str):

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

    with SeqIO.parse(file_path, file_format) as seq_records:
        for seq_record in seq_records:
            seq = str(seq_record.seq)
            gc_content = (seq.count('G') + seq.count('C')) / len(seq) * 100
            if gc_content > max_gc_content:
                max_gc_content = gc_content
                max_seq_id = seq_record.id

    return {max_seq_id: max_gc_content}

def kmer_count(seq: str, kmer: str):
    sequence = Seq(seq)
    return sequence.count(kmer)

def kmer_count_overlap(seq: str, kmer: str):
    sequence = Seq(seq)
    return sequence.count_overlap(kmer)

def fasta_to_pandas(file_path: str, file_format: str):
    with SeqIO.parse(file_path, file_format) as seq_records:
        fasta_dict = {seq_record.id: [len(seq_record.seq), gc_fraction(seq_record.seq) * 100] for seq_record in seq_records}
    return pd.DataFrame.from_dict(fasta_dict, orient='index', columns=['Length', 'GC content'])

def plot_fasta_gc_content(file_path: str, file_format: str):
    df = fasta_to_pandas(file_path, file_format)
    plt.figure(figsize=(18, 6))
    plt.title('GC content')
    plt.bar(df.index, df['GC content'])
    plt.show()

def reverse_complement(seq: str):
    sequence = Seq(seq)
    return sequence.reverse_complement()

def transcribe(seq: str):
    m_rna = Seq(seq)
    return m_rna.transcribe()

def translate(seq: str):
    protein = Seq(seq)
    return protein.translate()