import textwrap
from pathlib import Path

import pandas as pd
import pytest

from monthly_plan.week_one import fasta_functions as ff


def test_gc_content_basic():
    assert ff.gc_content("GCGC") == pytest.approx(100.0)
    assert ff.gc_content("ATAT") == pytest.approx(0.0)
    assert ff.gc_content("ATGC") == pytest.approx(50.0)


def test_gc_content_empty_raises():
    with pytest.raises(ValueError):
        ff.gc_content("")


def test_reverse_complement_basic():
    assert ff.reverse_complement("ATGC") == "GCAT"
    assert ff.reverse_complement("A") == "T"
    assert ff.reverse_complement("GGCC") == "GGCC"


def test_get_subsequence_default_end():
    seq = "ATGCGT"
    subseq = ff.get_subsequence(seq, 1, -1)
    assert subseq == "TGCGT"


def test_get_subsequence_invalid_indices():
    seq = "ATGC"
    with pytest.raises(ValueError):
        ff.get_subsequence(seq, 3, 1)  # start > end


def test_read_fasta_file_creates_dict(tmp_path: Path):
    # создаём мини-FASTA прямо в тесте
    fasta_content = textwrap.dedent(
        """
        >seq1
        ATGC
        >seq2
        GGGG
        """
    ).strip()

    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    result = ff.read_fasta_file(fasta_file)

    assert isinstance(result, dict)
    assert set(result.keys()) == {"seq1", "seq2"}
    assert result["seq1"] == "ATGC"
    assert result["seq2"] == "GGGG"


def test_fasta_to_pandas(tmp_path: Path):
    fasta_content = textwrap.dedent(
        """
        >seq1
        ATGC
        >seq2
        GGGG
        """
    ).strip()

    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)

    fasta_dict = ff.read_fasta_file(fasta_file)
    df = ff.fasta_to_pandas(fasta_dict)

    assert isinstance(df, pd.DataFrame)
    assert set(df.index) == {"seq1", "seq2"}
    assert list(df.columns) == ["Length", "GC content"]
    assert df.loc["seq1", "Length"] == 4
