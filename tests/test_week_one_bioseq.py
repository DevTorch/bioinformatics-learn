from monthly_plan.week_one import bioseq_functions as bf


def test_bio_gc_content_basic():
    assert bf.gc_content("GCGC") == 100.0
    assert bf.gc_content("ATAT") == 0.0


def test_kmer_count():
    seq = "ATATAT"
    assert bf.kmer_count(seq, "AT") == 3  # без перекрытий
    assert bf.kmer_count_overlap(seq, "ATA") == 2  # с перекрытиями


def test_reverse_complement_bio():
    assert str(bf.reverse_complement("ATGC")) == "GCAT"


def test_transcribe():
    assert str(bf.transcribe("ATGC")) == "AUGC"


def test_translate_simple():
    # ATG -> M (Met, старт-кодон)
    protein = bf.translate("ATG")
    assert str(protein) == "M"
