"""FASTA input."""

import bz2
import gzip
import io
import lzma

import pytest

from restrict_cut.reader import Record, SequenceError, read_one, read_records

FASTA = ">seq1 first\nGAATTC\nAAAA\n>seq2 second\nGGATCC\n"


def write(tmp_path, text, name="input.fasta", opener=None):
    path = tmp_path / name
    if opener is None:
        path.write_text(text)
    else:
        with opener(path, "wt") as handle:
            handle.write(text)
    return path


def test_reads_several_records(tmp_path):
    records = read_records(write(tmp_path, FASTA))
    assert [record.identifier for record in records] == ["seq1", "seq2"]
    assert records[0].sequence == "GAATTCAAAA"
    assert records[0].description == "first"


def test_sequence_is_uppercased(tmp_path):
    records = read_records(write(tmp_path, ">s\ngaattc\n"))
    assert records[0].sequence == "GAATTC"


def test_whitespace_inside_a_record_is_dropped(tmp_path):
    records = read_records(write(tmp_path, ">s\nGA AT\n\nTC\n"))
    assert records[0].sequence == "GAATTC"


def test_read_one_returns_a_single_record(tmp_path):
    record = read_one(write(tmp_path, ">s\nGAATTC\n"))
    assert record.identifier == "s"


def test_read_one_refuses_several(tmp_path):
    with pytest.raises(SequenceError, match="holds several"):
        read_one(write(tmp_path, FASTA))


def test_limit_stops_early(tmp_path):
    assert len(read_records(write(tmp_path, FASTA), limit=1)) == 1


@pytest.mark.parametrize(
    "opener,name",
    [(gzip.open, "a.fasta.gz"), (bz2.open, "b.fasta.bz2"), (lzma.open, "c.fasta.xz")],
)
def test_compression_is_recognised_from_content(tmp_path, opener, name):
    records = read_records(write(tmp_path, FASTA, name, opener))
    assert records[0].sequence == "GAATTCAAAA"


def test_compression_does_not_depend_on_the_file_name(tmp_path):
    """A gzip file named .fasta must still be read."""
    records = read_records(write(tmp_path, FASTA, "plain.fasta", gzip.open))
    assert len(records) == 2


def test_byte_order_mark_is_tolerated(tmp_path):
    records = read_records(write(tmp_path, "\ufeff" + FASTA))
    assert records[0].identifier == "seq1"


def test_ambiguity_codes_are_accepted(tmp_path):
    records = read_records(write(tmp_path, ">s\nACGTRYSWKMBDHVN\n"))
    assert records[0].sequence == "ACGTRYSWKMBDHVN"
    assert records[0].ambiguous_bases == 11


def test_protein_sequence_is_refused(tmp_path):
    with pytest.raises(SequenceError, match="not nucleotide sequence"):
        read_records(write(tmp_path, ">s\nMEEPQSDPSV\n"))


def test_prose_is_refused_rather_than_measured(tmp_path):
    with pytest.raises(SequenceError, match="not nucleotide sequence"):
        read_records(write(tmp_path, ">s\nNot Found: the page you requested\n"))


def test_sequence_before_a_header_is_refused(tmp_path):
    with pytest.raises(SequenceError, match="before any"):
        read_records(write(tmp_path, "GAATTC\n>s\nGAATTC\n"))


def test_empty_input_is_refused(tmp_path):
    with pytest.raises(SequenceError, match="no FASTA records"):
        read_records(write(tmp_path, ""))


def test_header_without_an_identifier_is_refused(tmp_path):
    with pytest.raises(SequenceError, match="no identifier"):
        read_records(write(tmp_path, ">\nGAATTC\n"))


def test_missing_file_is_reported(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_records(tmp_path / "absent.fasta")


def test_directory_is_reported(tmp_path):
    with pytest.raises(IsADirectoryError):
        read_records(tmp_path)


def test_a_handle_can_be_passed_directly():
    records = read_records(io.StringIO(FASTA))
    assert len(records) == 2


def test_record_header_and_length():
    record = Record("seq1", "first", "GAATTC")
    assert record.header == "seq1 first"
    assert len(record) == 6


def test_record_header_without_description():
    assert Record("seq1", "", "A").header == "seq1"


def _compressed(path, opener, text=FASTA):
    with opener(path, "wt") as handle:
        handle.write(text)
    return path.read_bytes()


@pytest.mark.parametrize(
    "opener,name",
    [(gzip.open, "a.gz"), (bz2.open, "b.bz2"), (lzma.open, "c.xz")],
)
def test_truncated_compressed_input_is_an_input_error(tmp_path, opener, name):
    """A half-written gzip raised EOFError, which reached the user raw.

    Decompression is lazy, so the failure happens mid-read rather than at open,
    and EOFError is not an error any caller would think to catch from a FASTA
    reader.
    """
    whole = _compressed(tmp_path / name, opener)
    truncated = tmp_path / f"cut-{name}"
    truncated.write_bytes(whole[: len(whole) // 2])
    with pytest.raises(SequenceError, match="truncated"):
        read_records(truncated)


def test_corrupt_compressed_input_is_an_input_error(tmp_path):
    whole = _compressed(tmp_path / "a.gz", gzip.open, FASTA * 40)
    damaged = bytearray(whole)
    for index in range(20, len(damaged) - 10, 7):
        damaged[index] ^= 0xFF
    path = tmp_path / "damaged.gz"
    path.write_bytes(bytes(damaged))
    with pytest.raises(SequenceError, match="cannot be read"):
        read_records(path)


def test_a_file_claiming_to_be_gzip_but_is_not(tmp_path):
    path = tmp_path / "lying.gz"
    path.write_bytes(b"\x1f\x8b" + b"not actually compressed" * 4)
    with pytest.raises(SequenceError, match="cannot be read"):
        read_records(path)


def test_utf16_input_is_reported_as_not_utf8(tmp_path):
    """A sequence saved from a Windows editor is the usual way to meet this."""
    path = tmp_path / "utf16.fasta"
    path.write_bytes(FASTA.encode("utf-16"))
    with pytest.raises(SequenceError, match="not UTF-8"):
        read_records(path)


def test_binary_input_is_reported_as_not_utf8(tmp_path):
    path = tmp_path / "binary.fasta"
    path.write_bytes(bytes(range(256)) * 4)
    with pytest.raises(SequenceError, match="not UTF-8"):
        read_records(path)


def test_intact_compressed_input_is_unaffected(tmp_path):
    for opener, name in [(gzip.open, "a.gz"), (bz2.open, "b.bz2"), (lzma.open, "c.xz")]:
        _compressed(tmp_path / name, opener)
        assert read_records(tmp_path / name)[0].sequence == "GAATTCAAAA"
