"""Ambiguity codes, complementation and the matching rule."""

import pytest

from restrict_cut import iupac


def test_every_code_has_bases():
    for code, bases in iupac.CODES.items():
        assert bases, code
        assert bases <= frozenset("ACGT")


def test_reverse_complement_round_trips():
    sequence = "ACGTRYSWKMBDHVN"
    assert iupac.reverse_complement(iupac.reverse_complement(sequence)) == sequence


def test_reverse_complement_known_sites():
    assert iupac.reverse_complement("GAATTC") == "GAATTC"
    assert iupac.reverse_complement("GGTCTC") == "GAGACC"
    assert iupac.reverse_complement("RGATCY") == "RGATCY"


def test_reverse_complement_preserves_case():
    assert iupac.reverse_complement("gaattc") == "gaattc"
    assert iupac.reverse_complement("AcGt") == "aCgT"


def test_ambiguity_codes_complement_consistently():
    """R (A or G) must complement to Y (C or T), base set for base set."""
    for code, bases in iupac.CODES.items():
        if code == "U":
            continue
        complemented = iupac.complement(code)
        expected = {iupac.complement(base) for base in bases}
        assert iupac.CODES[complemented] == frozenset(expected), code


@pytest.mark.parametrize(
    "pattern,expected",
    [
        ("GAATTC", True),
        ("RGATCY", True),
        ("GTYRAC", True),
        ("CASTG", True),
        ("GGTCTC", False),
        ("CCGC", False),
        ("GCANNNNNNTGC", True),
        ("CGANNNNNNTGC", False),
    ],
)
def test_palindrome_detection(pattern, expected):
    assert iupac.is_palindrome(pattern) is expected


def test_matches_is_containment_not_overlap():
    assert iupac.matches("N", "A")
    assert not iupac.matches("A", "N")
    assert iupac.matches("R", "G")
    assert not iupac.matches("R", "C")
    assert not iupac.matches("A", "R")


def test_n_in_sequence_satisfies_n_in_site():
    """The bug this guards: an N in the sequence must match an N in the site."""
    assert iupac.matches("N", "N")
    assert iupac.matches("N", "R")
    assert iupac.matches("R", "R")


def test_could_match_is_overlap():
    assert iupac.could_match("A", "N")
    assert iupac.could_match("R", "N")
    assert not iupac.could_match("A", "C")
    assert not iupac.could_match("R", "Y")


def test_matches_implies_could_match():
    for site in iupac.CODES:
        for sequence in iupac.CODES:
            if iupac.matches(site, sequence):
                assert iupac.could_match(site, sequence), (site, sequence)


def test_expand_rejects_non_codes():
    with pytest.raises(KeyError):
        iupac.expand("Z")


def test_expand_is_case_insensitive():
    assert iupac.expand("n") == iupac.expand("N")


def test_uracil_reads_as_thymine():
    assert iupac.expand("U") == iupac.expand("T")


def test_has_ambiguity_agrees_with_the_code_table():
    """The fast path must agree with expand() on every single character.

    has_ambiguity decides whether uncertain positions are looked for at all, so
    a False where it should be True does not slow the tool down, it silently
    stops it reporting uncertainty.
    """
    for code in iupac.CODES:
        expected = len(iupac.expand(code)) > 1
        assert iupac.has_ambiguity(code) is expected, code
        assert iupac.has_ambiguity(code.lower()) is expected, code


def test_has_ambiguity_on_sequences():
    assert not iupac.has_ambiguity("ACGTACGT")
    assert not iupac.has_ambiguity("acgtacgt")
    assert not iupac.has_ambiguity("")
    assert iupac.has_ambiguity("ACGTN")
    assert iupac.has_ambiguity("ACGTr")


def test_rna_is_not_treated_as_ambiguous():
    """U names one base, so an RNA sequence carries no uncertainty."""
    assert not iupac.has_ambiguity("AUGCAUGC")


def test_has_ambiguity_matches_a_per_character_scan():
    sequence = "ACGTURYSWKMBDHVNacgtn"
    expected = any(len(iupac.expand(base)) > 1 for base in sequence)
    assert iupac.has_ambiguity(sequence) is expected
