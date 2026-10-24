"""Figures computed from the bundled reference sequences.

These pin the tool against real, publicly deposited DNA rather than against
contrived input. Each number here was checked against the published record or
derived independently of the library before being written down; recompute
before quoting any of them elsewhere.
"""

import re

import pytest

from restrict_cut import enzymes
from restrict_cut.digest import digest
from restrict_cut.reader import read_one

#: pUC19 polylinker single cutters, as the top-strand nick each leaves. The set
#: and its order along the sequence are the defining property of the pUC19
#: polylinker, which is what makes them worth pinning.
PUC19_POLYLINKER = {
    "EcoRI": 396,
    "SacI": 406,
    "KpnI": 412,
    "XmaI": 412,
    "SmaI": 414,
    "BamHI": 417,
    "XbaI": 423,
    "SalI": 429,
    "PstI": 439,
    "SphI": 445,
    "HindIII": 447,
}

#: HindIII fragments of the deposited lambda sequence. Note that this is seven
#: fragments, not the eight of the commercial lambda size marker: see
#: test_lambda_lacks_the_marker_site_that_would_split_6682.
LAMBDA_HINDIII = [23130, 9416, 6682, 4361, 2322, 2027, 564]


@pytest.fixture(scope="module")
def puc19(examples):
    return read_one(examples / "puc19.fasta")


@pytest.fixture(scope="module")
def lam(examples):
    return read_one(examples / "lambda.fasta")


def test_example_lengths_match_the_records(examples):
    for filename, length in [
        ("puc19.fasta", 2686),
        ("pbr322.fasta", 4361),
        ("lambda.fasta", 48502),
        ("aav2.fasta", 4679),
    ]:
        assert len(read_one(examples / filename)) == length, filename


def test_examples_are_unmodified_records(examples):
    """Each example keeps the accession it was downloaded under."""
    for filename, accession in [
        ("puc19.fasta", "L09137.2"),
        ("pbr322.fasta", "J01749.1"),
        ("lambda.fasta", "NC_001416.1"),
        ("aav2.fasta", "NC_001401.2"),
    ]:
        assert read_one(examples / filename).identifier == accession


@pytest.mark.parametrize("name,position", sorted(PUC19_POLYLINKER.items()))
def test_puc19_polylinker_enzymes_cut_once(puc19, name, position):
    result = digest(puc19.sequence, [enzymes.get(name)], circular=True)
    assert [cut.top for cut in result.cuts] == [position]


def test_puc19_polylinker_order_is_preserved(puc19):
    """The polylinker runs EcoRI through HindIII, in that order."""
    order = [name for name, _ in sorted(PUC19_POLYLINKER.items(), key=lambda kv: kv[1])]
    assert order[0] == "EcoRI"
    assert order[-1] == "HindIII"


def test_puc19_pvuii_gives_the_published_pair(puc19):
    """PvuII cuts pUC19 twice, into 322 bp and 2364 bp."""
    result = digest(puc19.sequence, [enzymes.get("PvuII")], circular=True)
    assert sorted((f.length for f in result.fragments), reverse=True) == [2364, 322]
    assert result.total_length == 2686


def test_puc19_single_cutters_really_cut_once(puc19):
    """Every enzyme reported as a single cutter must cut exactly once."""
    panel = enzymes.select(available_only=True, min_length=6)
    from restrict_cut.selection import survey

    found = survey(puc19.sequence, panel, circular=True)
    for name in found.cutting(1):
        assert len(digest(puc19.sequence, [enzymes.get(name)], circular=True).cuts) == 1


def test_lambda_hindiii_fragments(lam):
    result = digest(lam.sequence, [enzymes.get("HindIII")], circular=False)
    assert sorted((f.length for f in result.fragments), reverse=True) == LAMBDA_HINDIII
    assert result.total_length == 48502


def test_lambda_hindiii_site_count_agrees_with_a_plain_text_search(lam):
    """Cross-check against a search that uses none of this library's code."""
    plain = [match.start() + 1 for match in re.finditer("(?=AAGCTT)", lam.sequence)]
    result = digest(lam.sequence, [enzymes.get("HindIII")])
    assert len(plain) == 6
    assert [site.start for site in result.sites] == plain


def test_lambda_lacks_the_marker_site_that_would_split_6682(lam):
    """The deposited sequence is not the commercial marker's lambda.

    The lambda-HindIII size standard is quoted as eight fragments including a
    6557 and a 125, which sum to 6682. Reaching those needs a seventh HindIII
    site near position 37584, and the deposited sequence reads AAGCTC there, so
    six sites and seven fragments is the correct answer for this file. The test
    records the difference rather than papering over it.
    """
    assert lam.sequence[37583: 37589] == "AAGCTC"
    assert 6557 + 125 == 6682
    assert 6682 in LAMBDA_HINDIII


def test_aav2_is_left_alone_by_some_puc19_single_cutters(examples, puc19):
    """The selection a cloning job actually needs: cuts vector, spares insert."""
    from restrict_cut.selection import choose

    aav2 = read_one(examples / "aav2.fasta")
    panel = enzymes.select(available_only=True, min_length=6)
    candidates = choose(
        puc19.sequence,
        panel,
        cuts=1,
        circular=True,
        protect=[("aav2", aav2.sequence, False)],
    ).candidates
    assert candidates
    for candidate in candidates:
        assert len(digest(puc19.sequence, [candidate.enzyme], circular=True).cuts) == 1
        assert digest(aav2.sequence, [candidate.enzyme]).cuts == []


def test_fragment_sums_hold_on_every_example(examples):
    """Whatever the enzyme or topology, the fragments must account for the whole."""
    panel = enzymes.resolve(["EcoRI", "BamHI", "HindIII", "MboI", "BsaI", "PvuII"])
    for filename in ["puc19.fasta", "pbr322.fasta", "aav2.fasta"]:
        record = read_one(examples / filename)
        for circular in (False, True):
            result = digest(record.sequence, panel, circular=circular)
            assert result.total_length == len(record), (filename, circular)
