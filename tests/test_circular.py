"""Properties that must hold however a circular molecule is written down.

Where the origin of a plasmid sequence falls is an arbitrary choice made by
whoever deposited it. Nothing a digest reports may depend on it, so these tests
rotate a real plasmid and require the answers to come back unchanged. They are
the sharpest check available on the wrapping arithmetic: an off-by-one or a
sign error in the origin handling survives a fixed test case far more easily
than it survives being asked the same question from eight different origins.
"""

from collections import Counter

import pytest

from restrict_cut import enzymes
from restrict_cut.digest import digest
from restrict_cut.reader import read_one

PANEL = [
    "EcoRI",
    "BamHI",
    "PstI",
    "SmaI",
    "BsaI",
    "HindIII",
    "MboI",
    "HauII",
    "BcgI",
    "KpnI",
    "SapI",
    "AatII",
]


@pytest.fixture(scope="module")
def puc19(examples):
    return read_one(examples / "puc19.fasta").sequence


def fingerprint(sequence, panel):
    """Everything about a digest that the choice of origin must not change."""
    result = digest(sequence, panel, circular=True)
    fragments = Counter(fragment.length for fragment in result.fragments)
    ends = Counter(
        (cut.enzyme, cut.end_type, cut.overhang(sequence, True)) for cut in result.cuts
    )
    return fragments, ends


@pytest.mark.parametrize("shift", [1, 7, 100, 395, 396, 397, 1343, 2685])
def test_digest_is_invariant_under_rotation(puc19, shift):
    panel = enzymes.resolve(PANEL)
    rotated = puc19[shift:] + puc19[:shift]
    assert fingerprint(rotated, panel) == fingerprint(puc19, panel)


@pytest.mark.parametrize("shift", [1, 396, 2685])
def test_cut_count_is_invariant_under_rotation(puc19, shift):
    panel = enzymes.resolve(PANEL)
    rotated = puc19[shift:] + puc19[:shift]
    assert digest(rotated, panel, circular=True).cut_counts == (
        digest(puc19, panel, circular=True).cut_counts
    )


def test_overhang_length_always_matches_the_signed_overhang(puc19):
    """Guards the origin-spanning overhang bug across the whole catalogue."""
    problems = []
    for enzyme in enzymes.all_enzymes().values():
        for circular in (False, True):
            result = digest(puc19, [enzyme], circular=circular)
            for cut in result.cuts:
                overhang = cut.overhang(puc19, circular)
                if len(overhang) != abs(cut.overhang_length):
                    problems.append((enzyme.name, circular, cut))
    assert problems == []


@pytest.mark.parametrize(
    "name,overhang",
    [
        ("EcoRI", "AATT"),
        ("BamHI", "GATC"),
        ("PstI", "TGCA"),
        ("KpnI", "GTAC"),
        ("HindIII", "AGCT"),
        ("SalI", "TCGA"),
        ("SphI", "CATG"),
        ("XbaI", "CTAG"),
    ],
)
def test_fixed_overhang_is_read_off_the_sequence(puc19, name, overhang):
    """What the site promises must be what the sequence yields, either topology."""
    for circular in (False, True):
        result = digest(puc19, [enzymes.get(name)], circular=circular)
        assert {cut.overhang(puc19, circular) for cut in result.cuts} == {overhang}


def test_a_cut_spanning_the_origin_still_gives_a_five_prime_overhang():
    """EcoRI nicks at 10 and 2 of a 12mer circle: four 5' bases, not eight 3'."""
    sequence = "TTCAAAAAAGAA"
    result = digest(sequence, [enzymes.get("EcoRI")], circular=True)
    cut = result.cuts[0]
    assert (cut.top, cut.bottom) == (10, 2)
    assert cut.overhang_length == 4
    assert cut.end_type == "5' overhang"
    assert cut.overhang(sequence, circular=True) == "AATT"


def test_circular_fragment_count_equals_cut_count(puc19):
    for name in ["EcoRI", "PvuII", "MboI", "HaeIII"]:
        result = digest(puc19, [enzymes.get(name)], circular=True)
        distinct = len({cut.top for cut in result.cuts})
        if distinct:
            assert len(result.fragments) == distinct
        assert result.total_length == len(puc19)


def test_linear_fragment_count_is_one_more_than_cut_count(puc19):
    for name in ["EcoRI", "PvuII", "MboI", "HaeIII"]:
        result = digest(puc19, [enzymes.get(name)], circular=False)
        distinct = len({cut.top for cut in result.cuts})
        assert len(result.fragments) == distinct + 1
        assert result.total_length == len(puc19)
