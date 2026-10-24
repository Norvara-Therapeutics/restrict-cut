"""Choosing enzymes by how often they cut."""

import pytest

from restrict_cut import enzymes
from restrict_cut.digest import digest
from restrict_cut.reader import read_one
from restrict_cut.selection import choose, compatible_with, survey


@pytest.fixture(scope="module")
def puc19(examples):
    return read_one(examples / "puc19.fasta").sequence


@pytest.fixture(scope="module")
def panel():
    return enzymes.select(available_only=True, min_length=6)


def test_survey_counts_agree_with_a_digest(puc19, panel):
    found = survey(puc19, panel[:40], circular=True)
    for enzyme in panel[:40]:
        result = digest(puc19, [enzyme], circular=True)
        assert found.counts[enzyme.name] == len(result.cuts), enzyme.name


def test_survey_positions_are_sorted(puc19, panel):
    found = survey(puc19, panel[:40], circular=True)
    for positions in found.positions.values():
        assert list(positions) == sorted(positions)


def test_cutting_selects_by_exact_count(puc19, panel):
    found = survey(puc19, panel, circular=True)
    for name in found.cutting(1):
        assert found.counts[name] == 1
    for name in found.cutting(0):
        assert found.counts[name] == 0


def test_choose_returns_only_matching_enzymes(puc19, panel):
    candidates = choose(puc19, panel, cuts=1, circular=True).candidates
    assert candidates
    for candidate in candidates:
        assert candidate.cuts == 1
        assert len(candidate.positions) == 1


def test_choose_finds_the_polylinker_single_cutters(puc19, panel):
    """Every unique site of the pUC19 polylinker must come back."""
    names = {candidate.enzyme.name for candidate in choose(puc19, panel, 1, circular=True).candidates}
    assert {
        "EcoRI",
        "SacI",
        "KpnI",
        "SmaI",
        "BamHI",
        "XbaI",
        "SalI",
        "PstI",
        "SphI",
        "HindIII",
    } <= names


def test_choose_orders_rarest_and_most_available_first(puc19, panel):
    candidates = choose(puc19, panel, cuts=1, circular=True).candidates
    frequencies = [candidate.enzyme.expected_frequency for candidate in candidates]
    assert frequencies == sorted(frequencies, reverse=True)


def test_choose_zero_cuts_finds_non_cutters(puc19, panel):
    candidates = choose(puc19, panel, cuts=0, circular=True).candidates
    assert candidates
    for candidate in candidates:
        assert candidate.cuts == 0
        assert digest(puc19, [candidate.enzyme], circular=True).cuts == []


def test_protect_excludes_enzymes_cutting_the_other_sequence(puc19, panel, examples):
    aav2 = read_one(examples / "aav2.fasta").sequence
    unprotected = {
        candidate.enzyme.name for candidate in choose(puc19, panel, 1, circular=True).candidates
    }
    protected = {
        candidate.enzyme.name
        for candidate in choose(
            puc19, panel, 1, circular=True, protect=[("aav2", aav2, False)]
        ).candidates
    }
    assert protected < unprotected
    for name in protected:
        assert digest(aav2, [enzymes.get(name)]).cuts == []


def test_max_cuts_accepts_a_range(puc19, panel):
    candidates = choose(puc19, panel, cuts=None, max_cuts=2, circular=True).candidates
    assert candidates
    assert {candidate.cuts for candidate in candidates} <= {0, 1, 2}


def test_choose_rejects_negative_counts(puc19, panel):
    with pytest.raises(ValueError):
        choose(puc19, panel[:5], cuts=-1).candidates
    with pytest.raises(ValueError):
        choose(puc19, panel[:5], max_cuts=-1).candidates


def test_uncertain_positions_are_flagged():
    sequence = "GAATTC" + "A" * 40 + "GANTTC" + "A" * 40
    candidates = choose(sequence, enzymes.resolve(["EcoRI"]), cuts=1).candidates
    assert candidates
    assert candidates[0].uncertain == 1
    assert not candidates[0].is_clean


def test_compatible_with_matches_gatc_overhangs():
    sequence = "AAAGGATCCAAAAAAAGATCTAAA"
    bam = digest(sequence, [enzymes.get("BamHI")]).cuts[0]
    bgl = digest(sequence, [enzymes.get("BglII")]).cuts[0]
    assert compatible_with(sequence, bam, [(bgl, "BglII")]) == ["BglII"]


def test_compatible_with_separates_end_types():
    sequence = "AAACTGCAGAAAAGGATCCAAA"
    pst = digest(sequence, [enzymes.get("PstI")]).cuts[0]
    bam = digest(sequence, [enzymes.get("BamHI")]).cuts[0]
    assert compatible_with(sequence, pst, [(bam, "BamHI")]) == []


def test_choose_returns_a_selection_describing_what_it_searched(puc19, panel):
    """The caller needs the exclusions, not only the survivors."""
    from restrict_cut.selection import Selection

    result = choose(puc19, panel, cuts=1, circular=True, name="pUC19")
    assert isinstance(result, Selection)
    assert result.target == "pUC19"
    assert result.considered == len(panel)
    assert result.protected == ()
    assert result.excluded_uncertain == ()


def test_protected_uncertainty_excludes_rather_than_approves():
    """No confirmed site is not the same as no site."""
    enzymes_used = enzymes.resolve(["EcoRI"])
    result = choose(
        "GAATTC", enzymes_used, cuts=1, protect=[("guard", "GANTTC", False)]
    )
    assert result.candidates == []
    assert result.excluded_uncertain == ("EcoRI",)


def test_protected_confirmed_cut_is_not_reported_as_uncertain():
    """A real site excludes on its own account, not through the uncertain list."""
    result = choose(
        "GAATTC", enzymes.resolve(["EcoRI"]), cuts=1,
        protect=[("guard", "GAATTC", False)],
    )
    assert result.candidates == []
    assert result.excluded_uncertain == ()


def test_protection_names_the_sequences_it_checked():
    result = choose(
        "GAATTC", enzymes.resolve(["EcoRI"]), cuts=1,
        protect=[("insert", "TTTTTT", False)],
    )
    assert result.protected == ("insert",)
    assert [item.enzyme.name for item in result.candidates] == ["EcoRI"]


def test_uncertainty_in_the_target_does_not_exclude(puc19):
    """Target uncertainty is reported; only protected uncertainty excludes."""
    sequence = "GAATTC" + "A" * 30 + "GANTTC" + "A" * 30
    result = choose(sequence, enzymes.resolve(["EcoRI"]), cuts=1)
    assert [item.enzyme.name for item in result.candidates] == ["EcoRI"]
    assert result.candidates[0].uncertain == 1
    assert not result.candidates[0].is_clean


def test_end_signature_is_read_from_the_molecule_the_cut_is_in():
    """Comparing a cut against a sequence it is not in invents a pairing.

    BauI at a forward site and BsiI at a reverse site leave complementary
    overhangs, which do not anneal. Reading the second cut's overhang off the
    first molecule returns the first molecule's bases and reports a match.
    """
    from restrict_cut.iupac import reverse_complement
    from restrict_cut.selection import end_signature, ends_compatible

    baui, bsii = enzymes.get("BauI"), enzymes.get("BsiI")
    forward = "AAAA" + baui.pattern + "AAAA"
    reverse = "AAAA" + reverse_complement(bsii.pattern) + "AAAA"
    first = digest(forward, [baui]).cuts[0]
    second = digest(reverse, [bsii]).cuts[0]

    assert end_signature(forward, first) == ("5' overhang", "ACGA")
    assert end_signature(reverse, second) == ("5' overhang", "TCGT")
    assert not ends_compatible(
        end_signature(forward, first), end_signature(reverse, second)
    )


def test_ends_compatible_matches_the_classic_families():
    from restrict_cut.selection import end_signature, ends_compatible

    sequence = "AAAGGATCCAAAAAAAGATCTAAA"
    bam = digest(sequence, [enzymes.get("BamHI")]).cuts[0]
    bgl = digest(sequence, [enzymes.get("BglII")]).cuts[0]
    assert ends_compatible(
        end_signature(sequence, bam), end_signature(sequence, bgl)
    )


def test_ends_compatible_separates_the_two_strands():
    from restrict_cut.selection import ends_compatible

    assert not ends_compatible(("5' overhang", "GATC"), ("3' overhang", "GATC"))


def test_survey_records_the_end_observed_not_the_one_assumed(puc19):
    """BspQI cuts outside its site, so only the sequence gives its overhang."""
    found = survey(puc19, enzymes.resolve(["BspQI"]), circular=True)
    assert found.overhangs["BspQI"] == (("5' overhang", "GCT"),)
    assert enzymes.get("BspQI").fixed_overhang is None


def test_candidates_carry_the_observed_ends(puc19, panel):
    candidates = choose(puc19, panel, cuts=1, circular=True).candidates
    for candidate in candidates:
        assert len(candidate.overhangs) >= 1
        for end_type, overhang in candidate.overhangs:
            assert end_type in ("5' overhang", "blunt", "3' overhang")
            if end_type == "blunt":
                assert overhang == ""


def test_case_does_not_make_two_compatible_ends_incompatible():
    """Site matching ignores case, so end comparison must too.

    A soft-masked region is written in lower case; BamHI's gatc and BglII's
    GATC are the same four bases and will ligate.
    """
    from restrict_cut.selection import end_signature, ends_compatible

    sequence = "AAAggatccAAAAAAAGATCTAAA"
    bam = digest(sequence, [enzymes.get("BamHI")]).cuts[0]
    bgl = digest(sequence, [enzymes.get("BglII")]).cuts[0]
    assert end_signature(sequence, bam) == ("5' overhang", "GATC")
    assert end_signature(sequence, bgl) == ("5' overhang", "GATC")
    assert ends_compatible(
        end_signature(sequence, bam), end_signature(sequence, bgl)
    )
    assert compatible_with(sequence, bam, [(bgl, "BglII")]) == ["BglII"]


def test_normalising_case_does_not_merge_the_two_strands():
    from restrict_cut.selection import ends_compatible

    assert not ends_compatible(("5' overhang", "GATC"), ("3' overhang", "GATC"))


def test_a_fully_lowercase_sequence_surveys_the_same_as_upper():
    upper = "AAAGGATCCAAAAAAAGATCTAAA"
    lower = upper.lower()
    panel = enzymes.resolve(["BamHI", "BglII"])
    assert survey(upper, panel).overhangs == survey(lower, panel).overhangs
