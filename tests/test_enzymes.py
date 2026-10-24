"""The generated catalogue and the properties derived from it."""

import pytest

from restrict_cut import enzymes
from restrict_cut.enzyme_data import ENZYMES
from restrict_cut.iupac import reverse_complement


def test_catalogue_is_internally_consistent():
    assert enzymes.validate_catalogue() == []


def test_catalogue_size_is_pinned():
    """A regenerated table should be noticed, not absorbed silently."""
    assert len(enzymes.all_enzymes()) == 720
    assert enzymes.REBASE_RELEASE == "610"


def test_most_enzymes_are_purchasable():
    available = [item for item in enzymes.all_enzymes().values() if item.is_available]
    assert len(available) == 581


def test_lookup_is_case_insensitive():
    assert enzymes.get("ecori").name == "EcoRI"
    assert enzymes.get("EcoRI").name == "EcoRI"
    assert enzymes.get(" BamHI ").name == "BamHI"


def test_unknown_enzyme_suggests_alternatives():
    with pytest.raises(KeyError, match="EcoRI"):
        enzymes.get("EcoR1")


def test_resolve_keeps_order_and_drops_repeats():
    chosen = enzymes.resolve(["PstI", "EcoRI", "psti"])
    assert [item.name for item in chosen] == ["PstI", "EcoRI"]


@pytest.mark.parametrize(
    "name,notation",
    [
        ("EcoRI", "G^AATTC"),
        ("PstI", "CTGCA^G"),
        ("SmaI", "CCC^GGG"),
        ("XmaI", "C^CCGGG"),
        ("BsaI", "GGTCTC(1/5)"),
        ("SapI", "GCTCTTC(1/4)"),
        ("BcgI", "(10/12)CGANNNNNNTGC(12/10)"),
        ("BstYI", "R^GATCY"),
    ],
)
def test_notation_matches_rebase(name, notation):
    assert enzymes.get(name).notation == notation


@pytest.mark.parametrize(
    "name,end_type,overhang",
    [
        ("EcoRI", "5' overhang", "AATT"),
        ("BamHI", "5' overhang", "GATC"),
        ("PstI", "3' overhang", "TGCA"),
        ("KpnI", "3' overhang", "GTAC"),
        ("SmaI", "blunt", ""),
        ("EcoRV", "blunt", ""),
    ],
)
def test_fixed_overhangs(name, end_type, overhang):
    enzyme = enzymes.get(name)
    assert enzyme.end_type() == end_type
    assert enzyme.fixed_overhang == overhang


def test_overhang_is_variable_when_the_enzyme_cuts_outside_its_site():
    assert enzymes.get("BsaI").fixed_overhang is None
    assert enzymes.get("BsaI").cuts_outside_site


def test_overhang_sign_convention():
    """bottom - top: positive is a 5' overhang, negative a 3' overhang."""
    assert enzymes.get("EcoRI").overhang_length() == 4
    assert enzymes.get("PstI").overhang_length() == -4
    assert enzymes.get("SmaI").overhang_length() == 0


def test_two_sided_enzymes_excise_their_site():
    assert enzymes.get("BcgI").excises_site
    assert len(enzymes.get("BcgI").cuts) == 2
    assert not enzymes.get("EcoRI").excises_site


def test_one_strand_search_requires_more_than_a_palindrome():
    """HauII reads the same both ways but cuts only downstream of its site."""
    hauii = enzymes.get("HauII")
    assert hauii.is_palindrome
    assert not hauii.searches_one_strand
    assert enzymes.get("EcoRI").searches_one_strand


def test_one_strand_flag_agrees_with_the_flip_it_describes():
    for enzyme in enzymes.all_enzymes().values():
        flipped = {
            (enzyme.length - bottom, enzyme.length - top)
            for top, bottom in enzyme.cuts
        }
        expected = enzyme.is_palindrome and flipped == set(enzyme.cuts)
        assert enzyme.searches_one_strand is expected, enzyme.name


def test_expected_frequency_accounts_for_ambiguity_and_strand():
    assert enzymes.get("EcoRI").expected_frequency == pytest.approx(4096)
    assert enzymes.get("NotI").expected_frequency == pytest.approx(65536)
    # BstYI's RGATCY has four fixed bases and two two-fold codes.
    assert enzymes.get("BstYI").expected_frequency == pytest.approx(1024)
    # A non-palindromic site is matched on both strands, so it comes up twice.
    assert enzymes.get("BsaI").expected_frequency == pytest.approx(2048)


def test_spacer_ns_do_not_add_specificity():
    """AhdI's site is eleven bases but only six of them are fixed."""
    ahdi = enzymes.get("AhdI")
    assert ahdi.length == 11
    assert ahdi.specific_bases == 6
    assert ahdi.expected_frequency == pytest.approx(enzymes.get("EcoRI").expected_frequency)


def test_isoschizomers_share_site_and_cut():
    smai = enzymes.get("SmaI")
    for name in enzymes.isoschizomers(smai):
        other = enzymes.get(name)
        assert other.pattern == smai.pattern
        assert other.cuts == smai.cuts


def test_neoschizomers_share_site_but_not_cut():
    smai = enzymes.get("SmaI")
    assert "XmaI" in enzymes.neoschizomers(smai)
    for name in enzymes.neoschizomers(smai):
        other = enzymes.get(name)
        assert other.pattern == smai.pattern
        assert other.cuts != smai.cuts


def test_an_enzyme_is_not_its_own_isoschizomer():
    for enzyme in list(enzymes.all_enzymes().values())[:50]:
        assert enzyme.name not in enzymes.isoschizomers(enzyme)
        assert enzyme.name not in enzymes.neoschizomers(enzyme)


def test_isoschizomers_are_symmetric():
    for name in ["SmaI", "EcoRI", "BsaI", "KpnI"]:
        enzyme = enzymes.get(name)
        for other in enzymes.isoschizomers(enzyme):
            assert name in enzymes.isoschizomers(enzymes.get(other))


def test_compatible_groups_find_the_classic_gatc_family():
    chosen = enzymes.select(available_only=True)
    groups = enzymes.compatible_groups(chosen)
    gatc = groups.get("5' overhang GATC", ())
    assert {"BamHI", "BglII", "BclI"} <= set(gatc)


def test_compatible_groups_separate_five_and_three_prime():
    """A 5'-GATC and a 3'-GATC expose the same bases but do not anneal."""
    groups = enzymes.compatible_groups(list(enzymes.all_enzymes().values()))
    for key, names in groups.items():
        if key == "blunt":
            continue
        end_type = "5' overhang" if key.startswith("5") else "3' overhang"
        for name in names:
            assert enzymes.get(name).end_type() == end_type


def test_pstI_and_sbfI_leave_compatible_ends():
    groups = enzymes.compatible_groups(list(enzymes.all_enzymes().values()))
    tgca = groups["3' overhang TGCA"]
    assert "PstI" in tgca and "SbfI" in tgca


def test_select_filters_combine():
    chosen = enzymes.select(available_only=True, min_length=8, end_type="blunt")
    assert chosen
    for enzyme in chosen:
        assert enzyme.is_available
        assert enzyme.length >= 8
        assert enzyme.end_type() == "blunt"


def test_select_rejects_unknown_filters():
    with pytest.raises(ValueError, match="end type"):
        enzymes.select(end_type="sticky")
    with pytest.raises(ValueError, match="supplier"):
        enzymes.select(supplier="ZZ")


def test_every_pattern_is_at_least_as_long_as_its_data():
    for name, (pattern, cuts, _) in ENZYMES.items():
        assert pattern
        assert cuts
        for top, bottom in cuts:
            assert isinstance(top, int) and isinstance(bottom, int), name


def test_palindromic_patterns_really_are_palindromes():
    for enzyme in enzymes.all_enzymes().values():
        if enzyme.is_palindrome:
            assert reverse_complement(enzyme.pattern) == enzyme.pattern, enzyme.name


def test_rarity_counts_distinct_sites_not_twice_one_strand():
    """Hin4I's site overlaps its own reverse complement.

    Doubling the one-strand count charges the overlapping sequences twice and
    makes the enzyme look about half as rare as it is.
    """
    hin4i = enzymes.get("Hin4I")
    assert hin4i.matching_sequences == 8192
    assert hin4i.expected_frequency == pytest.approx(512.0)


def test_rarity_agrees_with_counting_every_matching_sequence():
    """Check the closed form against brute force on short sites."""
    import itertools

    from restrict_cut.iupac import expand, reverse_complement

    for enzyme in enzymes.all_enzymes().values():
        if enzyme.length > 6 or enzyme.matching_sequences > 4096:
            continue
        forward = {
            "".join(combination)
            for combination in itertools.product(*[sorted(expand(b)) for b in enzyme.pattern])
        }
        flipped = enzyme.pattern if enzyme.is_palindrome else reverse_complement(enzyme.pattern)
        reverse = {
            "".join(combination)
            for combination in itertools.product(*[sorted(expand(b)) for b in flipped])
        }
        assert enzyme.matching_sequences == len(forward | reverse), enzyme.name


def test_rarity_of_the_familiar_enzymes_is_unchanged():
    assert enzymes.get("EcoRI").expected_frequency == pytest.approx(4096)
    assert enzymes.get("NotI").expected_frequency == pytest.approx(65536)
    assert enzymes.get("BstYI").expected_frequency == pytest.approx(1024)
    assert enzymes.get("BsaI").expected_frequency == pytest.approx(2048)


def test_a_non_palindromic_sticky_site_has_a_direction_dependent_overhang():
    """BauI leaves 5'-ACGA forward and 5'-TCGT reverse; they do not anneal."""
    baui = enzymes.get("BauI")
    assert not baui.is_palindrome
    assert baui.fixed_overhang == "ACGA"
    assert baui.overhang_depends_on_orientation


def test_palindromic_and_blunt_enzymes_do_not_depend_on_orientation():
    for name in ["EcoRI", "BamHI", "SmaI", "PstI"]:
        assert not enzymes.get(name).overhang_depends_on_orientation


def test_static_groups_exclude_direction_dependent_overhangs():
    """Grouping BauI with BsiI promises a pairing half the sites will not make."""
    groups = enzymes.compatible_groups(list(enzymes.all_enzymes().values()))
    assert not any({"BauI", "BsiI"} <= set(names) for names in groups.values())
    for names in groups.values():
        for name in names:
            assert not enzymes.get(name).overhang_depends_on_orientation, name


def test_static_groups_still_hold_the_palindromic_families():
    groups = enzymes.compatible_groups(enzymes.select(available_only=True))
    assert {"BamHI", "BglII", "BclI"} <= set(groups["5' overhang GATC"])
    assert {"PstI", "SbfI"} <= set(groups["3' overhang TGCA"])
