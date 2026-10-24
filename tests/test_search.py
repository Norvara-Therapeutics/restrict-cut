"""Finding recognition sites, including across the origin and on both strands."""

import pytest

from restrict_cut import enzymes
from restrict_cut.search import count_uncertain, find_sites


def sites(sequence, name, circular=False):
    return find_sites(sequence, enzymes.get(name), circular)


def test_single_site_is_found_with_its_span():
    found = sites("AAAGAATTCAAA", "EcoRI")
    assert len(found) == 1
    assert (found[0].start, found[0].end) == (4, 9)


def test_palindromic_site_is_reported_once():
    found = sites("GAATTCGAATTC", "EcoRI")
    assert [site.start for site in found] == [1, 7]
    assert all(site.strand == "=" for site in found)


def test_non_palindromic_site_is_found_on_both_strands():
    found = sites("AAGGTCTCAAAAAGAGACCTT", "BsaI")
    assert [(site.start, site.strand) for site in found] == [(3, "+"), (14, "-")]


def test_site_spanning_the_origin_needs_circular():
    assert sites("TTCAAAAAAGAA", "EcoRI") == []
    found = sites("TTCAAAAAAGAA", "EcoRI", circular=True)
    assert len(found) == 1
    assert (found[0].start, found[0].end) == (10, 3)
    assert found[0].spans_origin


def test_circular_search_does_not_duplicate_interior_sites():
    found = sites("AAAGAATTCAAA", "EcoRI", circular=True)
    assert len(found) == 1


def test_overlapping_sites_are_all_found():
    """CNNG tiles along a sequence; a non-overlapping scan would miss some."""
    found = sites("CAAGAAGAAG", "SgeI")
    assert len(found) >= 2


def test_ambiguity_in_the_site_widens_what_matches():
    assert len(sites("AGATCT", "BstYI")) == 1
    assert len(sites("GGATCC", "BstYI")) == 1
    assert sites("CGATCG", "BstYI") == []


def test_ambiguity_in_the_sequence_is_not_called_a_site():
    assert sites("AAAGANTTCAAA", "EcoRI") == []


def test_ambiguity_in_the_sequence_is_counted_as_uncertain():
    assert count_uncertain("AAAGANTTCAAA", enzymes.get("EcoRI")) == 1


def test_certain_sites_are_not_counted_as_uncertain():
    assert count_uncertain("AAAGAATTCAAA", enzymes.get("EcoRI")) == 0


def test_an_n_in_the_sequence_matches_an_n_in_the_site():
    """Guards the character-class bug: site N must admit sequence N."""
    found = sites("AACNNGAA", "SgeI")
    assert found


def test_no_sites_in_an_empty_or_short_sequence():
    assert sites("", "EcoRI") == []
    assert sites("GAATT", "EcoRI") == []


def test_sites_are_sorted_by_position():
    found = sites("GAATTCAAAAAGAATTC", "EcoRI")
    assert [site.start for site in found] == sorted(site.start for site in found)


def test_lowercase_sequence_is_matched():
    assert len(sites("aaagaattcaaa", "EcoRI")) == 1


@pytest.mark.parametrize("name", ["EcoRI", "BsaI", "BcgI", "HauII", "SgeI", "MboI"])
def test_circular_never_finds_fewer_sites_than_linear(name):
    sequence = "GGATCCGAATTCGGTCTCTGGCCACAAGCGANNNNNNTGCAAAGATC"
    linear = sites(sequence, name)
    circular = sites(sequence, name, circular=True)
    assert len(circular) >= len(linear)


def test_strand_labels_are_always_from_the_documented_set():
    """STRANDS is the documented domain of Site.strand; hold it to that."""
    from restrict_cut.search import STRANDS

    sequence = "GGATCCGAATTCGGTCTCTGGCCACAAGCGANNNNNNTGCAAAGATCCASTGCNNG"
    for enzyme in enzymes.all_enzymes().values():
        for circular in (False, True):
            for site in find_sites(sequence, enzyme, circular):
                assert site.strand in STRANDS, (enzyme.name, site)


def test_one_strand_enzymes_are_labelled_as_such():
    """A site read identically both ways is '=', never '+' or '-'."""
    sequence = "GGATCCGAATTCAAGCTTCTGCAG"
    for enzyme in enzymes.all_enzymes().values():
        if not enzyme.searches_one_strand:
            continue
        for site in find_sites(sequence, enzyme):
            assert site.strand == "="
