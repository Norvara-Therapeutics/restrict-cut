"""Cut placement and fragment assembly."""

import pytest

from restrict_cut import enzymes
from restrict_cut.digest import digest, slice_sequence, span_length


def cut_of(sequence, name, circular=False):
    return digest(sequence, [enzymes.get(name)], circular=circular)


def test_nick_positions_of_a_known_cut():
    """EcoRI G^AATTC at site 4..9 nicks after base 4 and after base 8."""
    result = cut_of("AAAGAATTCAAA", "EcoRI")
    assert len(result.cuts) == 1
    cut = result.cuts[0]
    assert (cut.top, cut.bottom) == (4, 8)
    assert cut.overhang("AAAGAATTCAAA") == "AATT"
    assert cut.end_type == "5' overhang"


def test_three_prime_overhang_has_the_nicks_the_other_way_round():
    sequence = "AAACTGCAGAAA"
    result = cut_of(sequence, "PstI")
    cut = result.cuts[0]
    assert cut.bottom < cut.top
    assert cut.end_type == "3' overhang"
    assert cut.overhang(sequence) == "TGCA"


def test_blunt_cut_has_equal_nicks():
    result = cut_of("AAACCCGGGAAA", "SmaI")
    cut = result.cuts[0]
    assert cut.top == cut.bottom
    assert cut.end_type == "blunt"
    assert cut.overhang("AAACCCGGGAAA") == ""


def test_reverse_strand_site_gives_the_same_overhang_type():
    """BsaI on either strand must still leave a 5' overhang."""
    forward = cut_of("AAGGTCTCAAAAAAAA", "BsaI")
    reverse = cut_of("AAAAAAAAGAGACCTT", "BsaI")
    assert forward.cuts[0].end_type == "5' overhang"
    assert reverse.cuts[0].end_type == "5' overhang"
    assert forward.cuts[0].strand == "+"
    assert reverse.cuts[0].strand == "-"


def test_linear_n_cuts_give_n_plus_one_fragments():
    result = digest(
        "GAATTCAAAAAAAAAAGAATTCAAAAAAAAAA", [enzymes.get("EcoRI")], circular=False
    )
    assert len(result.cuts) == 2
    assert len(result.fragments) == 3


def test_circular_n_cuts_give_n_fragments():
    result = digest(
        "GAATTCAAAAAAAAAAGAATTCAAAAAAAAAA", [enzymes.get("EcoRI")], circular=True
    )
    assert len(result.cuts) == 2
    assert len(result.fragments) == 2


def test_one_cut_linearises_a_circle():
    result = cut_of("AAAGAATTCAAA", "EcoRI", circular=True)
    assert len(result.fragments) == 1
    assert result.fragments[0].length == 12
    assert result.fragments[0].spans_origin


def test_uncut_circle_is_reported_as_one_circular_fragment():
    result = cut_of("AAAAAAAAAAAA", "EcoRI", circular=True)
    assert result.is_uncut
    assert len(result.fragments) == 1
    assert result.fragments[0].circular


def test_uncut_linear_sequence_is_one_fragment():
    result = cut_of("AAAAAAAAAAAA", "EcoRI")
    assert len(result.fragments) == 1
    assert not result.fragments[0].circular


def test_fragment_lengths_sum_to_the_molecule():
    sequence = "GAATTCAAAGGATCCAAACTGCAGAAACCCGGGAAAGAATTCAAA"
    for circular in (False, True):
        result = digest(
            sequence,
            enzymes.resolve(["EcoRI", "BamHI", "PstI", "SmaI"]),
            circular=circular,
        )
        assert result.total_length == len(sequence)


def test_fragments_are_sorted_longest_first():
    sequence = "GAATTCAAAAAAAAAAAAAAAAAAAAGAATTCAAA"
    result = digest(sequence, [enzymes.get("EcoRI")])
    lengths = [fragment.length for fragment in result.fragments]
    assert lengths == sorted(lengths, reverse=True)


def test_linear_ends_are_labelled_as_termini():
    result = cut_of("AAAGAATTCAAA", "EcoRI")
    assert result.fragments[0].left in {"end", "EcoRI"}
    assert {fragment.left for fragment in result.fragments} == {"end", "EcoRI"}


def test_a_type_iis_site_too_near_a_linear_end_does_not_cut():
    """BsaI cuts seven bases downstream, so a site at the very end cannot."""
    sequence = "AAAAAAAAGGTCTC"
    result = cut_of(sequence, "BsaI")
    assert result.sites_for("BsaI")
    assert result.cuts == []
    assert result.uncut_sites["BsaI"] == 1
    assert result.partial_sites == {}


def test_the_same_site_cuts_when_the_molecule_is_circular():
    sequence = "AAAAAAAAGGTCTC"
    result = cut_of(sequence, "BsaI", circular=True)
    assert result.cuts
    assert not result.uncut_sites


def test_a_two_sided_site_near_a_linear_end_is_cut_but_not_excised():
    """One of BcgI's two breaks lands and the other does not.

    Counting that site as uncut would contradict the cut reported from it, so
    the two cases are kept apart: no break at all is an uncut site, fewer breaks
    than the enzyme makes is a partial one.
    """
    sequence = "CGAAAAAAATGC" + "A" * 40
    result = cut_of(sequence, "BcgI")
    assert len(result.sites_for("BcgI")) == 1
    assert len(result.cuts) == 1
    assert result.partial_sites == {"BcgI": 1}
    assert result.uncut_sites == {}


def test_two_sided_enzyme_makes_two_cuts_per_site():
    sequence = "A" * 30 + "CGAAAAAAATGC" + "A" * 30
    result = cut_of(sequence, "BcgI")
    assert len(result.sites_for("BcgI")) == 1
    assert len(result.cuts) == 2
    assert len(result.fragments) == 3


def test_a_double_digest_combines_both_enzymes():
    sequence = "GAATTCAAAAAAAAAAGGATCCAAAAAAAAAA"
    single = digest(sequence, [enzymes.get("EcoRI")])
    double = digest(sequence, enzymes.resolve(["EcoRI", "BamHI"]))
    assert len(double.cuts) == len(single.cuts) + 1
    assert double.cut_counts == {"EcoRI": 1, "BamHI": 1}


def test_cut_counts_include_enzymes_that_do_not_cut():
    result = digest("AAAAAAAAAAAA", enzymes.resolve(["EcoRI", "BamHI"]))
    assert result.cut_counts == {"EcoRI": 0, "BamHI": 0}


def test_identical_breaks_are_not_counted_twice():
    """A symmetric site reached from both orientations is one physical cut."""
    sequence = "AAAA" + "TGGCCA" + "AAAA"
    result = cut_of(sequence, "HauII", circular=True)
    pairs = {(cut.top, cut.bottom) for cut in result.cuts}
    assert len(pairs) == len(result.cuts)


def test_overhang_wraps_across_the_origin():
    sequence = "TTCAAAAAAGAA"
    result = cut_of(sequence, "EcoRI", circular=True)
    assert result.cuts[0].overhang(sequence, circular=True) == "AATT"


def test_slice_sequence_wraps_only_when_circular():
    assert slice_sequence("ABCDEF", 2, 4) == "BCD"
    assert slice_sequence("ABCDEF", 5, 2, circular=True) == "EFAB"
    with pytest.raises(ValueError):
        slice_sequence("ABCDEF", 5, 2)


def test_span_length_handles_wrapping():
    assert span_length(2, 4, 10) == 3
    assert span_length(9, 2, 10) == 4


def test_empty_sequence_yields_nothing():
    result = digest("", [enzymes.get("EcoRI")])
    assert result.fragments == []
    assert result.cuts == []


def test_topology_is_named():
    assert digest("AAAA", [enzymes.get("EcoRI")], circular=True).topology == "circular"
    assert digest("AAAA", [enzymes.get("EcoRI")]).topology == "linear"


def test_two_enzymes_breaking_one_bond_share_a_boundary():
    """AccII and Alw26I both cut after base 3 of pUC19, differently.

    One top-strand bond is broken, so it is one fragment boundary and the
    fragment count follows the boundaries rather than the cuts. Crediting only
    the first enzyme would drop the other from the fragment it helped make.
    """
    import pathlib

    from restrict_cut.reader import read_one

    examples = pathlib.Path(__file__).resolve().parent.parent / "examples"
    sequence = read_one(examples / "puc19.fasta").sequence
    result = digest(sequence, enzymes.resolve(["AccII", "Alw26I"]), circular=True)

    assert len(result.cuts) == 14
    assert result.boundaries == 13
    assert len(result.fragments) == 13
    assert result.shared_boundaries == {3: ("AccII", "Alw26I")}
    assert result.total_length == len(sequence)

    labels = {fragment.left for fragment in result.fragments}
    labels |= {fragment.right for fragment in result.fragments}
    assert "AccII+Alw26I" in labels


def test_a_shared_boundary_names_every_enzyme_once():
    sequence = "AAAACGCGAAAA"
    result = digest(sequence, enzymes.resolve(["AccII"]), circular=True)
    for fragment in result.fragments:
        for label in (fragment.left, fragment.right):
            names = label.split("+")
            assert len(names) == len(set(names))


def test_fragment_count_follows_boundaries_not_cuts():
    """The invariant that makes 'more cuts than fragments' consistent."""
    import pathlib

    from restrict_cut.reader import read_one

    examples = pathlib.Path(__file__).resolve().parent.parent / "examples"
    sequence = read_one(examples / "puc19.fasta").sequence
    panel = enzymes.resolve(["AccII", "Alw26I", "MboI", "HaeIII", "EcoRI"])
    for circular in (False, True):
        result = digest(sequence, panel, circular=circular)
        expected = result.boundaries if circular else result.boundaries + 1
        assert len(result.fragments) == (expected if result.boundaries else 1)
        assert result.total_length == len(sequence)


def test_overhang_longer_than_the_molecule_still_agrees_with_its_length():
    """A five-base circle cut by TscAI wraps the origin almost twice.

    The molecule is nonsense, but a nine-base overhang must not be reported
    with a four-base sequence.
    """
    sequence = "CACTG"
    result = digest(sequence, [enzymes.get("TscAI")], circular=True)
    for cut in result.cuts:
        overhang = cut.overhang(sequence, circular=True)
        assert len(overhang) == abs(cut.overhang_length)


def test_overhangs_agree_with_their_length_on_every_small_circle():
    """Swept across short circles, where wrapping arithmetic goes wrong."""
    for size in range(4, 20):
        sequence = ("CACTGGATCCGAATTC" * 3)[:size]
        for enzyme in enzymes.resolve(["TscAI", "BcgI", "BsaI", "EcoRI", "HauII"]):
            result = digest(sequence, [enzyme], circular=True)
            for cut in result.cuts:
                assert len(cut.overhang(sequence, True)) == abs(cut.overhang_length)
