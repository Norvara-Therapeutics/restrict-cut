"""Dam, Dcm and CpG overlap with recognition sites."""

import pytest

from restrict_cut import enzymes
from restrict_cut.methylation import (
    DEFAULT_SYSTEMS,
    SYSTEMS,
    context,
    exposure,
    find_targets,
    resolve_systems,
    site_overlaps,
)
from restrict_cut.search import find_sites

DAM = SYSTEMS["dam"]
DCM = SYSTEMS["dcm"]
CPG = SYSTEMS["cpg"]


def test_default_systems_are_the_e_coli_ones():
    assert DEFAULT_SYSTEMS == ("dam", "dcm")


def test_system_labels_mark_the_modified_base():
    assert DAM.label == "G<A>TC"
    assert DCM.label == "C<C>WGG"


@pytest.mark.parametrize(
    "name,expected",
    [
        ("MboI", "always"),
        ("DpnI", "always"),
        ("BamHI", "always"),
        ("BclI", "always"),
        ("BstYI", "always"),
        ("XbaI", "possible"),
        ("ClaI", "possible"),
        ("TaqI", "possible"),
        ("EcoRI", "never"),
        ("SmaI", "never"),
        ("NheI", "never"),
    ],
)
def test_dam_exposure_classes(name, expected):
    assert exposure(enzymes.get(name), DAM) == expected


@pytest.mark.parametrize(
    "name,expected",
    [
        ("EcoRII", "always"),
        ("PspGI", "always"),
        ("StuI", "possible"),
        ("EcoRI", "never"),
        ("MboI", "never"),
    ],
)
def test_dcm_exposure_classes(name, expected):
    assert exposure(enzymes.get(name), DCM) == expected


def test_always_means_every_site_really_does_carry_the_target():
    """GGATCC contains GATC whatever surrounds it."""
    sequence = "AAAAGGATCCAAAA"
    enzyme = enzymes.get("BamHI")
    overlaps = site_overlaps(
        sequence, enzyme, find_sites(sequence, enzyme), [DAM]
    )
    # Both strands of the GATC are marked, so both modified bases are reported.
    assert len(overlaps) == 2
    assert {overlap.system for overlap in overlaps} == {"dam"}
    assert {overlap.position for overlap in overlaps} == {7, 8}


def test_never_means_no_overlap_is_ever_reported():
    sequence = "GATCGAATTCGATC"
    enzyme = enzymes.get("EcoRI")
    assert exposure(enzyme, DAM) == "never"
    assert site_overlaps(sequence, enzyme, find_sites(sequence, enzyme), [DAM]) == []


def test_possible_depends_on_the_context():
    """ClaI's ATCGAT overlaps a Dam GATC only when the next base is a C."""
    enzyme = enzymes.get("ClaI")
    with_target = "AAAATCGATCAAA"
    without = "AAAATCGATAAAA"
    assert site_overlaps(
        with_target, enzyme, find_sites(with_target, enzyme), [DAM]
    )
    assert (
        site_overlaps(without, enzyme, find_sites(without, enzyme), [DAM]) == []
    )


def test_find_targets_reports_the_modified_base_on_each_strand():
    """Dam marks the adenine of GATC on both strands, at two forward positions."""
    assert find_targets("GGATCCAAA", DAM) == [(2, 3), (2, 4)]


def test_cpg_is_modified_on_both_strands():
    """CG is its own reverse complement, but the modified C is not centred."""
    assert find_targets("ACGT", CPG) == [(2, 2), (2, 3)]


def test_dcm_marks_the_inner_cytosine_of_both_strands():
    assert find_targets("CCAGG", DCM) == [(1, 2), (1, 4)]


def test_find_targets_wraps_when_circular():
    assert find_targets("TCAAAAAAAGA", DAM) == []
    assert find_targets("TCAAAAAAAGA", DAM, circular=True) == [(10, 1), (10, 11)]


def test_dcm_target_accepts_both_ambiguity_resolutions():
    assert find_targets("CCAGG", DCM)
    assert find_targets("CCTGG", DCM)
    assert find_targets("CCCGG", DCM) == []


def test_resolve_systems_defaults_and_names():
    assert [system.name for system in resolve_systems(None)] == list(DEFAULT_SYSTEMS)
    assert [system.name for system in resolve_systems(["cpg"])] == ["cpg"]
    assert resolve_systems([]) == []


def test_resolve_systems_rejects_unknown():
    with pytest.raises(KeyError, match="unknown methylation system"):
        resolve_systems(["ecoki"])


def test_context_window_wraps_on_a_circle():
    assert context("ACGTACGT", 1, 2, circular=True) == "GTACG"
    assert context("ACGTACGT", 1, 2, circular=False) == "ACG"


def test_exposure_never_is_consistent_with_an_exhaustive_check():
    """A 'never' enzyme must have no overlap on heavily targeted sequence."""
    sequence = "GATCCCAGGGATCCCTGG" * 6
    for enzyme in list(enzymes.all_enzymes().values())[:120]:
        for system in (DAM, DCM):
            if exposure(enzyme, system) != "never":
                continue
            found = site_overlaps(
                sequence, enzyme, find_sites(sequence, enzyme), [system]
            )
            assert found == [], (enzyme.name, system.name)


def test_exposure_only_ever_returns_a_documented_class():
    """EXPOSURE is the documented domain of exposure(); hold it to that."""
    from restrict_cut.methylation import EXPOSURE

    for enzyme in enzymes.all_enzymes().values():
        for system in SYSTEMS.values():
            assert exposure(enzyme, system) in EXPOSURE, (enzyme.name, system.name)


def _loop_targets(sequence, system, circular=False):
    """The per-position loop the regex scan replaced, kept as a reference."""
    from restrict_cut.iupac import matches

    length, width = len(sequence), len(system.target)
    if not length or width > length:
        return []
    reach = min(width - 1, length) if circular else 0
    window = sequence + sequence[:reach]
    found = set()
    for pattern, offset in system.readings:
        for index in range(len(window) - width + 1):
            if index >= length:
                break
            if all(matches(pattern[step], window[index + step]) for step in range(width)):
                found.add((index % length + 1, (index + offset) % length + 1))
    return sorted(found)


@pytest.mark.parametrize("system", list(SYSTEMS.values()), ids=list(SYSTEMS))
@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize(
    "sequence",
    ["GATC" * 20, "CCAGGCCTGG" * 12, "N" * 60, "CG", "AAAAGATC", "ACGTRYSWKMNNGATCCCAGG"],
)
def test_regex_target_scan_matches_the_loop_it_replaced(system, circular, sequence):
    assert find_targets(sequence, system, circular) == _loop_targets(
        sequence, system, circular
    )


def test_memoised_targets_do_not_leak_between_topologies():
    """The cache is keyed on topology as well as sequence."""
    sequence = "TCAAAAAAAGA"
    linear = find_targets(sequence, DAM, circular=False)
    circular = find_targets(sequence, DAM, circular=True)
    assert linear != circular
    assert find_targets(sequence, DAM, circular=False) == linear


def test_context_window_never_shrinks_as_it_widens():
    """A circular window wider than the molecule once wrapped past its own start.

    Asking a five-base circle for seven bases of context returned two, because
    the slice wrapped the origin only once.
    """
    sequence = "GATCG"
    widths = [len(context(sequence, 3, flank, True)) for flank in range(0, 12)]
    assert widths == sorted(widths)
    assert max(widths) == len(sequence)


def test_context_clamps_to_the_whole_circle():
    assert context("GATCG", 3, 10, circular=True) == "GATCG"
    assert len(context("GATCG", 3, 3, circular=True)) == 5


def test_context_centres_the_position_on_a_small_circle():
    window = context("GATCG", 3, 10, circular=True)
    assert len(window) == 5
    assert set(window) == set("GATCG")


def test_context_on_a_linear_sequence_clamps_to_the_ends():
    assert context("GATCG", 1, 10, circular=False) == "GATCG"
    assert context("GATCG", 1, 1, circular=False) == "GA"
