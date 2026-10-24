"""A brute-force reference implementation, compared against the library.

The library finds sites with compiled character classes over a window that is
extended to let a site span the origin, and places cuts with offset arithmetic
that swaps the two nicks for a site read on the reverse strand. Each of those is
a place where an off-by-one or a sign error hides easily and reads plausibly.

So the same question is answered a second time by the slowest, most obvious code
that could work: try every start position, compare base by base, and place the
cuts from the definitions directly. The reference shares no code with the
library beyond the enzyme table itself, and the two must agree on every cut, for
every enzyme in the catalogue, on real plasmids and on sequence contrived to be
awkward, under both topologies.

This is what found the character-class bug: recognition patterns were compiled
into classes of the *bases* a site code stands for, so an ``N`` in the sequence
failed to match an ``N`` in the site, and SgeI (CNNG) and Hpy188I (TCNGA)
silently lost sites on ambiguous input.
"""

import random

import pytest

from restrict_cut import enzymes
from restrict_cut.digest import digest
from restrict_cut.reader import read_one

#: Spelled out rather than imported, so a change to the library's own table of
#: codes cannot quietly change what the reference believes.
BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S",
    "W": "W", "K": "M", "M": "K", "B": "V", "D": "H", "H": "D", "V": "B",
    "N": "N",
}


def reverse_complement(sequence):
    return "".join(COMPLEMENT[base] for base in reversed(sequence))


def reference_cuts(sequence, enzyme, circular):
    """Every cut, found the slow and obvious way."""
    length = len(sequence)
    found = set()
    if enzyme.searches_one_strand:
        readings = [("=", enzyme.pattern)]
    else:
        readings = [("+", enzyme.pattern), ("-", reverse_complement(enzyme.pattern))]

    for strand, pattern in readings:
        width = len(pattern)
        for start in range(1, length + 1):
            if not circular and start + width - 1 > length:
                continue
            if any(
                set(BASES[sequence[(start - 1 + step) % length]])
                - set(BASES[pattern[step]])
                for step in range(width)
            ):
                continue
            end = start + width - 1
            for top_offset, bottom_offset in enzyme.cuts:
                if strand == "-":
                    top, bottom = end - bottom_offset, end - top_offset
                else:
                    top, bottom = start + top_offset - 1, start + bottom_offset - 1
                if circular:
                    top = (top - 1) % length + 1
                    bottom = (bottom - 1) % length + 1
                elif not (1 <= top < length and 1 <= bottom < length):
                    continue
                found.add((top, bottom))
    return found


@pytest.fixture(scope="module")
def sequences(examples):
    rng = random.Random(20261024)
    cases = {
        "puc19": read_one(examples / "puc19.fasta").sequence,
        "pbr322": read_one(examples / "pbr322.fasta").sequence,
        "aav2": read_one(examples / "aav2.fasta").sequence,
        "random": "".join(rng.choice("ACGT") for _ in range(400)),
        "ambiguous": "".join(rng.choice("ACGTNRYSWKM") for _ in range(300)),
        "all_n": "N" * 80,
        "homopolymer": "A" * 120,
        "exact_site": "GAATTC",
    }
    return cases


@pytest.mark.parametrize("case", ["random", "ambiguous", "all_n", "homopolymer"])
@pytest.mark.parametrize("circular", [False, True])
def test_whole_catalogue_agrees_with_the_reference(sequences, case, circular):
    sequence = sequences[case]
    disagreements = []
    for enzyme in enzymes.all_enzymes().values():
        if enzyme.length > len(sequence):
            continue
        expected = reference_cuts(sequence, enzyme, circular)
        got = {(cut.top, cut.bottom) for cut in digest(sequence, [enzyme], circular).cuts}
        if expected != got:
            disagreements.append((enzyme.name, sorted(expected ^ got)[:4]))
    assert disagreements == []


@pytest.mark.parametrize("case", ["puc19", "pbr322", "aav2"])
@pytest.mark.parametrize("circular", [False, True])
def test_real_plasmids_agree_with_the_reference(sequences, case, circular):
    """A sample of the catalogue over real sequence, kept quick enough to run."""
    sequence = sequences[case]
    names = sorted(enzymes.all_enzymes())
    chosen = names[::7] + [
        "EcoRI", "BsaI", "SmaI", "XmaI", "HauII", "SgeI", "BcgI", "BsaXI",
        "PstI", "BstYI", "MboI", "SapI", "AlfI", "Hpy188I", "CspCI", "BaeI",
    ]
    disagreements = []
    for name in dict.fromkeys(chosen):
        enzyme = enzymes.get(name)
        expected = reference_cuts(sequence, enzyme, circular)
        got = {(cut.top, cut.bottom) for cut in digest(sequence, [enzyme], circular).cuts}
        if expected != got:
            disagreements.append((name, sorted(expected ^ got)[:4]))
    assert disagreements == []


def test_the_reference_disagrees_with_a_deliberately_broken_rule(sequences):
    """Confirm the comparison can fail, so a passing run means something.

    If the reference were insensitive to the matching rule, the tests above
    would pass whatever the library did. Compiling the site code's bases rather
    than the codes it admits is exactly the bug that was found; reproducing it
    here must change the answer on ambiguous sequence.
    """
    sequence = sequences["ambiguous"]
    enzyme = enzymes.get("SgeI")
    correct = reference_cuts(sequence, enzyme, False)

    def broken(sequence, enzyme):
        found = set()
        width = enzyme.length
        for start in range(1, len(sequence) - width + 2):
            window = sequence[start - 1: start - 1 + width]
            if all(base in BASES[code] for base, code in zip(window, enzyme.pattern)):
                for top_offset, bottom_offset in enzyme.cuts:
                    found.add((start + top_offset - 1, start + bottom_offset - 1))
        return found

    assert broken(sequence, enzyme) != correct
