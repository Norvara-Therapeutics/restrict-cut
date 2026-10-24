"""IUPAC nucleotide codes, matching and reverse complement.

A recognition site is written in IUPAC codes, and so, sometimes, is the DNA
being searched. The two meet under a deliberately asymmetric rule, set out in
:func:`matches`: an ambiguous base in the *site* widens what the enzyme accepts,
while an ambiguous base in the *sequence* means the sequence is not known, and
cannot be reported as a site the enzyme will cut.
"""

from __future__ import annotations

#: The bases each IUPAC code stands for.
CODES: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

#: Complement of every code, including the ambiguity codes.
_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)

#: Codes that stand for exactly one base.
UNAMBIGUOUS: frozenset[str] = frozenset("ACGT")

#: Deleting every code that stands for a single base leaves only the ambiguous
#: ones, which makes :func:`has_ambiguity` a C-speed scan rather than a Python
#: loop. The difference matters: the check decides whether a second pass over
#: the sequence is needed at all, and it is asked once per enzyme, so a
#: per-character loop would walk the whole molecule hundreds of times over on a
#: catalogue-wide digest. U is deleted alongside T because it also names one
#: base; an RNA sequence is not an uncertain one.
_AMBIGUITY_SIEVE = str.maketrans("", "", "ACGTUacgtu")


def has_ambiguity(sequence: str) -> bool:
    """Whether any residue stands for more than one base.

    Used to decide whether uncertain positions need to be looked for at all, so
    it must never report False for a sequence that does contain an ambiguity
    code: that would silently switch off the uncertainty reporting rather than
    change its speed.
    """
    return bool(sequence.translate(_AMBIGUITY_SIEVE))


def complement(sequence: str) -> str:
    """Complement a sequence, preserving case and ambiguity codes."""
    return sequence.translate(_COMPLEMENT)


def reverse_complement(sequence: str) -> str:
    """Reverse complement a sequence, preserving case and ambiguity codes."""
    return sequence.translate(_COMPLEMENT)[::-1]


def is_palindrome(pattern: str) -> bool:
    """Whether a recognition pattern reads the same on both strands."""
    return pattern.upper() == reverse_complement(pattern.upper())


def expand(code: str) -> frozenset[str]:
    """The bases an IUPAC code stands for.

    Raises:
        KeyError: if the character is not an IUPAC nucleotide code.
    """
    try:
        return CODES[code.upper()]
    except KeyError:
        raise KeyError(f"not an IUPAC nucleotide code: {code!r}") from None


def matches(site_code: str, sequence_code: str) -> bool:
    """Whether a sequence base certainly satisfies a site base.

    The test is containment, not overlap: the bases the sequence might be must
    all be bases the enzyme accepts. So ``N`` in a recognition site is
    satisfied by any base, while ``N`` in the sequence satisfies only ``N``.

    The asymmetry is the point. An enzyme with an ambiguous site really does
    cut every sequence its code covers, so widening there is correct. An
    ambiguous base in the sequence is missing information, and reporting a cut
    the enzyme may well not make is worse than reporting none: it sends someone
    to the bench with a digest that cannot work. Those positions are counted and
    reported separately as uncertain rather than silently dropped.
    """
    return expand(sequence_code) <= expand(site_code)


def could_match(site_code: str, sequence_code: str) -> bool:
    """Whether a sequence base might satisfy a site base.

    True when the two codes have any base in common, so an ``N`` in the
    sequence might match anything. Used to count the sites that cannot be
    called either way.
    """
    return bool(expand(site_code) & expand(sequence_code))
