"""Finding recognition sites in a sequence.

Coordinates are 1-based and inclusive on the forward strand, and a site is
always reported by its forward-strand span whichever strand it was matched on.

Three things here are easy to get wrong and are handled explicitly:

* A site that is not a palindrome occurs on either strand, so both are
  searched. One that is a palindrome *and* is cut the same way in either
  orientation is searched once, because the second pass would only rediscover
  the same breaks.
* On a circular molecule a site may straddle the origin, so the search runs
  over the sequence with its own start appended.
* An ambiguity code in the *sequence* means the base is not known. Such a
  position cannot be called a site, but neither is it certainly not one, so it
  is counted separately rather than silently dropped.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Callable

from restrict_cut.enzymes import Enzyme
from restrict_cut.iupac import (
    CODES,
    could_match,
    has_ambiguity,
    matches,
    reverse_complement,
)

#: Strand labels. ``=`` marks a site read identically on both strands, so that
#: a palindromic site is not reported twice.
STRANDS = ("+", "-", "=")


@dataclass(frozen=True)
class Site:
    """One occurrence of a recognition sequence.

    Attributes:
        enzyme: name of the enzyme that recognises it.
        start: 1-based first base of the site on the forward strand.
        end: 1-based last base of the site on the forward strand, inclusive.
        strand: ``+``, ``-``, or ``=`` for a site that reads the same either
            way. On a circular molecule ``end`` may be smaller than ``start``,
            which means the site spans the origin.
    """

    enzyme: str
    start: int
    end: int
    strand: str

    @property
    def spans_origin(self) -> bool:
        """Whether the site wraps past the end of a circular sequence."""
        return self.end < self.start


def _character_class(site_code: str, admit: Callable[[str, str], bool]) -> str:
    """Every sequence code a site code accepts, as a regex character class.

    The class is built by asking the matching predicate about each of the
    sixteen IUPAC codes rather than by listing bases. Listing bases would be
    wrong: the sequence may itself be written in ambiguity codes, and an ``N``
    in the sequence is a character the class has to admit explicitly when the
    site is ``N`` - it is not one of ``A``, ``C``, ``G``, ``T``. Deriving the
    class from the predicate keeps this fast path and
    :func:`restrict_cut.iupac.matches` from drifting apart.
    """
    allowed = sorted(code for code in CODES if admit(site_code, code))
    both_cases = "".join(allowed) + "".join(allowed).lower()
    return f"[{both_cases}]"


def _compile(pattern: str, admit: Callable[[str, str], bool]) -> re.Pattern[str]:
    """Compile a recognition pattern under one matching rule.

    The lookahead keeps overlapping occurrences, which matter for sites such as
    SgeI's CNNG that can tile along a sequence.
    """
    classes = "".join(_character_class(code, admit) for code in pattern)
    return re.compile(f"(?=({classes}))")


def _pattern_regex(pattern: str) -> re.Pattern[str]:
    """Compile a pattern that matches only certain sites."""
    return _compile(pattern, matches)


def compile_site(pattern: str) -> re.Pattern[str]:
    """Compile a pattern that matches only where the sequence certainly satisfies it.

    Shared with :mod:`restrict_cut.methylation`, which has the same job to do
    for a methyltransferase target as this module has for a recognition site.
    Both must agree with :func:`restrict_cut.iupac.matches`, so both go through
    one compiler rather than each walking the sequence their own way.
    """
    return _pattern_regex(pattern)


def _possible_regex(pattern: str) -> re.Pattern[str]:
    """Compile a pattern that also matches where the sequence is ambiguous.

    Used only to count positions that cannot be called either way: every
    sequence code sharing any base with the site code is admitted.
    """
    return _compile(pattern, could_match)


def _search_window(sequence: str, length: int, circular: bool) -> str:
    """The text to search, extended so a circular site can span the origin.

    A site of length ``L`` can begin at most ``L - 1`` bases before the origin
    and still cross it, so that much of the start is appended.
    """
    if not circular or length < 2 or not sequence:
        return sequence
    overlap = min(length - 1, len(sequence))
    return sequence + sequence[:overlap]


def find_sites(
    sequence: str,
    enzyme: Enzyme,
    circular: bool = False,
) -> list[Site]:
    """Every certain occurrence of one enzyme's recognition site.

    Sites are returned sorted by forward-strand start. On a circular molecule a
    site spanning the origin is reported with ``end < start``.
    """
    length = len(sequence)
    if not length or enzyme.length > length:
        return []

    window = _search_window(sequence, enzyme.length, circular)
    found: list[Site] = []
    strands: tuple[tuple[str, str], ...]
    if enzyme.searches_one_strand:
        strands = (("=", enzyme.pattern),)
    else:
        strands = (("+", enzyme.pattern), ("-", reverse_complement(enzyme.pattern)))

    for strand, pattern in strands:
        for match in _pattern_regex(pattern).finditer(window):
            offset = match.start()
            if offset >= length:
                # Already reported from the unwrapped part of the window.
                break
            start = offset + 1
            end = offset + enzyme.length
            if end > length:
                if not circular:
                    continue
                end -= length
            found.append(Site(enzyme.name, start, end, strand))

    found.sort(key=lambda site: (site.start, site.strand))
    return found


def count_uncertain(sequence: str, enzyme: Enzyme, circular: bool = False) -> int:
    """How many positions might be sites but cannot be called.

    A position counts when every base is compatible with the recognition
    sequence but at least one of them is an ambiguity code in the *sequence*,
    so whether the enzyme cuts there depends on information the file does not
    carry. Positions that are certain sites are excluded.
    """
    if not has_ambiguity(sequence):
        return 0
    length = len(sequence)
    if not length or enzyme.length > length:
        return 0

    window = _search_window(sequence, enzyme.length, circular)
    certain = {(site.strand, site.start) for site in find_sites(sequence, enzyme, circular)}
    strands: tuple[tuple[str, str], ...]
    if enzyme.searches_one_strand:
        strands = (("=", enzyme.pattern),)
    else:
        strands = (("+", enzyme.pattern), ("-", reverse_complement(enzyme.pattern)))

    uncertain = 0
    for strand, pattern in strands:
        for match in _possible_regex(pattern).finditer(window):
            offset = match.start()
            if offset >= length:
                break
            if not circular and offset + enzyme.length > length:
                continue
            if (strand, offset + 1) in certain:
                continue
            uncertain += 1
    return uncertain
