"""The enzyme catalogue: specificities, filters and relationships.

The specificities live in the generated :mod:`restrict_cut.enzyme_data` module,
built from REBASE by ``scripts/build_enzyme_table.py``. This module wraps each
row in an :class:`Enzyme` that can describe itself, and derives the
relationships REBASE leaves implicit: which enzymes share a specificity, and
which share one but cut it in a different place.
"""

from __future__ import annotations

import difflib
from dataclasses import dataclass
from functools import lru_cache

from restrict_cut.enzyme_data import ENZYMES, REBASE_RELEASE, SUPPLIERS
from restrict_cut.iupac import CODES, expand, is_palindrome, reverse_complement

__all__ = [
    "END_TYPES",
    "Enzyme",
    "REBASE_RELEASE",
    "SUPPLIERS",
    "all_enzymes",
    "compatible_groups",
    "get",
    "isoschizomers",
    "neoschizomers",
    "resolve",
    "select",
    "validate_catalogue",
]

#: End left by a cut, named for the strand that overhangs.
END_TYPES = ("5' overhang", "blunt", "3' overhang")


@dataclass(frozen=True)
class Enzyme:
    """One restriction enzyme specificity.

    Attributes:
        name: the enzyme name as REBASE writes it, such as ``BamHI``.
        pattern: recognition sequence in IUPAC codes.
        cuts: one ``(top, bottom)`` offset pair per double-strand break, each
            counting the recognition-site bases lying 5' of the cut. ``0`` is
            immediately before the site, so a negative offset reaches upstream.
            Both are in the top-strand frame, so ``bottom - top`` is the
            overhang.
        suppliers: REBASE one-letter codes of the companies selling the enzyme.
    """

    name: str
    pattern: str
    cuts: tuple[tuple[int, int], ...]
    suppliers: str = ""

    @property
    def length(self) -> int:
        """Length of the recognition sequence."""
        return len(self.pattern)

    @property
    def is_palindrome(self) -> bool:
        """Whether the recognition sequence reads the same on both strands."""
        return is_palindrome(self.pattern)

    @property
    def is_ambiguous(self) -> bool:
        """Whether the recognition sequence contains an ambiguity code."""
        return any(base not in "ACGT" for base in self.pattern)

    @property
    def specific_bases(self) -> int:
        """How many site bases are fixed to a single possibility.

        This, not the site length, governs how often a site is expected to
        occur whenever the site contains ambiguity codes or spacer Ns.
        """
        return sum(1 for base in self.pattern if len(expand(base)) == 1)

    @property
    def matching_sequences(self) -> int:
        """How many distinct sequences of this length the enzyme recognises.

        Counted over both strands, which is what governs how often a site
        occurs. Doubling the one-strand count is wrong whenever the pattern and
        its reverse complement accept some of the same sequences: Hin4I's
        GAYNNNNNVTC overlaps its own reverse complement, so doubling counts
        those sequences twice and makes the enzyme look twice as common as it
        is. The overlap is subtracted instead, position by position.
        """
        one_strand = 1
        for base in self.pattern:
            one_strand *= len(expand(base))
        if self.is_palindrome:
            return one_strand
        flipped = reverse_complement(self.pattern)
        both = 1
        for base, mirrored in zip(self.pattern, flipped):
            both *= len(expand(base) & expand(mirrored))
        return 2 * one_strand - both

    @property
    def expected_frequency(self) -> float:
        """Expected bases between sites in random, equal-composition DNA.

        The reciprocal of the chance that any given position starts a site,
        which is :attr:`matching_sequences` out of the ``4 ** n`` sequences of
        that length.
        """
        return 4.0 ** self.length / self.matching_sequences

    @property
    def cuts_outside_site(self) -> bool:
        """Whether any cut falls outside the recognition sequence.

        True for the Type IIS enzymes used in Golden Gate assembly, whose
        overhang is spelled by the DNA rather than by the enzyme.
        """
        return any(
            top < 0 or top > self.length or bottom < 0 or bottom > self.length
            for top, bottom in self.cuts
        )

    @property
    def excises_site(self) -> bool:
        """Whether the enzyme cuts on both sides, excising its own site."""
        return len(self.cuts) > 1

    def overhang_length(self, cut: int = 0) -> int:
        """Signed overhang of one cut: positive 5', negative 3', zero blunt."""
        top, bottom = self.cuts[cut]
        return bottom - top

    def end_type(self, cut: int = 0) -> str:
        """Name the end one cut leaves."""
        overhang = self.overhang_length(cut)
        if overhang > 0:
            return "5' overhang"
        if overhang < 0:
            return "3' overhang"
        return "blunt"

    @property
    def fixed_overhang(self) -> str | None:
        """The overhang left at a site read in the pattern's own orientation.

        ``None`` when the enzyme cuts outside its recognition site, or when the
        overhang falls on ambiguity codes, because then the bases left single
        stranded depend on the DNA rather than on the enzyme. The empty string
        means a blunt end.

        For a site that is not a palindrome this is only half the story: the
        same enzyme meeting its site on the other strand leaves the reverse
        complement of this sequence. BauI reads CACGAG and leaves 5'-ACGA at a
        forward site and 5'-TCGT at a reverse one, and those two ends do not
        anneal to each other. See :attr:`overhang_depends_on_orientation`.
        """
        if self.excises_site or self.cuts_outside_site:
            return None
        top, bottom = self.cuts[0]
        overhang = self.pattern[min(top, bottom): max(top, bottom)]
        if any(base not in "ACGT" for base in overhang):
            return None
        return overhang

    @property
    def overhang_depends_on_orientation(self) -> bool:
        """Whether the overhang differs between a forward and a reverse site.

        True for a sticky-ended enzyme whose recognition sequence is not a
        palindrome. Such an enzyme leaves one of two complementary overhangs
        depending on which strand carries the site, so which ends it produces
        cannot be settled without a sequence.
        """
        overhang = self.fixed_overhang
        return bool(overhang) and not self.is_palindrome

    @property
    def searches_one_strand(self) -> bool:
        """Whether searching the forward strand alone finds every cut.

        Reverse complementing a site of length ``L`` maps an offset pair
        ``(top, bottom)`` to ``(L - bottom, L - top)``. When the whole set of
        pairs survives that map unchanged, binding in either orientation breaks
        the same bonds and a second pass would only find duplicates.

        This is not the same question as whether the site is a palindrome.
        HauII reads TGGCCA on both strands but cuts only downstream of it, so it
        makes a different break in each orientation and both passes are needed.
        """
        if not self.is_palindrome:
            return False
        flipped = {(self.length - bottom, self.length - top) for top, bottom in self.cuts}
        return flipped == set(self.cuts)

    @property
    def notation(self) -> str:
        """The specificity in REBASE's notation, such as ``G^AATTC``.

        A mirrored cut inside a palindromic site is written with a caret.
        Anything else is written as offsets from the ends of the site, because
        a caret cannot say where the bottom strand breaks.
        """
        top, bottom = self.cuts[0]
        if self.searches_one_strand and len(self.cuts) == 1 and 0 <= top <= self.length:
            return f"{self.pattern[:top]}^{self.pattern[top:]}"
        if len(self.cuts) == 1:
            return f"{self.pattern}({top - self.length}/{bottom - self.length})"
        second_top, second_bottom = self.cuts[1]
        return (
            f"({-top}/{-bottom}){self.pattern}"
            f"({second_top - self.length}/{second_bottom - self.length})"
        )

    @property
    def supplier_names(self) -> tuple[str, ...]:
        """Names of the companies selling the enzyme."""
        return tuple(SUPPLIERS[code] for code in self.suppliers if code in SUPPLIERS)

    @property
    def is_available(self) -> bool:
        """Whether any company sells the enzyme."""
        return bool(self.suppliers)


@lru_cache(maxsize=1)
def all_enzymes() -> dict[str, Enzyme]:
    """Every enzyme in the catalogue, by name."""
    return {
        name: Enzyme(name, pattern, cuts, suppliers)
        for name, (pattern, cuts, suppliers) in ENZYMES.items()
    }


@lru_cache(maxsize=1)
def _folded() -> dict[str, str]:
    """Case-insensitive index of enzyme names, for forgiving lookup."""
    return {name.casefold(): name for name in ENZYMES}


def get(name: str) -> Enzyme:
    """Look up one enzyme by name, ignoring case.

    Raises:
        KeyError: naming the closest alternatives, if the name is unknown.
    """
    catalogue = all_enzymes()
    if name in catalogue:
        return catalogue[name]
    folded = _folded().get(name.strip().casefold())
    if folded is not None:
        return catalogue[folded]
    close = difflib.get_close_matches(name, catalogue, n=3, cutoff=0.6)
    hint = f"; did you mean {', '.join(close)}?" if close else ""
    raise KeyError(f"unknown enzyme: {name!r}{hint}")


def resolve(names: list[str]) -> list[Enzyme]:
    """Look up several enzymes, keeping the order given and dropping repeats."""
    seen: dict[str, Enzyme] = {}
    for name in names:
        enzyme = get(name)
        seen.setdefault(enzyme.name, enzyme)
    return list(seen.values())


def _specificity(enzyme: Enzyme) -> tuple[str, tuple[tuple[int, int], ...]]:
    return enzyme.pattern, enzyme.cuts


@lru_cache(maxsize=1)
def _by_specificity() -> dict[tuple[str, tuple[tuple[int, int], ...]], tuple[str, ...]]:
    groups: dict[tuple[str, tuple[tuple[int, int], ...]], list[str]] = {}
    for enzyme in all_enzymes().values():
        groups.setdefault(_specificity(enzyme), []).append(enzyme.name)
    return {key: tuple(sorted(names)) for key, names in groups.items()}


@lru_cache(maxsize=1)
def _by_pattern() -> dict[str, tuple[str, ...]]:
    groups: dict[str, list[str]] = {}
    for enzyme in all_enzymes().values():
        groups.setdefault(enzyme.pattern, []).append(enzyme.name)
    return {key: tuple(sorted(names)) for key, names in groups.items()}


def isoschizomers(enzyme: Enzyme) -> tuple[str, ...]:
    """Other enzymes recognising the same site and cutting it the same way.

    These are interchangeable at the bench, so this is the list to reach for
    when the enzyme a protocol names is not the one in the freezer.
    """
    shared = _by_specificity().get(_specificity(enzyme), ())
    return tuple(name for name in shared if name != enzyme.name)


def neoschizomers(enzyme: Enzyme) -> tuple[str, ...]:
    """Enzymes recognising the same site but cutting it in a different place.

    SmaI and XmaI both recognise CCCGGG. SmaI cuts it in the middle for a blunt
    end, XmaI after the first C for a 5' overhang, so they are not
    interchangeable even though they cut at the same place in a map.
    """
    same_site = _by_pattern().get(enzyme.pattern, ())
    cuts = enzyme.cuts
    catalogue = all_enzymes()
    return tuple(name for name in same_site if catalogue[name].cuts != cuts)


def compatible_groups(enzymes: list[Enzyme]) -> dict[str, tuple[str, ...]]:
    """Group enzymes by the end they leave, where the enzyme alone fixes it.

    Two things keep an enzyme out of these groups, and both are the same
    objection: the end it leaves is not a property of the enzyme.

    An enzyme that cuts outside its recognition site, such as BsaI, takes its
    overhang from whatever DNA happens to sit there. An enzyme with a sticky end
    and a non-palindromic site, such as BauI, leaves 5'-ACGA where its site
    reads forward and 5'-TCGT where it reads backward, and those two do not
    anneal; grouping it with BsiI on the strength of the forward reading would
    promise a pairing that half the sites will not deliver.

    Either way the question can only be settled against a sequence, which is
    what :func:`restrict_cut.selection.compatible_with` is for.
    """
    groups: dict[str, list[str]] = {}
    for enzyme in enzymes:
        overhang = enzyme.fixed_overhang
        if overhang is None or enzyme.overhang_depends_on_orientation:
            continue
        key = f"{enzyme.end_type()} {overhang}" if overhang else "blunt"
        groups.setdefault(key, []).append(enzyme.name)
    return {key: tuple(sorted(names)) for key, names in sorted(groups.items())}


def select(
    available_only: bool = False,
    min_length: int = 0,
    max_length: int = 0,
    end_type: str | None = None,
    supplier: str | None = None,
) -> list[Enzyme]:
    """Filter the catalogue.

    Args:
        available_only: keep only enzymes a company sells.
        min_length: keep sites at least this many bases long, 0 for no bound.
        max_length: keep sites at most this many bases long, 0 for no bound.
        end_type: keep enzymes leaving this end, as named by :data:`END_TYPES`.
        supplier: keep enzymes sold by this REBASE supplier code.
    """
    if end_type is not None and end_type not in END_TYPES:
        raise ValueError(
            f"unknown end type: {end_type!r}; choose from {', '.join(END_TYPES)}"
        )
    if supplier is not None and supplier not in SUPPLIERS:
        known = ", ".join(sorted(SUPPLIERS))
        raise ValueError(f"unknown supplier code: {supplier!r}; known codes are {known}")
    chosen = []
    for enzyme in all_enzymes().values():
        if available_only and not enzyme.is_available:
            continue
        if min_length and enzyme.length < min_length:
            continue
        if max_length and enzyme.length > max_length:
            continue
        if end_type is not None and enzyme.end_type() != end_type:
            continue
        if supplier is not None and supplier not in enzyme.suppliers:
            continue
        chosen.append(enzyme)
    return sorted(chosen, key=lambda enzyme: enzyme.name)


def validate_catalogue() -> list[str]:
    """Check the generated table for internal contradictions.

    Returns a list of problems, empty when the table is sound. Kept beside the
    data rather than only in the tests, so a regenerated table can be checked
    without running the suite.
    """
    problems = []
    for name, enzyme in all_enzymes().items():
        if not enzyme.pattern:
            problems.append(f"{name}: empty recognition sequence")
        for base in enzyme.pattern:
            if base not in CODES:
                problems.append(f"{name}: {base!r} is not an IUPAC code")
        if not enzyme.cuts:
            problems.append(f"{name}: no cut positions")
        for code in enzyme.suppliers:
            if code not in SUPPLIERS:
                problems.append(f"{name}: unknown supplier code {code!r}")
    return problems
