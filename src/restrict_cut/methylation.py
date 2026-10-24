"""Whether a site is methylated in the DNA you are about to cut.

Plasmid grown in an ordinary laboratory *E. coli* strain is methylated before
you ever see it. Dam methyltransferase marks the adenine in GATC, Dcm marks the
inner cytosine in CCWGG, and DNA from a mammalian source carries 5-methyl\u00adcytosine
at CpG. An enzyme whose recognition site covers one of those bases may cut
poorly or not at all, which is a common and avoidable reason for a digest to
fail with no other explanation.

What this module reports, and what it does not
----------------------------------------------
It reports *overlap*: whether a methylated base falls inside the recognition
site, worked out from the sequence in hand and its flanking context. That is a
fact about the DNA, and it is the part a sequence can settle.

It deliberately does not claim whether a given enzyme is *blocked*. Overlap and
blocking are different questions, and the difference is not something the
sequence can answer: BamHI's site GGATCC always contains a Dam GATC and BamHI
cuts it regardless, MboI's site is GATC and Dam methylation blocks it outright,
and DpnI requires the methylation before it will cut at all. Sensitivity is an
empirical property of each enzyme, so this module points at the overlap and
leaves the verdict to the supplier's table.
"""

from __future__ import annotations

from dataclasses import dataclass

from functools import lru_cache

from restrict_cut.enzymes import Enzyme
from restrict_cut.iupac import could_match, expand
from restrict_cut.search import Site, compile_site


@dataclass(frozen=True)
class System:
    """A methyltransferase activity.

    Attributes:
        name: short name used on the command line.
        target: recognition sequence of the methyltransferase, in IUPAC codes.
        position: 0-based index within ``target`` of the base it modifies.
        base: the modified base as it is usually written.
        description: what carries this activity.
    """

    name: str
    target: str
    position: int
    base: str
    description: str

    @property
    def label(self) -> str:
        """Target with the modified base marked, as ``G<A>TC``."""
        return (
            f"{self.target[: self.position]}"
            f"<{self.target[self.position]}>"
            f"{self.target[self.position + 1:]}"
        )

    @property
    def readings(self) -> tuple[tuple[str, int], ...]:
        """The target as read on each strand, with the base it modifies.

        Both strands are always listed, even when the target is a palindrome,
        because the *modified base* need not sit symmetrically within it. Dam
        marks the adenine of GATC on both strands: in forward coordinates that
        is the second base for the top strand and the third for the bottom, so
        a site covering only the third base still carries a methylated base. A
        reading is dropped only when it would repeat another exactly.
        """
        from restrict_cut.iupac import reverse_complement

        width = len(self.target)
        seen: list[tuple[str, int]] = [(self.target, self.position)]
        mirrored = (reverse_complement(self.target), width - 1 - self.position)
        if mirrored not in seen:
            seen.append(mirrored)
        return tuple(seen)


#: The activities a sequence can be checked against.
SYSTEMS: dict[str, System] = {
    "dam": System(
        "dam",
        "GATC",
        1,
        "N6-methyladenine",
        "Dam methyltransferase, active in most laboratory E. coli strains",
    ),
    "dcm": System(
        "dcm",
        "CCWGG",
        1,
        "5-methylcytosine",
        "Dcm methyltransferase, active in most laboratory E. coli strains",
    ),
    "cpg": System(
        "cpg",
        "CG",
        0,
        "5-methylcytosine",
        "CpG methylation, carried by DNA from a mammalian source",
    ),
}

#: Activities checked unless the caller asks for others. Dam and Dcm are the
#: two that act on plasmid grown in a standard cloning strain; CpG is left off
#: because it does not apply to DNA prepared from E. coli.
DEFAULT_SYSTEMS: tuple[str, ...] = ("dam", "dcm")


def resolve_systems(names: list[str] | None) -> list[System]:
    """Look up methylation systems by name.

    Raises:
        KeyError: if a name is not one of :data:`SYSTEMS`.
    """
    if names is None:
        names = list(DEFAULT_SYSTEMS)
    chosen = []
    for name in names:
        key = name.strip().lower()
        if key not in SYSTEMS:
            known = ", ".join(sorted(SYSTEMS))
            raise KeyError(f"unknown methylation system: {name!r}; known: {known}")
        chosen.append(SYSTEMS[key])
    return chosen


#: How a recognition site relates to a methylation target, before any sequence
#: is considered.
EXPOSURE = ("always", "possible", "never")


def exposure(enzyme: Enzyme, system: System) -> str:
    """Whether the site can carry a methylated base, ignoring any sequence.

    Returns:
        ``always`` when every sequence the enzyme recognises contains a
        methylated base inside the site, whatever surrounds it; ``never`` when
        no arrangement of flanking DNA can put one there; ``possible``
        otherwise, meaning it depends on the context of each individual site.
    """
    pattern = enzyme.pattern
    span = len(pattern)

    # 'always': on some strand an alignment lies wholly inside the site, and
    # every base the pattern permits there is one the target accepts, so the
    # modified base is present whatever surrounds the site.
    for target, _ in system.readings:
        for start in range(span - len(target) + 1):
            if all(
                expand(pattern[start + index]) <= expand(target[index])
                for index in range(len(target))
            ):
                return "always"

    # 'possible': some alignment, perhaps hanging off the end into flanking DNA,
    # is compatible with the pattern and puts the modified base inside the site.
    for target, position in system.readings:
        for start in range(-len(target) + 1, span):
            modified = start + position
            if not 0 <= modified < span:
                continue
            overlaps = range(max(0, start), min(span, start + len(target)))
            if all(
                could_match(target[index - start], pattern[index])
                for index in overlaps
            ):
                return "possible"

    return "never"


@dataclass(frozen=True)
class Overlap:
    """A methylated base falling inside one recognition site.

    Attributes:
        enzyme: the enzyme whose site is affected.
        site_start: 1-based start of that site on the forward strand.
        system: name of the methyltransferase activity.
        position: 1-based position of the modified base.
        target_start: 1-based start of the methyltransferase's own target.
    """

    enzyme: str
    site_start: int
    system: str
    position: int
    target_start: int


@lru_cache(maxsize=16)
def _targets(
    sequence: str, system: System, circular: bool
) -> tuple[tuple[int, int], ...]:
    """Find a methyltransferase's targets, remembering the last few answers.

    Where the targets fall depends only on the sequence and the activity, never
    on the enzyme being checked against them, but the natural way to write the
    caller asks once per enzyme. Without this the same scan runs several
    hundred times over for a catalogue-wide digest, and it dominates the run.
    The cache is small and keyed on the sequence itself, so it holds a reference
    to at most a handful of molecules the caller already has in hand.
    """
    length = len(sequence)
    width = len(system.target)
    if not length or width > length:
        return ()

    reach = min(width - 1, length) if circular else 0
    window = sequence + sequence[:reach]
    found: set[tuple[int, int]] = set()

    for pattern, offset in system.readings:
        for match in compile_site(pattern).finditer(window):
            index = match.start()
            if index >= length:
                break
            if not circular and index + width > length:
                continue
            found.add((index % length + 1, (index + offset) % length + 1))

    return tuple(sorted(found))


def find_targets(
    sequence: str, system: System, circular: bool = False
) -> list[tuple[int, int]]:
    """Every methyltransferase target, as ``(target start, modified base)``.

    Both are 1-based forward-strand positions. The target is read on both
    strands, so a modified base that sits asymmetrically within an otherwise
    palindromic target - as Dam's adenine does within GATC - is reported at its
    position on each strand.
    """
    return list(_targets(sequence, system, circular))


def site_overlaps(
    sequence: str,
    enzyme: Enzyme,
    sites: list[Site],
    systems: list[System],
    circular: bool = False,
) -> list[Overlap]:
    """Methylated bases falling inside the given recognition sites."""
    length = len(sequence)
    overlaps: list[Overlap] = []
    for system in systems:
        if exposure(enzyme, system) == "never":
            continue
        targets = find_targets(sequence, system, circular)
        if not targets:
            continue
        modified_positions = {modified: start for start, modified in targets}
        for site in sites:
            covered = _site_positions(site, enzyme.length, length, circular)
            for position in covered:
                if position in modified_positions:
                    overlaps.append(
                        Overlap(
                            enzyme=enzyme.name,
                            site_start=site.start,
                            system=system.name,
                            position=position,
                            target_start=modified_positions[position],
                        )
                    )
    overlaps.sort(key=lambda item: (item.site_start, item.system, item.position))
    return overlaps


def _site_positions(
    site: Site, width: int, length: int, circular: bool
) -> list[int]:
    """The 1-based positions a recognition site occupies."""
    if circular:
        return [(site.start - 1 + step) % length + 1 for step in range(width)]
    return [site.start + step for step in range(width) if site.start + step <= length]


def context(sequence: str, position: int, flank: int, circular: bool) -> str:
    """A window of sequence around one position, for showing a methylated base.

    On a circular molecule the window wraps past the origin, and is clamped to
    the whole molecule once it would ask for more bases than there are: without
    that, a window wider than the sequence wraps past its own start and comes
    back *shorter* than a narrower one, which is how a five-base circle asked
    for seven bases of context returned two.
    """
    length = len(sequence)
    if not length:
        return ""
    if circular:
        if 2 * flank + 1 >= length:
            # Rotate the whole molecule so the position sits in the middle,
            # which is as much context as a circle this small can offer.
            start = (position - 1 - length // 2) % length
            return sequence[start:] + sequence[:start]
        start = (position - 1 - flank) % length
        return "".join(
            sequence[(start + step) % length] for step in range(2 * flank + 1)
        )
    first = max(1, position - flank)
    last = min(length, position + flank)
    return sequence[first - 1: last]
