"""Choosing an enzyme for a job.

Mapping a digest answers 'where does this enzyme cut'. Cloning asks the
question the other way round: *which* enzyme cuts here once and nowhere else,
or cuts this vector but leaves my insert alone. That search runs the whole
catalogue against one or more sequences, so it is kept apart from
:mod:`restrict_cut.digest`, which applies a named set of enzymes to one.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from restrict_cut.digest import Cut, cuts_from_site
from restrict_cut.enzymes import Enzyme
from restrict_cut.iupac import has_ambiguity
from restrict_cut.search import count_uncertain, find_sites


@dataclass
class Survey:
    """How often each enzyme cuts one sequence.

    Attributes:
        name: identifier of the sequence surveyed.
        length: length of the molecule.
        circular: whether it was treated as circular.
        counts: cuts made, by enzyme name.
        positions: top-strand nick positions, by enzyme name.
        overhangs: the distinct ends observed, by enzyme name, each as
            ``(end type, overhang)``. Read from the sequence rather than from
            the recognition site, which is the only way to be right about an
            enzyme whose site is not a palindrome: the same enzyme leaves
            complementary overhangs at a forward and a reverse site.
        uncertain: positions that could not be called, by enzyme name.
    """

    name: str
    length: int
    circular: bool
    counts: dict[str, int] = field(default_factory=dict)
    positions: dict[str, tuple[int, ...]] = field(default_factory=dict)
    overhangs: dict[str, tuple[tuple[str, str], ...]] = field(default_factory=dict)
    uncertain: dict[str, int] = field(default_factory=dict)

    def cutting(self, times: int) -> list[str]:
        """Enzymes cutting exactly this many times, in catalogue order."""
        return sorted(name for name, count in self.counts.items() if count == times)


def survey(
    sequence: str,
    enzymes: list[Enzyme],
    circular: bool = False,
    name: str = "sequence",
) -> Survey:
    """Count how often every enzyme given cuts one sequence."""
    length = len(sequence)
    result = Survey(name=name, length=length, circular=circular)
    ambiguous = has_ambiguity(sequence)

    for enzyme in enzymes:
        sites = find_sites(sequence, enzyme, circular)
        nicks: set[tuple[int, int]] = set()
        ends: set[tuple[str, str]] = set()
        for site in sites:
            for cut in cuts_from_site(length, enzyme, site, circular):
                if (cut.top, cut.bottom) in nicks:
                    continue
                nicks.add((cut.top, cut.bottom))
                ends.add(end_signature(sequence, cut, circular))
        result.counts[enzyme.name] = len(nicks)
        result.positions[enzyme.name] = tuple(sorted(top for top, _ in nicks))
        result.overhangs[enzyme.name] = tuple(sorted(ends))
        if ambiguous:
            uncertain = count_uncertain(sequence, enzyme, circular)
            if uncertain:
                result.uncertain[enzyme.name] = uncertain

    return result


@dataclass(frozen=True)
class Candidate:
    """One enzyme that meets a selection.

    Attributes:
        enzyme: the enzyme itself.
        cuts: how often it cuts the target.
        positions: top-strand nick positions in the target.
        overhangs: the distinct ends this enzyme leaves in the target, read
            from the sequence rather than assumed from the recognition site.
        uncertain: target positions that could not be called.
        protected_cuts: cuts found in each sequence that had to be left alone.
            Always zero for a candidate that is returned, since any cut
            disqualifies the enzyme; kept so a caller can see what was checked.
    """

    enzyme: Enzyme
    cuts: int
    positions: tuple[int, ...]
    overhangs: tuple[tuple[str, str], ...] = ()
    uncertain: int = 0
    protected_cuts: tuple[tuple[str, int], ...] = ()

    @property
    def is_clean(self) -> bool:
        """Whether no ambiguous position could have changed the count."""
        return not self.uncertain


@dataclass
class Selection:
    """The outcome of searching a catalogue against one target.

    Attributes:
        target: identifier of the sequence being cut.
        considered: how many enzymes were searched.
        candidates: the enzymes that met the selection.
        protected: names of the sequences that had to be left untouched.
        excluded_uncertain: enzymes rejected *only* because a protected
            sequence carries an ambiguity code where they might cut. They are
            not safe to recommend, but neither are they known to cut, so they
            are named rather than silently dropped.
    """

    target: str
    considered: int
    candidates: list[Candidate] = field(default_factory=list)
    protected: tuple[str, ...] = ()
    excluded_uncertain: tuple[str, ...] = ()


def choose(
    sequence: str,
    enzymes: list[Enzyme],
    cuts: int | None = 1,
    max_cuts: int | None = None,
    circular: bool = False,
    name: str = "sequence",
    protect: list[tuple[str, str, bool]] | None = None,
) -> Selection:
    """Find the enzymes that cut a target a given number of times.

    Args:
        sequence: the molecule to cut, usually the vector.
        enzymes: the catalogue to search.
        cuts: required number of cuts, or None to accept any number.
        max_cuts: accept any count up to this, instead of an exact ``cuts``.
        circular: treat the target as circular.
        name: identifier recorded on the result.
        protect: sequences that must *not* be cut, each as
            ``(name, sequence, circular)``. An insert belongs here: an enzyme
            that cuts it is no use for moving it, however well it cuts the
            vector.

    An enzyme is kept out of the result if a protected sequence carries a
    confirmed site *or* a position that cannot be called because the sequence
    is ambiguous there. Protection is a guarantee, so the benefit of the doubt
    has to run the other way from the target: an enzyme that might cut the
    insert is not one to hand someone, and a protected sequence written with an
    N where the site would be is exactly the case that looks safe and is not.

    Returns:
        A :class:`Selection`, whose candidates are sorted by cut count, then by how rare the site is in random
        DNA, then by how many companies sell the enzyme, then by name. Site
        length is deliberately not the sort key: AhdI's GACNNNNNGTC is eleven
        bases long but only six of them are fixed, so it is no rarer than
        EcoRI's GAATTC and is harder to buy.
    """
    if cuts is not None and cuts < 0:
        raise ValueError(f"cut count must not be negative: {cuts}")
    if max_cuts is not None and max_cuts < 0:
        raise ValueError(f"maximum cut count must not be negative: {max_cuts}")

    target = survey(sequence, enzymes, circular, name)

    protected: list[tuple[str, Survey]] = []
    for protect_name, protect_sequence, protect_circular in protect or []:
        protected.append(
            (
                protect_name,
                survey(protect_sequence, enzymes, protect_circular, protect_name),
            )
        )

    candidates: list[Candidate] = []
    uncertain_elsewhere: list[str] = []
    for enzyme in enzymes:
        count = target.counts[enzyme.name]
        if max_cuts is not None:
            if count > max_cuts:
                continue
        elif cuts is not None and count != cuts:
            continue

        elsewhere = tuple(
            (protect_name, survey_result.counts[enzyme.name])
            for protect_name, survey_result in protected
        )
        if any(number for _, number in elsewhere):
            continue
        if any(
            survey_result.uncertain.get(enzyme.name, 0)
            for _, survey_result in protected
        ):
            uncertain_elsewhere.append(enzyme.name)
            continue

        candidates.append(
            Candidate(
                enzyme=enzyme,
                cuts=count,
                positions=target.positions[enzyme.name],
                overhangs=target.overhangs[enzyme.name],
                uncertain=target.uncertain.get(enzyme.name, 0),
                protected_cuts=elsewhere,
            )
        )

    candidates.sort(
        key=lambda item: (
            item.cuts,
            -item.enzyme.expected_frequency,
            -len(item.enzyme.suppliers),
            item.enzyme.name,
        )
    )
    return Selection(
        target=name,
        considered=len(enzymes),
        candidates=candidates,
        protected=tuple(protect_name for protect_name, _ in protected),
        excluded_uncertain=tuple(sorted(uncertain_elsewhere)),
    )


def end_signature(sequence: str, cut: Cut, circular: bool = False) -> tuple[str, str]:
    """The end a cut leaves, as ``(end type, overhang)``.

    This pair is everything that decides whether two ends will anneal, and it
    is read from the molecule the cut is actually in. Reducing a cut to its
    signature before comparing is what lets ends from two different molecules -
    a vector and an insert, which is the case that matters - be compared at all.

    The overhang is upper-cased. Site matching already ignores case, so a
    sequence written with soft-masked regions in lower case finds its sites
    either way; comparing the overhangs as written would then report BamHI's
    ``gatc`` and BglII's ``GATC`` as incompatible, which is a fact about the
    file's capitalisation rather than about the DNA.
    """
    return cut.end_type, cut.overhang(sequence, circular).upper()


def ends_compatible(first: tuple[str, str], second: tuple[str, str]) -> bool:
    """Whether two end signatures will ligate.

    Ends anneal when they expose the same bases on the same strand, so both the
    overhang text and the end type have to agree. Comparing the forward-strand
    text between the nicks is what makes BamHI and BglII interchangeable at a
    junction despite recognising different sites, and what keeps a 5' GATC
    apart from a 3' GATC.
    """
    return first == second


def compatible_with(
    sequence: str,
    cut: Cut,
    others: list[tuple[Cut, str]],
    circular: bool = False,
) -> list[str]:
    """Which of the other cuts leave an end that will ligate to this one.

    Every cut in ``others`` must be a cut in ``sequence``; their overhangs are
    read from it. To compare ends across two molecules - a vector end against an
    insert end - take :func:`end_signature` of each cut against its own
    sequence and pass the pair to :func:`ends_compatible`. Reading a cut's
    overhang off a molecule it does not belong to silently returns whatever
    bases happen to sit at those coordinates, which for a non-palindromic site
    read on the other strand is exactly the complementary overhang, and so
    reports a pairing that will not happen.
    """
    reference = end_signature(sequence, cut, circular)
    return [
        label
        for other, label in others
        if ends_compatible(end_signature(sequence, other, circular), reference)
    ]
