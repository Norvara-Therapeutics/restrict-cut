"""Turning recognition sites into cuts, and cuts into fragments.

A cut is recorded as the pair of nick positions it leaves, each named by the
base immediately 5' of it on the forward strand: a nick written ``n`` falls
between base ``n`` and base ``n + 1``. Keeping both nicks rather than one
position is what lets the overhang be reported, and lets a 3' overhang be told
from a 5' one.

Fragment lengths are measured between the *top-strand* nicks. For sticky ends
the two strands of a fragment differ in length, so some convention is needed;
this one has the property that the fragment lengths of a complete digest sum to
the length of the molecule, which is what makes the output checkable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from restrict_cut.enzymes import Enzyme
from restrict_cut.iupac import has_ambiguity
from restrict_cut.search import Site, count_uncertain, find_sites

#: How the end of a fragment arose.
TERMINUS = "end"


def slice_sequence(sequence: str, start: int, end: int, circular: bool = False) -> str:
    """Take a 1-based inclusive span, wrapping past the origin when circular."""
    length = len(sequence)
    if not length:
        return ""
    if start <= end:
        return sequence[start - 1: end]
    if not circular:
        raise ValueError(f"span {start}..{end} runs backwards on a linear sequence")
    return sequence[start - 1:] + sequence[:end]


def span_length(start: int, end: int, molecule: int) -> int:
    """Length of a 1-based inclusive span that may wrap past the origin."""
    if start <= end:
        return end - start + 1
    return molecule - start + 1 + end


@dataclass(frozen=True)
class Cut:
    """One double-strand break.

    Attributes:
        enzyme: name of the enzyme making the cut.
        top: forward-strand nick, as the base 5' of it; the cut falls between
            ``top`` and ``top + 1``.
        bottom: the complementary-strand nick, in the same forward-strand
            frame.
        overhang_length: signed overhang, positive for a 5' overhang, negative
            for a 3' one, zero for a blunt end. Carried rather than recomputed
            from the two nicks, because on a circular molecule a cut near the
            origin has its nicks wrapped independently and their difference no
            longer gives the overhang: an EcoRI cut with nicks at 10 and 2 of a
            12 base circle leaves a four base 5' overhang, not an eight base
            3' one.
        site_start: 1-based start of the recognition site that directed it.
        strand: strand the recognition site was read on.
    """

    enzyme: str
    top: int
    bottom: int
    overhang_length: int
    site_start: int
    strand: str

    @property
    def end_type(self) -> str:
        """Name the end this cut leaves."""
        if self.overhang_length > 0:
            return "5' overhang"
        if self.overhang_length < 0:
            return "3' overhang"
        return "blunt"

    def overhang(self, sequence: str, circular: bool = False) -> str:
        """The single-stranded bases the cut exposes, read 5' to 3'.

        Given as the forward-strand text between the two nicks, which is the
        form in which two ends can be compared: ends anneal when this text and
        the end type both agree. A 5' overhang begins just after the top nick
        and a 3' overhang just after the bottom one, each running for the length
        of the overhang, so the span is correct even where it crosses the origin.
        """
        width = abs(self.overhang_length)
        if not width or not sequence:
            return ""
        anchor = self.top if self.overhang_length > 0 else self.bottom
        length = len(sequence)
        # Walked base by base rather than sliced, because an overhang can be
        # longer than the molecule it sits on. A five-base circle cut by an
        # enzyme that leaves a nine-base overhang wraps past the origin twice,
        # and a single wrapped slice would silently return four bases for a
        # nine-base overhang. Such a molecule is biological nonsense, but the
        # length and the sequence reported for it should at least agree.
        return "".join(sequence[(anchor + step) % length] for step in range(width))


@dataclass(frozen=True)
class Fragment:
    """One product of a digest.

    Attributes:
        start: 1-based first base on the forward strand.
        end: 1-based last base, inclusive; smaller than ``start`` when the
            fragment spans the origin of a circular molecule.
        length: bases between the flanking top-strand nicks.
        left: name of the enzyme cutting the 5' end, or ``end`` for the native
            end of a linear molecule. Two enzymes that break the same
            top-strand bond share one boundary and are both named, joined by
            ``+``.
        right: as ``left``, for the 3' end.
        circular: whether the fragment is itself an uncut circle.
    """

    start: int
    end: int
    length: int
    left: str
    right: str
    circular: bool = False

    @property
    def spans_origin(self) -> bool:
        """Whether the fragment wraps past the origin."""
        return self.end < self.start and not self.circular


@dataclass
class Digest:
    """The result of cutting one sequence with one set of enzymes.

    Attributes:
        name: identifier of the sequence digested.
        length: length of the molecule in bases.
        circular: whether it was treated as circular.
        enzymes: the enzymes applied, in the order given.
        sites: every recognition site found, sorted by position.
        cuts: every double-strand break, sorted by top-strand nick.
        fragments: the products, sorted longest first.
        uncertain: per enzyme, how many positions could not be called because
            the sequence carries an ambiguity code there.
        uncut_sites: per enzyme, how many recognition sites were found that
            produce no cut at all, because every nick would fall beyond the end
            of a linear molecule.
        considered: how many enzymes were applied, when that is more than the
            number listed in ``enzymes``. ``digest --all`` applies the whole
            catalogue and then drops the enzymes with nothing to report, so the
            two differ and both are worth knowing. Zero means every enzyme
            applied is listed.
        partial_sites: per enzyme, how many sites produce some but not all of
            the cuts the enzyme makes. Only an enzyme that cuts on both sides of
            its site, such as BcgI, can do this: near the end of a linear
            molecule one of its two breaks lands and the other does not, so the
            site is cut but does not excise.
    """

    name: str
    length: int
    circular: bool
    enzymes: list[Enzyme]
    sites: list[Site] = field(default_factory=list)
    cuts: list[Cut] = field(default_factory=list)
    fragments: list[Fragment] = field(default_factory=list)
    uncertain: dict[str, int] = field(default_factory=dict)
    uncut_sites: dict[str, int] = field(default_factory=dict)
    partial_sites: dict[str, int] = field(default_factory=dict)
    considered: int = 0

    @property
    def topology(self) -> str:
        """``circular`` or ``linear``."""
        return "circular" if self.circular else "linear"

    def sites_for(self, enzyme: str) -> list[Site]:
        """Every site of one enzyme."""
        return [site for site in self.sites if site.enzyme == enzyme]

    def cuts_for(self, enzyme: str) -> list[Cut]:
        """Every cut made by one enzyme."""
        return [cut for cut in self.cuts if cut.enzyme == enzyme]

    @property
    def cut_counts(self) -> dict[str, int]:
        """How many times each enzyme cuts, including enzymes that do not."""
        counts = {enzyme.name: 0 for enzyme in self.enzymes}
        for cut in self.cuts:
            counts[cut.enzyme] = counts.get(cut.enzyme, 0) + 1
        return counts

    @property
    def is_uncut(self) -> bool:
        """Whether no enzyme cut at all."""
        return not self.cuts

    @property
    def total_length(self) -> int:
        """Sum of the fragment lengths, which must equal ``length``."""
        return sum(fragment.length for fragment in self.fragments)

    @property
    def boundaries(self) -> int:
        """Distinct top-strand nicks, which is what sets the fragment count.

        Lower than the number of cuts when two enzymes break the same
        top-strand bond, which is why a digest can report more cuts than
        fragments without contradicting itself.
        """
        return len({cut.top for cut in self.cuts})

    @property
    def shared_boundaries(self) -> dict[int, tuple[str, ...]]:
        """Top-strand nicks broken by more than one enzyme, by position."""
        sharing: dict[int, list[str]] = {}
        for cut in self.cuts:
            names = sharing.setdefault(cut.top, [])
            if cut.enzyme not in names:
                names.append(cut.enzyme)
        return {
            position: tuple(sorted(names))
            for position, names in sorted(sharing.items())
            if len(names) > 1
        }


def cuts_from_site(sequence_length: int, enzyme: Enzyme, site: Site, circular: bool) -> list[Cut]:
    """Place one site's cuts on the forward strand.

    A site read on the reverse strand has the enzyme working in the opposite
    direction, so its offsets are measured leftwards from the far end of the
    site, and the roles of the two nicks swap: what the enzyme calls its
    top-strand cut falls on this sequence's bottom strand.
    """
    unwrapped_end = site.start + enzyme.length - 1
    placed: list[Cut] = []

    for top_offset, bottom_offset in enzyme.cuts:
        if site.strand == "-":
            top = unwrapped_end - bottom_offset
            bottom = unwrapped_end - top_offset
        else:
            top = site.start + top_offset - 1
            bottom = site.start + bottom_offset - 1
        # Taken before any wrapping, while the two nicks are still in one frame.
        overhang = bottom - top

        if circular:
            top = (top - 1) % sequence_length + 1
            bottom = (bottom - 1) % sequence_length + 1
        elif not (1 <= top < sequence_length and 1 <= bottom < sequence_length):
            # The recognition site is present but at least one nick would fall
            # beyond the end of the molecule, so no double-strand break is
            # made. Counted as an unusable site rather than reported as a cut.
            continue

        placed.append(
            Cut(enzyme.name, top, bottom, overhang, site.start, site.strand)
        )

    return placed


def build_fragments(
    cuts: list[Cut], sequence_length: int, circular: bool
) -> list[Fragment]:
    """Assemble fragments from the top-strand nicks of a set of cuts.

    On a linear molecule ``n`` cuts give ``n + 1`` fragments. On a circular one
    they give ``n``, because the molecule has no ends of its own to contribute;
    a single cut therefore yields one full-length linear fragment rather than
    two pieces. Getting this wrong by one is the classic error in reading a
    plasmid digest, so the two cases are kept visibly separate.
    """
    if not sequence_length:
        return []

    # More than one enzyme can break the same top-strand bond while nicking the
    # bottom strand in different places: AccII and Alw26I both cut after base 3
    # of pUC19, one bluntly and one leaving a four-base overhang. That is one
    # fragment boundary, so the fragment count follows the number of distinct
    # top-strand nicks rather than the number of cuts, but crediting only the
    # first enzyme would drop the other from the fragment it helped make.
    sharing: dict[int, list[str]] = {}
    for cut in sorted(cuts, key=lambda item: (item.top, item.enzyme)):
        names = sharing.setdefault(cut.top, [])
        if cut.enzyme not in names:
            names.append(cut.enzyme)
    boundaries = {position: "+".join(names) for position, names in sharing.items()}

    positions = sorted(boundaries)

    if not positions:
        return [
            Fragment(
                start=1,
                end=sequence_length,
                length=sequence_length,
                left=TERMINUS,
                right=TERMINUS,
                circular=circular,
            )
        ]

    fragments: list[Fragment] = []

    if circular:
        for index, position in enumerate(positions):
            following = positions[(index + 1) % len(positions)]
            start = position % sequence_length + 1
            end = following
            fragments.append(
                Fragment(
                    start=start,
                    end=end,
                    length=span_length(start, end, sequence_length),
                    left=boundaries[position],
                    right=boundaries[following],
                )
            )
    else:
        edges = [0, *positions, sequence_length]
        for index in range(len(edges) - 1):
            start = edges[index] + 1
            end = edges[index + 1]
            left = TERMINUS if index == 0 else boundaries[edges[index]]
            right = TERMINUS if index == len(edges) - 2 else boundaries[edges[index + 1]]
            fragments.append(
                Fragment(
                    start=start,
                    end=end,
                    length=end - start + 1,
                    left=left,
                    right=right,
                )
            )

    fragments.sort(key=lambda fragment: (-fragment.length, fragment.start))
    return fragments


def digest(
    sequence: str,
    enzymes: list[Enzyme],
    circular: bool = False,
    name: str = "sequence",
) -> Digest:
    """Cut a sequence with every enzyme given.

    Args:
        sequence: the molecule, as IUPAC nucleotide codes.
        enzymes: the enzymes to apply together, as in a double digest.
        circular: treat the sequence as a closed circle, so that sites may span
            the origin and ``n`` cuts give ``n`` fragments.
        name: identifier recorded on the result.
    """
    length = len(sequence)
    result = Digest(name=name, length=length, circular=circular, enzymes=list(enzymes))

    # Settled once for the whole sequence rather than once per enzyme: the test
    # walks the sequence, so repeating it across a catalogue-wide digest costs
    # more than finding the sites does.
    ambiguous = has_ambiguity(sequence)

    seen: set[tuple[str, int, int]] = set()
    for enzyme in enzymes:
        sites = find_sites(sequence, enzyme, circular)
        result.sites.extend(sites)

        if ambiguous:
            uncertain = count_uncertain(sequence, enzyme, circular)
            if uncertain:
                result.uncertain[enzyme.name] = uncertain

        uncut = 0
        partial = 0
        for site in sites:
            placed = cuts_from_site(length, enzyme, site, circular)
            if not placed:
                uncut += 1
            elif len(placed) < len(enzyme.cuts):
                partial += 1
            for cut in placed:
                # One physical break can be reached from two site orientations
                # when the recognition sequence is symmetric but the cut is
                # not, so identical breaks are collapsed.
                key = (cut.enzyme, cut.top, cut.bottom)
                if key in seen:
                    continue
                seen.add(key)
                result.cuts.append(cut)
        if uncut:
            result.uncut_sites[enzyme.name] = uncut
        if partial:
            result.partial_sites[enzyme.name] = partial

    result.sites.sort(key=lambda site: (site.start, site.enzyme, site.strand))
    result.cuts.sort(key=lambda cut: (cut.top, cut.enzyme))
    result.fragments = build_fragments(result.cuts, length, circular)
    return result
