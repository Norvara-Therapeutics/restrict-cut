"""Laying fragments out by size, the way a gel would separate them.

This is a size ladder drawn on a logarithmic axis, not a simulation of
electrophoresis. Nothing here models agarose percentage, buffer, voltage or
run time, and no migration distance in millimetres is claimed. What it does
give is the thing a map alone will not: an ordering by size with the bands that
will *not* come apart marked, because a digest that looks informative in a
table of fragment lengths can be unreadable on a gel.

Two fragments are called unresolved when their sizes differ by less than a set
percentage. That is a rule of thumb rather than a measurement, and the
percentage is the caller's to set.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: Default percentage difference below which two bands are called unresolved.
DEFAULT_RESOLUTION = 5.0

@dataclass(frozen=True)
class Band:
    """One band on the ladder.

    Attributes:
        size: fragment length in bases.
        count: how many fragments of this size are present.
        labels: what produced them, for annotating the band.
        offset: position on the axis, 0.0 at the largest size and 1.0 at the
            smallest, spaced by the logarithm of the size.
        unresolved_with: sizes of other bands too close to separate from this
            one.
    """

    size: int
    count: int
    labels: tuple[str, ...]
    offset: float
    unresolved_with: tuple[int, ...] = ()

    @property
    def is_resolved(self) -> bool:
        """Whether this band separates from every other band."""
        return not self.unresolved_with


def _offsets(sizes: list[int]) -> dict[int, float]:
    """Place sizes on a logarithmic axis between 0.0 and 1.0."""
    if not sizes:
        return {}
    largest, smallest = max(sizes), min(sizes)
    if largest == smallest:
        return {size: 0.0 for size in sizes}
    high, low = math.log10(largest), math.log10(smallest)
    return {size: (high - math.log10(size)) / (high - low) for size in set(sizes)}


def unresolved_pairs(
    sizes: list[int], resolution: float = DEFAULT_RESOLUTION
) -> list[tuple[int, int]]:
    """Pairs of distinct sizes that differ by less than ``resolution`` percent.

    The comparison is relative, taken against the smaller of the two, because a
    40 bp difference is plain at 200 bp and invisible at 10 kb.
    """
    if resolution <= 0:
        return []
    ordered = sorted(set(sizes))
    close = []
    for index, size in enumerate(ordered):
        for other in ordered[index + 1:]:
            if (other - size) / size * 100.0 >= resolution:
                break
            close.append((size, other))
    return close


def lanes(
    fragments: list[int],
    labels: dict[int, tuple[str, ...]] | None = None,
    resolution: float = DEFAULT_RESOLUTION,
) -> list[Band]:
    """Collect fragment sizes into bands, largest first.

    Fragments of identical size share one band, because they would run as one,
    and the band records how many fragments it stands for. A digest that gives
    two fragments of the same length shows a single band of twice the
    intensity, which is a reading a table of lengths will not warn about.
    """
    if not fragments:
        return []

    tally: dict[int, int] = {}
    for size in fragments:
        tally[size] = tally.get(size, 0) + 1

    offsets = _offsets(list(tally))
    close = unresolved_pairs(list(tally), resolution)
    neighbours: dict[int, set[int]] = {size: set() for size in tally}
    for first, second in close:
        neighbours[first].add(second)
        neighbours[second].add(first)

    bands = [
        Band(
            size=size,
            count=count,
            labels=(labels or {}).get(size, ()),
            offset=offsets[size],
            unresolved_with=tuple(sorted(neighbours[size])),
        )
        for size, count in tally.items()
    ]
    bands.sort(key=lambda band: -band.size)
    return bands


def _nearest_free(preferred: int, rows: int, used: set[int]) -> int:
    """The free row closest to the one a band would prefer.

    Searching outward in both directions is guaranteed to find a row whenever
    there are at least as many rows as bands, which :func:`render` ensures. An
    earlier version gave up after one step in each direction and returned -1,
    which silently dropped a band from the lane while rows sat empty.
    """
    for step in range(rows):
        for option in ((preferred,) if not step else (preferred + step, preferred - step)):
            if 0 <= option < rows and option not in used:
                return option
    raise ValueError(f"no free row among {rows} for {len(used)} already placed")


def _rows_for(
    bands: list[Band], rows: int, offsets: dict[int, float] | None = None
) -> dict[int, int]:
    """Give each band a distinct row, as near its size position as possible.

    Rounding the logarithmic offset alone puts two bands of different size on
    one row of a short lane, which would hide a difference the digest does
    resolve. Bands are placed largest first and nudged to the nearest free row.

    ``offsets`` overrides each band's own offset, which is how a sample and a
    ladder are made to share one axis: a band's position has to mean the same
    thing in both lanes, or reading one against the other gives the wrong size.
    """
    taken: dict[int, int] = {}
    used: set[int] = set()
    for band in sorted(bands, key=lambda item: -item.size):
        offset = band.offset if offsets is None else offsets[band.size]
        preferred = min(rows - 1, max(0, round(offset * (rows - 1))))
        row = _nearest_free(preferred, rows, used)
        used.add(row)
        taken[band.size] = row
    return taken


def render(
    bands: list[Band],
    width: int = 28,
    ladder: list[Band] | None = None,
    ladder_name: str = "",
) -> list[str]:
    """Draw the bands as text, largest at the top.

    Each band sits on a logarithmic size axis and is drawn as one row. When a
    ladder is given, both lanes are placed on a single axis spanning every size
    in either of them, so a sample band level with a ladder band really is the
    same size. A band that does not separate from a neighbour is marked with a
    tilde.
    """
    if not bands:
        return []

    everything = [band.size for band in bands] + [band.size for band in ladder or []]
    shared = _offsets(everything)
    rows = max(14, min(60, (len(bands) + len(ladder or [])) * 2 + 4))
    rows = max(rows, len(bands), len(ladder or []))

    sample_lane = [""] * rows
    ladder_lane = [""] * rows

    by_size = {band.size: band for band in bands}
    for size, row in _rows_for(bands, rows, shared).items():
        found = by_size[size]
        text = f"{found.size:,}"
        if found.count > 1:
            text += f" x{found.count}"
        if not found.is_resolved:
            text += " ~"
        sample_lane[row] = text

    if ladder:
        for size, row in _rows_for(ladder, rows, shared).items():
            ladder_lane[row] = f"{size:,}"

    lines = []
    if ladder:
        lines.append(f"{'ladder':>{width}}  sample")
    for row in range(rows):
        left = f"{ladder_lane[row]:>{width}}" if ladder else ""
        marker = "---" if sample_lane[row] else "   "
        lines.append(f"{left}  {marker} {sample_lane[row]}".rstrip())
    if ladder and ladder_name:
        lines.append(f"{'':>{width}}  ladder: {ladder_name}")
    return lines


def summary(bands: list[Band], resolution: float = DEFAULT_RESOLUTION) -> str:
    """One line on whether the digest can be read off a gel."""
    if not bands:
        return "no fragments"
    unresolved = [band for band in bands if not band.is_resolved]
    doubled = [band for band in bands if band.count > 1]
    parts = [f"{len(bands)} band{'' if len(bands) == 1 else 's'}"]
    if doubled:
        parts.append(
            f"{len(doubled)} carrying fragments of equal length"
        )
    if unresolved:
        parts.append(
            f"{len(unresolved)} within {resolution:g}% of another band"
        )
    else:
        parts.append(f"all separated by more than {resolution:g}%")
    return "; ".join(parts)
