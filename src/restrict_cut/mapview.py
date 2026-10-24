"""A text map of where the cuts fall.

The map is drawn to scale along the molecule, so the spacing between marks
carries information: a cluster of sites in a polylinker looks like a cluster.
Where several cuts land on the same column of a narrow terminal they are
counted rather than overprinted, because a map that silently loses a site is
worse than one that admits it ran out of room.
"""

from __future__ import annotations

from dataclasses import dataclass

from restrict_cut.text import thousands

#: Character marking a single cut on the ruler.
MARK = "|"

#: Character marking a column holding more than one cut.
CROWDED = "*"


@dataclass(frozen=True)
class Column:
    """One column of the map ruler.

    Attributes:
        index: 0-based column on the ruler.
        position: representative cut position for the column.
        enzymes: names of the enzymes cutting in this column.
    """

    index: int
    position: int
    enzymes: tuple[str, ...]
    positions: tuple[int, ...] = ()

    @property
    def crowded(self) -> bool:
        """Whether more than one distinct cut shares the column."""
        return len(self.positions) > 1


def columns(
    cuts: list[tuple[int, str]], length: int, width: int
) -> list[Column]:
    """Assign cuts to ruler columns, scaled to the molecule length."""
    if length <= 0 or width <= 0:
        return []
    grouped: dict[int, list[tuple[int, str]]] = {}
    for position, enzyme in sorted(cuts):
        index = min(width - 1, int((position - 1) / length * width))
        grouped.setdefault(index, []).append((position, enzyme))
    return [
        Column(
            index=index,
            position=entries[0][0],
            enzymes=tuple(dict.fromkeys(enzyme for _, enzyme in entries)),
            positions=tuple(dict.fromkeys(position for position, _ in entries)),
        )
        for index, entries in sorted(grouped.items())
    ]


def ruler(length: int, width: int, circular: bool) -> list[str]:
    """Draw the molecule as a scaled line with its end positions.

    The line is never drawn wider than the width asked for, which at a width of
    one or two leaves room for the end caps and nothing between them.
    """
    caps, fill = ("o", "=") if circular else ("[", "-")
    if width < 2:
        line = caps * max(0, width)
    else:
        line = caps + fill * (width - 2) + caps
    labels = f"1{' ' * max(1, width - len(thousands(length)) - 1)}{thousands(length)}"
    return [line, labels]


def render(
    cuts: list[tuple[int, str]],
    length: int,
    width: int = 64,
    circular: bool = False,
) -> list[str]:
    """Draw a scaled map of the cut positions.

    Each cut gets a mark on the ruler and a labelled line beneath giving its
    enzyme and position. Columns holding more than one cut are marked
    differently and list every enzyme that falls there.
    """
    if length <= 0:
        return []
    placed = columns(cuts, length, width)
    marks = [" "] * width
    for column in placed:
        marks[column.index] = CROWDED if column.crowded else MARK

    lines = ["".join(marks).rstrip()]
    lines.extend(ruler(length, width, circular))

    for column in placed:
        names = ", ".join(column.enzymes)
        marker = CROWDED if column.crowded else MARK
        positions = ", ".join(thousands(position) for position in column.positions)
        lines.append(f"  {marker} {names} at {positions}")
    return lines
