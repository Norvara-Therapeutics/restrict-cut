"""Number and sequence formatting shared by the terminal output and the report.

Only human-facing output is formatted here. The TSV and JSON outputs carry
unrounded values and full sequences, so a figure another program reads is never
the rounded one.
"""

from __future__ import annotations


def thousands(value: int) -> str:
    """Group an integer in thousands."""
    return f"{value:,}"


def plural(count: int, singular: str, suffix: str = "s") -> str:
    """Count something, agreeing the noun with the number."""
    return f"{thousands(count)} {singular}{'' if count == 1 else suffix}"


def number(value: float | None, places: int = 2, missing: str = "n/a") -> str:
    """Format an optional number to a fixed number of decimal places."""
    if value is None:
        return missing
    return f"{value:,.{places}f}"


def base_pairs(value: int) -> str:
    """Write a fragment length the way a gel is read.

    Lengths stay in bases up to ten kilobases and are given in kilobases above
    that, where a gel cannot resolve the last digits anyway. The exact value is
    always available in the TSV and JSON outputs.
    """
    if value < 10_000:
        return f"{value:,} bp"
    return f"{value / 1000:,.1f} kb"


def overhang_display(overhang: str, end_type: str) -> str:
    """Show an overhang with the strand it sits on.

    A 5' overhang is written ``5'-GATC``, a 3' overhang ``3'-TGCA``, and a
    blunt end ``blunt``, so that two ends can be told apart at a glance even
    when they expose the same four bases.
    """
    if not overhang:
        return "blunt"
    if end_type.startswith("5"):
        return f"5'-{overhang}"
    return f"3'-{overhang}"


def wrap_sequence(sequence: str, width: int = 60) -> list[str]:
    """Break a sequence into fixed-width lines."""
    return [sequence[index: index + width] for index in range(0, len(sequence), width)]
