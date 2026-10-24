"""Formatting helpers."""

import pytest

from restrict_cut.text import (
    base_pairs,
    number,
    overhang_display,
    plural,
    thousands,
    wrap_sequence,
)


def test_thousands_groups_digits():
    assert thousands(2686) == "2,686"
    assert thousands(7) == "7"


def test_plural_agrees_with_the_count():
    assert plural(1, "cut") == "1 cut"
    assert plural(2, "cut") == "2 cuts"
    assert plural(0, "cut") == "0 cuts"


def test_plural_takes_an_irregular_suffix():
    assert plural(2, "analys", "es") == "2 analyses"


def test_number_handles_missing_values():
    assert number(None) == "n/a"
    assert number(1.5, 1) == "1.5"


def test_base_pairs_switches_to_kilobases_above_ten_thousand():
    assert base_pairs(2686) == "2,686 bp"
    assert base_pairs(9999) == "9,999 bp"
    assert base_pairs(23130) == "23.1 kb"


@pytest.mark.parametrize(
    "overhang,end_type,expected",
    [
        ("GATC", "5' overhang", "5'-GATC"),
        ("TGCA", "3' overhang", "3'-TGCA"),
        ("", "blunt", "blunt"),
    ],
)
def test_overhang_display_names_the_strand(overhang, end_type, expected):
    assert overhang_display(overhang, end_type) == expected


def test_overhang_display_separates_identical_bases_on_different_strands():
    assert overhang_display("GATC", "5' overhang") != overhang_display(
        "GATC", "3' overhang"
    )


def test_wrap_sequence_splits_at_the_given_width():
    assert wrap_sequence("ACGTACGTAC", 4) == ["ACGT", "ACGT", "AC"]
    assert wrap_sequence("") == []
