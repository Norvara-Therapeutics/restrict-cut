"""Ordering fragments by size and judging whether they separate."""

import pytest

from restrict_cut.gel import (
    DEFAULT_RESOLUTION,
    lanes,
    render,
    summary,
    unresolved_pairs,
)


def test_bands_are_ordered_largest_first():
    bands = lanes([500, 3000, 1200])
    assert [band.size for band in bands] == [3000, 1200, 500]


def test_fragments_of_equal_size_share_one_band():
    bands = lanes([1000, 1000, 500])
    assert len(bands) == 2
    doubled = next(band for band in bands if band.size == 1000)
    assert doubled.count == 2


def test_offsets_run_from_largest_to_smallest():
    bands = lanes([100, 1000, 10000])
    offsets = {band.size: band.offset for band in bands}
    assert offsets[10000] == pytest.approx(0.0)
    assert offsets[100] == pytest.approx(1.0)
    assert 0.0 < offsets[1000] < 1.0


def test_offsets_are_logarithmic():
    """1000 sits halfway between 100 and 10000 on a log axis."""
    bands = lanes([100, 1000, 10000])
    offsets = {band.size: band.offset for band in bands}
    assert offsets[1000] == pytest.approx(0.5)


def test_a_single_band_sits_at_the_top():
    bands = lanes([1234])
    assert bands[0].offset == pytest.approx(0.0)
    assert bands[0].is_resolved


def test_close_sizes_are_called_unresolved():
    assert unresolved_pairs([1000, 1020]) == [(1000, 1020)]
    assert unresolved_pairs([1000, 2000]) == []


def test_resolution_is_relative_not_absolute():
    """A 40 bp gap is plain at 200 bp and invisible at 10 kb."""
    assert unresolved_pairs([200, 240]) == []
    assert unresolved_pairs([10000, 10040]) == [(10000, 10040)]


def test_resolution_threshold_is_configurable():
    assert unresolved_pairs([1000, 1200], resolution=10.0) == []
    assert unresolved_pairs([1000, 1200], resolution=30.0) == [(1000, 1200)]


def test_zero_resolution_calls_everything_resolved():
    assert unresolved_pairs([1000, 1001], resolution=0.0) == []


def test_unresolved_is_recorded_on_both_bands():
    bands = lanes([1000, 1020])
    for band in bands:
        assert not band.is_resolved
        assert band.unresolved_with


def test_summary_reports_doubled_and_unresolved_bands():
    text = summary(lanes([3000, 3000, 1000, 1020]))
    assert "equal length" in text
    assert "within 5%" in text


def test_summary_of_a_clean_digest():
    text = summary(lanes([4000, 2000, 500]))
    assert "all separated" in text


def test_summary_of_nothing():
    assert summary([]) == "no fragments"


def test_render_gives_one_row_per_band():
    bands = lanes([23130, 9416, 6682, 4361, 2322, 2027, 564])
    drawn = "\n".join(render(bands))
    for band in bands:
        assert f"{band.size:,}" in drawn


def test_render_separates_bands_that_resolve():
    """Two bands more than the threshold apart must not share a row."""
    bands = lanes([2322, 2027])
    rows = [line for line in render(bands) if line.strip()]
    assert len(rows) == 2


def test_render_marks_doubled_and_unresolved_bands():
    drawn = "\n".join(render(lanes([3000, 3000, 1000, 1020])))
    assert "x2" in drawn
    assert "~" in drawn


def test_render_of_nothing_is_empty():
    assert render([]) == []


def test_render_with_a_ladder_labels_both_lanes():
    sample = lanes([1500, 700])
    ladder = lanes([2000, 1000, 500])
    drawn = "\n".join(render(sample, ladder=ladder, ladder_name="given sizes"))
    assert "ladder" in drawn
    assert "sample" in drawn
    assert "1,500" in drawn
    assert "2,000" in drawn


def test_default_resolution_is_five_percent():
    assert DEFAULT_RESOLUTION == 5.0


def test_no_band_is_dropped_when_rows_are_free():
    """Sizes 100, 99, 98, 1 once lost the 98 band to a row index of -1."""
    bands = lanes([100, 99, 98, 1])
    drawn = "\n".join(render(bands))
    for size in (100, 99, 98, 1):
        assert f"{size:,}" in drawn, size


@pytest.mark.parametrize(
    "sizes",
    [
        [100, 99, 98, 1],
        [100, 99, 98, 97, 96, 95],
        [5000, 4999, 4998, 10, 9, 8],
        list(range(100, 130)),
    ],
)
def test_every_band_is_drawn_however_crowded(sizes):
    bands = lanes(sizes)
    drawn = "\n".join(render(bands))
    for band in bands:
        assert f"{band.size:,}" in drawn, band.size


def test_rows_are_distinct_for_distinct_bands():
    from restrict_cut.gel import _rows_for

    bands = lanes([100, 99, 98, 1])
    rows = _rows_for(bands, 14)
    assert len(set(rows.values())) == len(bands)
    assert all(row >= 0 for row in rows.values())


def test_sample_and_ladder_share_one_size_axis():
    """A sample band level with a ladder band must be the same size."""
    sample = lanes([1000, 500])
    ladder = lanes([2000, 1000, 500])
    drawn = render(sample, ladder=ladder, ladder_name="given sizes")

    def row_of(text, lane):
        for index, line in enumerate(drawn):
            if lane == "ladder" and line.strip().startswith(text):
                return index
            if lane == "sample" and line.rstrip().endswith(text) and "---" in line:
                return index
        raise AssertionError(f"{text} not found in {lane}")

    assert row_of("1,000", "sample") == row_of("1,000", "ladder")
    assert row_of("500", "sample") == row_of("500", "ladder")


def test_a_ladder_band_larger_than_the_sample_sits_above_it():
    sample = lanes([1000, 500])
    ladder = lanes([2000, 1000, 500])
    drawn = render(sample, ladder=ladder)
    top = next(i for i, line in enumerate(drawn) if "2,000" in line)
    aligned = next(i for i, line in enumerate(drawn) if line.rstrip().endswith("1,000"))
    assert top < aligned
