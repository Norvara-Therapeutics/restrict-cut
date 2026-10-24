"""The scaled text map."""

from restrict_cut.mapview import CROWDED, MARK, columns, render


def test_columns_scale_to_the_molecule():
    placed = columns([(1, "A"), (500, "B"), (1000, "C")], 1000, 10)
    assert [column.index for column in placed] == [0, 4, 9]


def test_cuts_sharing_a_column_are_grouped():
    placed = columns([(10, "A"), (12, "B")], 1000, 10)
    assert len(placed) == 1
    assert placed[0].enzymes == ("A", "B")
    assert placed[0].crowded


def test_a_single_cut_is_not_crowded():
    placed = columns([(10, "A")], 1000, 10)
    assert not placed[0].crowded


def test_two_enzymes_cutting_the_same_base_are_one_position():
    placed = columns([(10, "A"), (10, "B")], 1000, 10)
    assert placed[0].positions == (10,)
    assert placed[0].enzymes == ("A", "B")


def test_render_marks_every_cut():
    lines = render([(100, "EcoRI"), (900, "PstI")], 1000, width=40)
    assert lines[0].count(MARK) == 2
    assert "EcoRI at 100" in "\n".join(lines)
    assert "PstI at 900" in "\n".join(lines)


def test_crowded_columns_are_marked_differently():
    lines = render([(100, "EcoRI"), (105, "PstI")], 10000, width=20)
    assert CROWDED in lines[0]
    assert "EcoRI, PstI" in "\n".join(lines)


def test_circular_and_linear_rulers_differ():
    circular = render([(100, "EcoRI")], 1000, width=30, circular=True)
    linear = render([(100, "EcoRI")], 1000, width=30, circular=False)
    assert circular[1] != linear[1]
    assert circular[1].startswith("o")
    assert linear[1].startswith("[")


def test_ruler_shows_the_length():
    lines = render([(100, "EcoRI")], 2686, width=40, circular=True)
    assert "2,686" in lines[2]


def test_no_cuts_still_draws_the_molecule():
    lines = render([], 1000, width=30)
    assert len(lines) == 3


def test_zero_length_renders_nothing():
    assert render([(1, "EcoRI")], 0, width=30) == []


def test_every_cut_appears_even_when_the_map_is_narrow():
    """A narrow map may crowd marks, but must not lose an enzyme."""
    cuts = [(position, f"E{index}") for index, position in enumerate(range(1, 200, 7))]
    text = "\n".join(render(cuts, 1000, width=8))
    for _, name in cuts:
        assert name in text
