"""The generator that builds the enzyme table from REBASE.

The package's documentation makes a claim about its own data: that no cut
offset ships without agreeing with a second, independent statement of the same
fact. That claim is only as good as this script's refusal to write, so the
refusal is tested rather than assumed - including the case that matters most,
an entry the second source says nothing about. Silence is not agreement.
"""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "build_enzyme_table.py"


@pytest.fixture(scope="module")
def builder():
    if not SCRIPT.is_file():
        pytest.skip("the generator is not present in this build")
    spec = importlib.util.spec_from_file_location("build_enzyme_table", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_source(directory: Path, emboss: str, gtype2: str) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "link_emboss_e").write_text(f"# REBASE version 999\n{emboss}")
    (directory / "link_gtype2").write_text(gtype2)
    (directory / "link_emboss_s").write_text("# codes\nN New England Biolabs\n")
    return directory


ECORI_EMBOSS = "EcoRI\tGAATTC\t6\t2\t0\t1\t5\t0\t0\n"
ECORI_GTYPE2 = "<1>EcoRI\n<2>\n<3>G^AATTC\n<4>\n<5>N\n<6>1\n\n"


def test_a_sound_source_is_written(builder, tmp_path):
    source = write_source(tmp_path / "ok", ECORI_EMBOSS, ECORI_GTYPE2)
    out = tmp_path / "table.py"
    assert builder.main(["--source", str(source), "--output", str(out)]) == 0
    assert '"EcoRI"' in out.read_text()


def test_an_entry_the_second_source_cannot_confirm_is_refused(builder, tmp_path, capsys):
    """The failure that matters: not a disagreement, but silence."""
    source = write_source(
        tmp_path / "unverifiable",
        ECORI_EMBOSS + "FakeI\tGGGCCC\t6\t2\t0\t2\t4\t0\t0\n",
        ECORI_GTYPE2 + "<1>FakeI\n<2>\n<3>\n<4>\n<5>N\n<6>2\n\n",
    )
    out = tmp_path / "table.py"
    assert builder.main(["--source", str(source), "--output", str(out)]) == 1
    assert not out.exists()
    printed = capsys.readouterr().out
    assert "FakeI" in printed
    assert "could not be cross-checked" in printed


def test_an_entry_that_disagrees_is_refused(builder, tmp_path, capsys):
    source = write_source(
        tmp_path / "disagreeing",
        "EcoRI\tGAATTC\t6\t2\t0\t2\t4\t0\t0\n",
        ECORI_GTYPE2,
    )
    out = tmp_path / "table.py"
    assert builder.main(["--source", str(source), "--output", str(out)]) == 1
    assert not out.exists()
    assert "disagree" in capsys.readouterr().out


def test_nothing_is_written_when_the_gate_refuses(builder, tmp_path):
    """A refusal must not leave a partial table where a good one was."""
    out = tmp_path / "table.py"
    out.write_text("# the previous, trusted table\n")
    source = write_source(
        tmp_path / "bad",
        ECORI_EMBOSS + "FakeI\tGGGCCC\t6\t2\t0\t2\t4\t0\t0\n",
        ECORI_GTYPE2 + "<1>FakeI\n<2>\n<3>\n<4>\n<5>N\n<6>2\n\n",
    )
    assert builder.main(["--source", str(source), "--output", str(out)]) == 1
    assert out.read_text() == "# the previous, trusted table\n"


def test_a_missing_source_file_is_reported(builder, tmp_path, capsys):
    bare = tmp_path / "empty"
    bare.mkdir()
    assert builder.main(["--source", str(bare), "--output", str(tmp_path / "t.py")]) == 1
    assert "missing" in capsys.readouterr().err


@pytest.mark.parametrize(
    "notation,pattern,cuts",
    [
        ("G^AATTC", "GAATTC", ((1, 5),)),
        ("CTGCA^G", "CTGCAG", ((5, 1),)),
        ("GGTCTC(1/5)", "GGTCTC", ((7, 11),)),
        ("CCGC(-3/-1)", "CCGC", ((1, 3),)),
        ("(10/12)CGANNNNNNTGC(12/10)", "CGANNNNNNTGC", ((-10, -12), (24, 22))),
    ],
)
def test_notation_parser_reads_every_rebase_form(builder, notation, pattern, cuts):
    read_pattern, read_cuts = builder.parse_notation(notation)
    assert read_pattern == pattern
    assert read_cuts == cuts


def test_notation_parser_handles_the_padded_trailing_caret(builder):
    """REBASE spells some downstream cuts out in Ns and omits the bottom one."""
    pattern, cuts = builder.parse_notation("TGGCCANNNNNNNNNNN^")
    assert pattern == "TGGCCA"
    assert cuts == ((17, None),)


def test_emboss_upstream_coordinates_are_corrected(builder):
    """EMBOSS writes a cut ten bases upstream as -11, not -10."""
    assert builder._offset(-11) == -10
    assert builder._offset(1) == 1
    assert builder._offset(0) == 1


def test_the_shipped_table_came_from_a_gate_that_passed():
    """Whatever the generator does now, the shipped data must still be sound."""
    from restrict_cut import enzymes

    assert enzymes.validate_catalogue() == []
