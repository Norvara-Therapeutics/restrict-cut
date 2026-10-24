"""The command-line interface."""

import io
import json
import os

import pytest

from restrict_cut import __version__
from restrict_cut.cli import main
from restrict_cut.enzymes import all_enzymes


def run(argv):
    """Run the CLI, returning (exit code, stdout)."""
    out = io.StringIO()
    code = main(argv, stdout=out)
    return code, out.getvalue()


def fasta(tmp_path, sequence, name="seq", filename="input.fasta"):
    path = tmp_path / filename
    path.write_text(f">{name}\n{sequence}\n")
    return str(path)


PLASMID = "GAATTCAAAAAAAAAAGGATCCAAAAAAAAAACTGCAGAAAAAAAAAA"


def test_version_is_reported():
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])
    assert exit_info.value.code == 0


def test_no_command_is_an_error():
    with pytest.raises(SystemExit):
        main([])


def test_digest_table_names_the_enzymes(tmp_path):
    code, text = run(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI,BamHI"])
    assert code == 0
    assert "EcoRI" in text and "BamHI" in text
    assert "G^AATTC" in text


def test_digest_reports_fragments(tmp_path):
    _, text = run(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI"])
    assert "fragments" in text


def test_digest_circular_changes_the_fragment_count(tmp_path):
    path = fasta(tmp_path, PLASMID)
    _, linear = run(["digest", "-i", path, "-e", "EcoRI", "-f", "tsv"])
    _, circular = run(["digest", "-i", path, "-e", "EcoRI", "-f", "tsv", "--circular"])
    assert len(linear.splitlines()) == 3
    assert len(circular.splitlines()) == 2


def test_digest_tsv_has_a_header_and_one_row_per_fragment(tmp_path):
    _, text = run(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI", "-f", "tsv"])
    lines = text.splitlines()
    assert lines[0].split("\t")[:3] == ["sequence", "length", "start"]
    assert len(lines) == 3


def test_digest_json_is_valid_and_complete(tmp_path):
    _, text = run(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI", "-f", "json"])
    payload = json.loads(text)
    assert payload["tool"] == "restrict-cut"
    assert payload["version"] == __version__
    assert payload["rebase_release"]
    assert payload["sequence"]["length"] == len(PLASMID)
    assert payload["total_fragment_length"] == len(PLASMID)
    assert len(payload["cuts"]) == 1
    assert payload["cuts"][0]["overhang"] == "AATT"


def test_digest_none_format_prints_nothing(tmp_path):
    _, text = run(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI", "-f", "none"])
    assert text == ""


def test_digest_writes_a_fragment_file(tmp_path):
    out = tmp_path / "fragments.tsv"
    run([
        "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
        "-o", str(out), "-f", "none",
    ])
    assert out.read_text().splitlines()[0].startswith("sequence\t")


def test_digest_writes_a_cut_file(tmp_path):
    out = tmp_path / "cuts.tsv"
    run([
        "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
        "--cut-table", str(out), "-f", "none",
    ])
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[:4] == ["sequence", "enzyme", "top_nick", "bottom_nick"]
    assert len(lines) == 2


def test_digest_writes_an_html_report(tmp_path):
    out = tmp_path / "report.html"
    run([
        "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI,BamHI",
        "--report", str(out), "-f", "none",
    ])
    html = out.read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert "EcoRI" in html


def test_digest_all_applies_the_catalogue(tmp_path):
    _, text = run([
        "digest", "-i", fasta(tmp_path, PLASMID), "--all", "--min-site", "6",
        "--available", "-f", "tsv",
    ])
    assert len(text.splitlines()) > 2


def test_digest_requires_enzymes_or_all(tmp_path):
    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", fasta(tmp_path, PLASMID)])
    assert exit_info.value.code == 2


def test_digest_refuses_enzyme_and_all_together(tmp_path):
    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI", "--all"])
    assert exit_info.value.code == 2


def test_unknown_enzyme_is_reported(tmp_path, capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "NotAnEnzyme"])
    assert exit_info.value.code == 2
    assert "unknown enzyme" in capsys.readouterr().err


def test_missing_file_is_reported(tmp_path, capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", str(tmp_path / "absent.fasta"), "-e", "EcoRI"])
    assert exit_info.value.code == 2
    assert "no such file" in capsys.readouterr().err


def test_enzymes_may_be_repeated_or_comma_separated(tmp_path):
    path = fasta(tmp_path, PLASMID)
    _, joined = run(["digest", "-i", path, "-e", "EcoRI,BamHI", "-f", "tsv"])
    _, repeated = run(["digest", "-i", path, "-e", "EcoRI", "-e", "BamHI", "-f", "tsv"])
    assert joined == repeated


def test_methylation_can_be_turned_off(tmp_path):
    path = fasta(tmp_path, PLASMID)
    _, with_check = run(["digest", "-i", path, "-e", "BamHI"])
    _, without = run(["digest", "-i", path, "-e", "BamHI", "--methylation", "none"])
    assert "methylation" in with_check
    assert "methylation" not in without


def test_methylation_reports_the_dam_site_inside_bamhi(tmp_path):
    _, text = run(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "BamHI"])
    assert "dam" in text
    assert "not the same as blocked" in text


def test_unknown_methylation_system_is_reported(tmp_path, capsys):
    with pytest.raises(SystemExit):
        main(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
              "--methylation", "ecoki"])
    assert "unknown methylation system" in capsys.readouterr().err


def test_map_can_be_suppressed(tmp_path):
    path = fasta(tmp_path, PLASMID)
    _, with_map = run(["digest", "-i", path, "-e", "EcoRI"])
    _, without = run(["digest", "-i", path, "-e", "EcoRI", "--no-map"])
    assert len(without) < len(with_map)


def test_ladder_sizes_are_shown(tmp_path):
    _, text = run([
        "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
        "--ladder", "100,50,25",
    ])
    assert "ladder" in text


def test_bad_ladder_is_reported(tmp_path, capsys):
    with pytest.raises(SystemExit):
        main(["digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
              "--ladder", "abc"])
    assert "not a fragment size" in capsys.readouterr().err


def test_select_finds_single_cutters(tmp_path):
    _, text = run([
        "select", "-i", fasta(tmp_path, PLASMID), "--cuts", "1", "--min-site", "6",
    ])
    assert "EcoRI" in text
    assert "meet the selection" in text


def test_select_zero_cuts_finds_non_cutters(tmp_path):
    _, text = run([
        "select", "-i", fasta(tmp_path, PLASMID), "--cuts", "0", "-f", "tsv",
        "--min-site", "6", "--available",
    ])
    rows = text.splitlines()[1:]
    assert rows
    assert all(row.split("\t")[2] == "0" for row in rows)


def test_select_protect_excludes_enzymes_cutting_the_other_file(tmp_path):
    target = fasta(tmp_path, PLASMID, filename="target.fasta")
    insert = fasta(tmp_path, "GAATTCGGG", filename="insert.fasta")
    _, unprotected = run(["select", "-i", target, "--cuts", "1", "-f", "tsv"])
    _, protected = run([
        "select", "-i", target, "--cuts", "1", "-f", "tsv", "-p", insert,
    ])
    assert "EcoRI" in unprotected
    assert "EcoRI" not in protected


def test_select_json_lists_the_protected_sequences(tmp_path):
    target = fasta(tmp_path, PLASMID, filename="target.fasta")
    insert = fasta(tmp_path, "GAATTCGGG", name="myinsert", filename="insert.fasta")
    _, text = run(["select", "-i", target, "-f", "json", "-p", insert])
    payload = json.loads(text)
    assert payload["protected"] == ["myinsert"]
    assert payload["considered"] > 0


def test_select_limit_caps_the_output(tmp_path):
    _, text = run([
        "select", "-i", fasta(tmp_path, PLASMID), "--cuts", "1", "-f", "tsv",
        "--limit", "3",
    ])
    assert len(text.splitlines()) == 4


def test_select_reports_when_nothing_matches(tmp_path):
    _, text = run([
        "select", "-i", fasta(tmp_path, "AAAAAAAAAAAA"), "--cuts", "5",
    ])
    assert "no enzyme" in text


def test_enzymes_lists_the_catalogue():
    _, text = run(["enzymes", "--available", "--min-site", "8"])
    assert "ENZYME" in text
    assert "NotI" in text


def test_enzymes_describes_one_enzyme():
    _, text = run(["enzymes", "BsaI"])
    assert "GGTCTC(1/5)" in text
    assert "set by the DNA" in text
    assert "dam overlap" in text


def test_enzyme_detail_lists_isoschizomers():
    _, text = run(["enzymes", "SmaI"])
    assert "same site, other cut" in text
    assert "XmaI" in text


def test_enzymes_filters_by_site():
    _, text = run(["enzymes", "--site", "CCCGGG", "-f", "tsv"])
    rows = text.splitlines()[1:]
    assert rows
    assert all(row.split("\t")[1] == "CCCGGG" for row in rows)


def test_enzymes_rejects_a_non_iupac_site():
    with pytest.raises(SystemExit):
        main(["enzymes", "--site", "GAATTZ"])


def test_enzymes_compatible_groups_ligatable_ends():
    _, text = run(["enzymes", "--compatible", "--available"])
    assert "BamHI" in text
    assert "ligate" in text


def test_enzymes_json_carries_the_cut_offsets():
    _, text = run(["enzymes", "--site", "GAATTC", "-f", "json"])
    payload = json.loads(text)
    ecori = next(item for item in payload["enzymes"] if item["name"] == "EcoRI")
    assert ecori["cuts"] == [[1, 5]]
    assert ecori["overhang"] == "AATT"


def test_unknown_supplier_is_reported(capsys):
    with pytest.raises(SystemExit):
        main(["enzymes", "--supplier", "ZZ"])
    assert "unknown supplier code" in capsys.readouterr().err


def test_empty_sequence_is_reported(tmp_path, capsys):
    path = tmp_path / "empty.fasta"
    path.write_text(">s\n\n")
    with pytest.raises(SystemExit):
        main(["digest", "-i", str(path), "-e", "EcoRI"])
    assert "no sequence" in capsys.readouterr().err


def test_output_refuses_to_overwrite_its_own_input(tmp_path, capsys):
    """The sequence is already in memory by then, so this would destroy it."""
    path = fasta(tmp_path, PLASMID)
    before = (tmp_path / "input.fasta").read_bytes()
    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", path, "-e", "EcoRI", "-o", path, "-f", "none"])
    assert exit_info.value.code == 2
    assert "would overwrite an input file" in capsys.readouterr().err
    assert (tmp_path / "input.fasta").read_bytes() == before


@pytest.mark.parametrize("flag", ["-o", "--cut-table", "--report"])
def test_every_output_flag_refuses_the_input_path(tmp_path, capsys, flag):
    path = fasta(tmp_path, PLASMID)
    with pytest.raises(SystemExit):
        main(["digest", "-i", path, "-e", "EcoRI", flag, path, "-f", "none"])
    assert "would overwrite an input file" in capsys.readouterr().err


def test_two_outputs_may_not_share_a_path(tmp_path, capsys):
    out = str(tmp_path / "both.tsv")
    with pytest.raises(SystemExit) as exit_info:
        main([
            "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
            "-o", out, "--cut-table", out, "-f", "none",
        ])
    assert exit_info.value.code == 2
    assert "would both write to" in capsys.readouterr().err


def test_distinct_output_paths_are_allowed(tmp_path):
    fragments = tmp_path / "fragments.tsv"
    cuts = tmp_path / "cuts.tsv"
    report = tmp_path / "report.html"
    code, _ = run([
        "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
        "-o", str(fragments), "--cut-table", str(cuts), "--report", str(report),
        "-f", "none",
    ])
    assert code == 0
    assert fragments.exists() and cuts.exists() and report.exists()


def test_stdin_needs_no_output_path_check(tmp_path, monkeypatch):
    """'-' is not a file, so it cannot be overwritten and must not be refused."""
    out = tmp_path / "fragments.tsv"

    class Stdin:
        """Standard input as the reader sees it: a text stream over bytes."""

        def __init__(self, text):
            self.buffer = io.BytesIO(text.encode())

    monkeypatch.setattr("sys.stdin", Stdin(f">s\n{PLASMID}\n"))
    code, _ = run(["digest", "-i", "-", "-e", "EcoRI", "-o", str(out), "-f", "none"])
    assert code == 0
    assert out.exists()


def test_cuts_and_max_cuts_are_alternatives(tmp_path, capsys):
    """Silently preferring one over the other hid which selection was run."""
    with pytest.raises(SystemExit) as exit_info:
        main([
            "select", "-i", fasta(tmp_path, PLASMID), "--cuts", "1", "--max-cuts", "3",
        ])
    assert exit_info.value.code == 2
    assert "not allowed with argument" in capsys.readouterr().err


def test_select_still_defaults_to_one_cut(tmp_path):
    _, text = run(["select", "-i", fasta(tmp_path, PLASMID), "-f", "tsv"])
    rows = text.splitlines()[1:]
    assert rows
    assert all(row.split("\t")[2] == "1" for row in rows)


def test_naming_an_enzyme_refuses_catalogue_filters(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["enzymes", "EcoRI", "--min-site", "8"])
    assert exit_info.value.code == 2
    assert "narrow the catalogue" in capsys.readouterr().err


@pytest.mark.parametrize(
    "flag,value",
    [("--available", None), ("--min-site", "8"), ("--max-site", "4"),
     ("--supplier", "N"), ("--end-type", "blunt"), ("--site", "GAATTC")],
)
def test_each_filter_is_refused_beside_a_name(flag, value):
    argv = ["enzymes", "EcoRI", flag] + ([value] if value else [])
    with pytest.raises(SystemExit):
        main(argv)


def test_compatible_honours_the_output_format():
    _, table = run(["enzymes", "--compatible", "--available", "--min-site", "8"])
    _, tsv = run(["enzymes", "--compatible", "--available", "--min-site", "8", "-f", "tsv"])
    _, text = run(["enzymes", "--compatible", "--available", "--min-site", "8", "-f", "json"])
    _, nothing = run([
        "enzymes", "--compatible", "--available", "--min-site", "8", "-f", "none",
    ])
    assert "ligate" in table
    assert tsv.splitlines()[0].split("\t") == ["end", "enzymes"]
    payload = json.loads(text)
    assert payload["shared_ends"]
    assert nothing == ""


def test_compatible_counts_only_the_groups_it_prints():
    """The header used to count ends it then omitted for having one member."""
    _, text = run(["enzymes", "--compatible", "--available"])
    header = text.splitlines()[0]
    shown = sum(1 for line in text.splitlines() if line.startswith("  "))
    claimed = int(header.split()[0])
    assert claimed == shown


def test_compatible_json_and_table_agree_on_the_counts():
    _, text = run(["enzymes", "--compatible", "--available", "-f", "json"])
    payload = json.loads(text)
    _, table = run(["enzymes", "--compatible", "--available"])
    assert len(payload["shared_ends"]) == int(table.splitlines()[0].split()[0])
    if payload["ends_unique_to_one_enzyme"]:
        assert "only one enzyme each" in table


def test_partial_sites_are_reported_separately(tmp_path):
    """A BcgI site cut on one side only is not reported as uncut."""
    path = fasta(tmp_path, "CGAAAAAAATGC" + "A" * 40)
    _, text = run(["digest", "-i", path, "-e", "BcgI", "-f", "json"])
    payload = json.loads(text)
    assert payload["partial_sites"] == {"BcgI": 1}
    assert payload["uncut_sites"] == {}
    _, table = run(["digest", "-i", path, "-e", "BcgI"])
    assert "one side only" in table


def test_output_refuses_a_hard_link_to_the_input(tmp_path, capsys):
    """Comparing resolved paths misses this: a hard link is a second real name.

    The two paths differ and resolve differently, yet name one inode, so a
    path-only check lets the output through and the input is destroyed.
    """
    source = tmp_path / "input.fasta"
    source.write_text(f">s\n{PLASMID}\n")
    link = tmp_path / "link.fasta"
    os.link(source, link)
    before = source.read_bytes()

    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", str(source), "-e", "EcoRI", "-o", str(link), "-f", "none"])
    assert exit_info.value.code == 2
    message = capsys.readouterr().err
    assert "would overwrite an input file" in message
    assert "the same file as" in message
    assert source.read_bytes() == before


def test_output_refuses_a_symlink_to_the_input(tmp_path, capsys):
    source = tmp_path / "input.fasta"
    source.write_text(f">s\n{PLASMID}\n")
    link = tmp_path / "link.fasta"
    link.symlink_to(source)
    before = source.read_bytes()

    with pytest.raises(SystemExit):
        main(["digest", "-i", str(source), "-e", "EcoRI", "-o", str(link), "-f", "none"])
    assert "would overwrite an input file" in capsys.readouterr().err
    assert source.read_bytes() == before


def test_two_outputs_hard_linked_to_each_other_are_refused(tmp_path, capsys):
    first = tmp_path / "a.tsv"
    first.write_text("")
    second = tmp_path / "b.tsv"
    os.link(first, second)
    with pytest.raises(SystemExit):
        main([
            "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
            "-o", str(first), "--cut-table", str(second), "-f", "none",
        ])
    assert "would both write to" in capsys.readouterr().err


def test_protect_refuses_an_enzyme_that_might_cut_the_guarded_sequence(tmp_path):
    """GANTTC has no confirmed EcoRI site but one that cannot be ruled out.

    Protection is a guarantee, so an enzyme that might cut the insert must not
    come back as safe just because the ambiguity hid the site.
    """
    target = fasta(tmp_path, "GAATTC", filename="target.fasta")
    guard = fasta(tmp_path, "GANTTC", name="guard", filename="guard.fasta")

    _, unprotected = run(["select", "-i", target, "-f", "tsv"])
    assert "EcoRI" in unprotected

    _, text = run(["select", "-i", target, "-f", "json", "-p", guard])
    payload = json.loads(text)
    assert [item["name"] for item in payload["enzymes"]] == []
    assert "EcoRI" in payload["excluded_protected_uncertain"]


def test_protect_exclusions_are_explained_when_nothing_is_left(tmp_path):
    """An empty result must not read as 'nothing cuts here'."""
    target = fasta(tmp_path, "GAATTC", filename="target.fasta")
    guard = fasta(tmp_path, "GANTTC", name="guard", filename="guard.fasta")
    _, text = run(["select", "-i", target, "-p", guard])
    assert "ambiguity code" in text
    assert "EcoRI" in text


def test_an_unambiguous_protected_sequence_still_approves_enzymes(tmp_path):
    """The stricter rule must not reject everything out of hand."""
    target = fasta(tmp_path, PLASMID, filename="target.fasta")
    guard = fasta(tmp_path, "TTTTTTTTTTTT", name="guard", filename="guard.fasta")
    _, text = run(["select", "-i", target, "-f", "tsv", "-p", guard])
    assert len(text.splitlines()) > 1


def test_empty_protected_record_is_refused(tmp_path, capsys):
    """Skipping it would leave the requested protection simply unapplied."""
    target = fasta(tmp_path, PLASMID, filename="target.fasta")
    guard = tmp_path / "guard.fasta"
    guard.write_text(">guard\n\n")
    with pytest.raises(SystemExit) as exit_info:
        main(["select", "-i", target, "-p", str(guard)])
    assert exit_info.value.code == 2
    assert "nothing can be protected" in capsys.readouterr().err


def test_all_reports_only_enzymes_with_something_to_say(tmp_path):
    """--all applies the catalogue; reporting 712 silent enzymes is not useful."""
    _, text = run(["digest", "-i", fasta(tmp_path, "GAATTC"), "--all", "-f", "json"])
    payload = json.loads(text)
    assert payload["enzymes_applied"] == len(all_enzymes())
    assert payload["enzymes_reported"] < payload["enzymes_applied"]
    for entry in payload["enzymes"]:
        assert entry["cuts"] or entry["sites"], entry["name"]


def test_all_keeps_an_enzyme_whose_site_could_not_be_cut(tmp_path):
    """A site present but uncuttable is a finding, not silence."""
    _, text = run([
        "digest", "-i", fasta(tmp_path, "AAAAAAAAGGTCTC"), "--all", "-f", "json",
    ])
    payload = json.loads(text)
    assert "BsaI" in payload["uncut_sites"]
    assert any(entry["name"] == "BsaI" for entry in payload["enzymes"])


def test_named_enzymes_apply_without_the_all_filter(tmp_path):
    """Explicit -e still reports an enzyme that does not cut."""
    _, text = run([
        "digest", "-i", fasta(tmp_path, "GAATTC"), "-e", "BamHI", "-f", "json",
    ])
    payload = json.loads(text)
    assert [entry["name"] for entry in payload["enzymes"]] == ["BamHI"]
    assert payload["enzymes"][0]["cuts"] == 0


@pytest.mark.parametrize("fmt", ["json", "tsv", "none"])
def test_a_named_enzyme_honours_the_output_format(fmt):
    _, text = run(["enzymes", "EcoRI", "-f", fmt])
    if fmt == "none":
        assert text == ""
    elif fmt == "tsv":
        assert text.splitlines()[0].startswith("enzyme\t")
        assert text.splitlines()[1].startswith("EcoRI\t")
    else:
        payload = json.loads(text)
        assert [item["name"] for item in payload["enzymes"]] == ["EcoRI"]


def test_named_enzyme_json_carries_the_detail_the_table_shows():
    _, text = run(["enzymes", "SmaI", "-f", "json"])
    entry = json.loads(text)["enzymes"][0]
    assert entry["notation"] == "CCC^GGG"
    assert entry["overhang"] == ""
    assert "XmaI" in entry["neoschizomers"]
    assert entry["methylation_overlap"]["dam"] == "never"


def test_enzyme_list_marks_a_direction_dependent_overhang():
    """Printing one of two complementary overhangs as the answer misleads."""
    _, text = run(["enzymes", "--site", "CACGAG"])
    assert "ACGA/TCGT" in text
    assert "forward site" in text


def test_enzyme_list_does_not_mark_a_palindromic_overhang():
    _, text = run(["enzymes", "--site", "GAATTC"])
    assert "AATT" in text
    assert "/TTAA" not in text


def test_enzyme_tsv_carries_the_orientation_flag():
    _, text = run(["enzymes", "--site", "CACGAG", "-f", "tsv"])
    header = text.splitlines()[0].split("\t")
    row = dict(zip(header, text.splitlines()[1].split("\t")))
    assert row["overhang_depends_on_orientation"] == "yes"


def test_select_reports_the_end_read_from_the_sequence(tmp_path, examples):
    """BspQI cuts outside its site, so only the sequence gives its overhang."""
    _, text = run([
        "select", "-i", str(examples / "puc19.fasta"), "--circular",
        "--cuts", "1", "-f", "json",
    ])
    payload = json.loads(text)
    entry = next(item for item in payload["enzymes"] if item["name"] == "BspQI")
    assert entry["ends_observed"] == [{"end_type": "5' overhang", "overhang": "GCT"}]


def test_digest_explains_a_cut_count_above_the_fragment_count(examples):
    _, text = run([
        "digest", "-i", str(examples / "puc19.fasta"), "-e", "AccII,Alw26I",
        "--circular", "--no-map", "--no-gel",
    ])
    assert "14 cuts fall on 13 boundaries" in text
    assert "AccII+Alw26I at 3" in text


@pytest.mark.parametrize("value", ["nan", "inf", "NaN", "Infinity"])
def test_resolution_rejects_values_that_are_not_finite(tmp_path, value, capsys):
    """nan compares false against everything; inf exceeds every difference.

    Either one silently calls 100 bp and 10,000 bp unresolved.
    """
    with pytest.raises(SystemExit) as exit_info:
        main([
            "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
            "--resolution", value,
        ])
    assert exit_info.value.code == 2
    assert "finite" in capsys.readouterr().err


def test_a_negative_resolution_is_refused(tmp_path, capsys):
    """argparse reads a leading '-' as a flag, so this fails earlier."""
    with pytest.raises(SystemExit) as exit_info:
        main([
            "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
            "--resolution", "-5",
        ])
    assert exit_info.value.code == 2
    assert capsys.readouterr().err


def test_resolution_still_accepts_ordinary_values(tmp_path):
    for value in ["0", "5", "12.5", "100"]:
        code, _ = run([
            "digest", "-i", fasta(tmp_path, PLASMID), "-e", "EcoRI",
            "--resolution", value, "-f", "none",
        ])
        assert code == 0


def test_truncated_gzip_reports_an_input_error_not_a_traceback(tmp_path, capsys):
    import gzip as _gzip

    whole = _gzip.compress(f">s\n{PLASMID}\n".encode())
    path = tmp_path / "cut.fasta.gz"
    path.write_bytes(whole[: len(whole) // 2])
    with pytest.raises(SystemExit) as exit_info:
        main(["digest", "-i", str(path), "-e", "EcoRI", "-f", "none"])
    assert exit_info.value.code == 2
    assert "truncated" in capsys.readouterr().err
