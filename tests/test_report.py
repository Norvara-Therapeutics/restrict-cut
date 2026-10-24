"""The self-contained HTML report."""

from restrict_cut import __version__, enzymes
from restrict_cut.digest import digest
from restrict_cut.gel import lanes
from restrict_cut.methylation import resolve_systems, site_overlaps
from restrict_cut.report import render_html

SEQUENCE = "GAATTCAAAAAAAAAAGGATCCAAAAAAAAAACTGCAGAAAAAAAAAA"


def build(circular=False, with_extras=True):
    result = digest(
        SEQUENCE,
        enzymes.resolve(["EcoRI", "BamHI", "PstI"]),
        circular=circular,
        name="test",
    )
    systems = resolve_systems(None) if with_extras else []
    overlaps = []
    for enzyme in result.enzymes:
        overlaps.extend(
            site_overlaps(
                SEQUENCE, enzyme, result.sites_for(enzyme.name), systems, circular
            )
        )
    bands = lanes([fragment.length for fragment in result.fragments])
    return render_html(
        result,
        SEQUENCE,
        parameters={"Input": "test.fasta"},
        version=__version__,
        bands=bands if with_extras else (),
        overlaps=overlaps,
        systems=systems,
        map_lines=["map line"] if with_extras else (),
    )


def test_report_is_a_complete_document():
    html = build()
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")


def test_report_loads_nothing_from_the_network():
    html = build()
    assert "http://" not in html.replace("http://www.w3.org", "")
    assert "<script" not in html
    assert "src=" not in html


def test_report_embeds_its_own_style():
    assert "<style>" in build()


def test_report_is_deterministic():
    assert build() == build()


def test_report_names_the_enzymes_and_fragments():
    html = build()
    for name in ["EcoRI", "BamHI", "PstI"]:
        assert name in html
    assert "Fragments" in html
    assert "48" in html


def test_report_states_the_fragment_sum():
    assert f"{len(SEQUENCE)} bp, the length of the molecule" in build()


def test_report_explains_the_methylation_caveat():
    html = build()
    assert "not the same as being blocked" in html
    assert "DpnI" in html


def test_report_disclaims_simulating_electrophoresis():
    assert "not a" in build() and "simulation of electrophoresis" in build()


def test_report_records_the_version_and_source():
    html = build()
    assert __version__ in html
    assert "REBASE" in html


def test_report_works_without_the_optional_sections():
    html = build(with_extras=False)
    assert "<!DOCTYPE html>" in html
    assert "Fragments" in html


def test_circular_report_mentions_the_origin():
    assert "origin" in build(circular=True)


def test_report_escapes_its_input():
    result = digest("GAATTC", [enzymes.get("EcoRI")], name="<script>x</script>")
    html = render_html(result, "GAATTC", {"Input": "<b>"}, __version__)
    assert "<script>x</script>" not in html
    assert "&lt;script&gt;" in html
