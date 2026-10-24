"""Self-contained HTML reporting.

The document embeds its own style rules, uses no JavaScript and loads nothing
from the network, so a report stays readable from a local file or an email
attachment years later. Output is byte-deterministic for a given input and set
of parameters: nothing in it depends on the clock or the filesystem beyond the
sequence name passed in.
"""

from __future__ import annotations

from html import escape
from typing import Sequence

from restrict_cut.digest import Cut, Digest, Fragment
from restrict_cut.enzymes import Enzyme, isoschizomers, neoschizomers
from restrict_cut.gel import DEFAULT_RESOLUTION, Band, summary as gel_summary
from restrict_cut.methylation import Overlap, System
from restrict_cut.text import overhang_display, plural, thousands

_STYLE = """
html { -webkit-text-size-adjust: 100%; }
body { font: 16px/1.55 Georgia, 'Times New Roman', serif; color: #1a1a1a;
       background: #fff; margin: 0 auto; padding: 3rem 1.25rem 4rem;
       max-width: 48rem; }
h1 { font-size: 1.5rem; font-weight: normal; margin: 0 0 .25rem; }
h2, h3 { font-size: .8rem; font-variant: small-caps; letter-spacing: .08em;
     font-weight: normal; margin: 2.5rem 0 .6rem; border-bottom: 1px solid #1a1a1a;
     padding-bottom: .2rem; }
h3 { margin-top: 1.8rem; }
p { margin: .5rem 0; }
.lede { color: #555; font-size: .9rem; margin-bottom: .5rem; }
.note { color: #555; font-size: .85rem; }
code, .seq, td.num, th.num, .mono { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
table { border-collapse: collapse; width: 100%; margin: .6rem 0 .2rem;
        font-size: .9rem; }
th { text-align: left; font-weight: normal; font-variant: small-caps;
     letter-spacing: .04em; border-bottom: 1px solid #1a1a1a; padding: .25rem .5rem .25rem 0; }
td { padding: .22rem .5rem .22rem 0; border-bottom: 1px solid #e6e6e6;
     vertical-align: baseline; }
td.num, th.num { text-align: right; padding-right: 1rem; font-size: .85rem; }
tr:last-child td { border-bottom: none; }
dl.metrics { display: grid; grid-template-columns: max-content 1fr; gap: .15rem 1.2rem;
             margin: .6rem 0; font-size: .9rem; }
dl.metrics dt { color: #555; }
dl.metrics dd { margin: 0; font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace; }
.bar { display: block; height: 2px; background: #1a1a1a; margin: .3rem 0 .1rem; }
.tag { font-variant: small-caps; letter-spacing: .04em; font-size: .8rem; color: #0f6f6a; }
.tag.warn { color: #8a5a00; }
pre.map { font-family: ui-monospace, 'SF Mono', Menlo, Consolas, monospace;
     font-size: .78rem; line-height: 1.35; overflow-x: auto; background: #fafafa;
     border: 1px solid #e6e6e6; padding: .7rem .8rem; margin: .6rem 0; }
footer { margin-top: 3rem; padding-top: .5rem; border-top: 1px solid #e6e6e6;
         color: #555; font-size: .8rem; }
@media print { body { padding: 0; max-width: none; } h2, h3 { page-break-after: avoid; } }
"""

#: Bar width as a percentage of the longest fragment.
_BAR_WIDTH = 100.0


def _definitions(rows: Sequence[tuple[str, str]]) -> list[str]:
    parts = ['<dl class="metrics">']
    for name, value in rows:
        parts.append(f"<dt>{escape(name)}</dt><dd>{escape(value)}</dd>")
    parts.append("</dl>")
    return parts


def _enzyme_table(result: Digest, sequence: str) -> list[str]:
    parts = [
        "<table>",
        '<tr><th>Enzyme</th><th>Site</th><th class="num">Sites</th>'
        '<th class="num">Cuts</th><th>End</th><th>Overhang</th></tr>',
    ]
    counts = result.cut_counts
    for enzyme in result.enzymes:
        cuts = result.cuts_for(enzyme.name)
        overhang = (
            overhang_display(cuts[0].overhang(sequence, result.circular), cuts[0].end_type)
            if cuts
            else "-"
        )
        parts.append(
            f"<tr><td>{escape(enzyme.name)}</td>"
            f"<td class=\"mono\">{escape(enzyme.notation)}</td>"
            f'<td class="num">{len(result.sites_for(enzyme.name))}</td>'
            f'<td class="num">{counts[enzyme.name]}</td>'
            f"<td>{escape(enzyme.end_type())}</td>"
            f"<td class=\"mono\">{escape(overhang)}</td></tr>"
        )
    parts.append("</table>")
    return parts


def _cut_table(cuts: Sequence[Cut], sequence: str, circular: bool) -> list[str]:
    parts = [
        "<table>",
        '<tr><th>Enzyme</th><th class="num">Top nick</th>'
        '<th class="num">Bottom nick</th><th class="num">Site</th>'
        "<th>Strand</th><th>End</th></tr>",
    ]
    for cut in cuts:
        overhang = overhang_display(cut.overhang(sequence, circular), cut.end_type)
        parts.append(
            f"<tr><td>{escape(cut.enzyme)}</td>"
            f'<td class="num">{thousands(cut.top)}</td>'
            f'<td class="num">{thousands(cut.bottom)}</td>'
            f'<td class="num">{thousands(cut.site_start)}</td>'
            f"<td>{escape(cut.strand)}</td>"
            f"<td class=\"mono\">{escape(overhang)}</td></tr>"
        )
    parts.append("</table>")
    return parts


def _fragment_table(fragments: Sequence[Fragment], longest: int) -> list[str]:
    parts = [
        "<table>",
        '<tr><th class="num">Length</th><th class="num">Start</th>'
        '<th class="num">End</th><th>5&prime; end</th><th>3&prime; end</th>'
        "<th></th></tr>",
    ]
    for fragment in fragments:
        width = max(1.0, _BAR_WIDTH * fragment.length / longest) if longest else 0.0
        bar = f'<span class="bar" style="width: {width:.1f}%"></span>'
        note = " (spans origin)" if fragment.spans_origin else ""
        if fragment.circular:
            note = " (uncut circle)"
        parts.append(
            f'<tr><td class="num">{thousands(fragment.length)}{escape(note)}</td>'
            f'<td class="num">{thousands(fragment.start)}</td>'
            f'<td class="num">{thousands(fragment.end)}</td>'
            f"<td>{escape(fragment.left)}</td>"
            f"<td>{escape(fragment.right)}</td>"
            f"<td>{bar}</td></tr>"
        )
    parts.append("</table>")
    return parts


def _gel_table(bands: Sequence[Band], resolution: float) -> list[str]:
    parts = [
        f'<p class="note">{escape(gel_summary(list(bands), resolution))}</p>',
        "<table>",
        '<tr><th class="num">Size</th><th class="num">Fragments</th>'
        "<th>Separation</th><th></th></tr>",
    ]
    for band in bands:
        if band.is_resolved:
            separation = "resolved"
            tag = '<span class="tag">resolved</span>'
        else:
            close = ", ".join(thousands(size) for size in band.unresolved_with)
            separation = f"within {resolution:g}% of {close}"
            tag = f'<span class="tag warn">{escape(separation)}</span>'
        indent = band.offset * 70.0
        bar = (
            f'<span class="bar" style="width: 18%; margin-left: {indent:.1f}%"></span>'
        )
        parts.append(
            f'<tr><td class="num">{thousands(band.size)}</td>'
            f'<td class="num">{band.count}</td>'
            f"<td>{tag}</td><td>{bar}</td></tr>"
        )
    parts.append("</table>")
    return parts


def _methylation_section(
    overlaps: Sequence[Overlap], systems: Sequence[System]
) -> list[str]:
    parts = ["<h2>Methylation</h2>"]
    names = ", ".join(f"{system.name} ({system.label})" for system in systems)
    parts.append(
        f'<p class="note">Checked against {escape(names)}. A site listed here '
        "carries a methylated base inside its recognition sequence. That is not "
        "the same as being blocked: BamHI's site always contains a Dam GATC and "
        "BamHI cuts it anyway, Dam methylation stops MboI outright, and DpnI "
        "needs the methylation before it will cut at all. Check the enzyme "
        "against its supplier's sensitivity table before changing a plan.</p>"
    )
    if not overlaps:
        parts.append(
            "<p>No recognition site of the enzymes applied carries a methylated "
            "base from the activities checked.</p>"
        )
        return parts
    parts.extend(
        [
            "<table>",
            '<tr><th>Enzyme</th><th class="num">Site</th><th>Activity</th>'
            '<th class="num">Modified base</th></tr>',
        ]
    )
    for overlap in overlaps:
        parts.append(
            f"<tr><td>{escape(overlap.enzyme)}</td>"
            f'<td class="num">{thousands(overlap.site_start)}</td>'
            f"<td>{escape(overlap.system)}</td>"
            f'<td class="num">{thousands(overlap.position)}</td></tr>'
        )
    parts.append("</table>")
    return parts


def _relationships(enzymes: Sequence[Enzyme]) -> list[str]:
    rows: list[tuple[str, str]] = []
    for enzyme in enzymes:
        same = isoschizomers(enzyme)
        other = neoschizomers(enzyme)
        if same:
            rows.append((f"{enzyme.name} isoschizomers", ", ".join(same)))
        if other:
            rows.append((f"{enzyme.name} same site, different cut", ", ".join(other)))
    if not rows:
        return []
    parts = ["<h3>Related enzymes</h3>"]
    parts.extend(_definitions(rows))
    parts.append(
        '<p class="note">Isoschizomers recognise the same site and cut it in the '
        "same place, so they are interchangeable. An enzyme listed as cutting the "
        "same site differently leaves a different end, and is not.</p>"
    )
    return parts


def render_html(
    result: Digest,
    sequence: str,
    parameters: dict[str, str],
    version: str,
    bands: Sequence[Band] = (),
    resolution: float = DEFAULT_RESOLUTION,
    overlaps: Sequence[Overlap] = (),
    systems: Sequence[System] = (),
    map_lines: Sequence[str] = (),
) -> str:
    """Render a complete, self-contained HTML document for one digest."""
    longest = max((fragment.length for fragment in result.fragments), default=0)
    enzyme_names = ", ".join(enzyme.name for enzyme in result.enzymes)

    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>Restriction map of {escape(result.name)}</title>",
        f"<style>{_STYLE}</style>",
        "</head>",
        "<body>",
        f"<h1>Restriction map of {escape(result.name)}</h1>",
        f'<p class="lede">{escape(enzyme_names)} &middot; '
        f"{thousands(result.length)} bp &middot; {escape(result.topology)} &middot; "
        f"{escape(plural(len(result.cuts), 'cut'))} &middot; "
        f"{escape(plural(len(result.fragments), 'fragment'))}</p>",
        "<h2>Parameters</h2>",
    ]
    parts.extend(_definitions(list(parameters.items())))

    parts.append("<h2>Enzymes</h2>")
    parts.extend(_enzyme_table(result, sequence))

    if result.uncut_sites:
        listed = ", ".join(
            f"{name} ({plural(count, 'site')})"
            for name, count in sorted(result.uncut_sites.items())
        )
        parts.append(
            f'<p class="note">Sites present but not cut, because the enzyme cuts '
            f"outside its recognition sequence and every nick would fall beyond "
            f"the end of this linear molecule: {escape(listed)}.</p>"
        )
    if result.partial_sites:
        listed = ", ".join(
            f"{name} ({plural(count, 'site')})"
            for name, count in sorted(result.partial_sites.items())
        )
        parts.append(
            f'<p class="note">Sites cut on one side only: the enzyme cuts on both '
            f"sides of its recognition sequence, but here one of the two breaks "
            f"falls beyond the end of the molecule, so the site is cut without "
            f"being excised: {escape(listed)}.</p>"
        )
    if result.uncertain:
        listed = ", ".join(
            f"{name} ({plural(count, 'position')})"
            for name, count in sorted(result.uncertain.items())
        )
        parts.append(
            f'<p class="note">Positions that could not be called either way, '
            f"because the sequence carries an ambiguity code there: "
            f"{escape(listed)}. They are not counted as sites.</p>"
        )

    if map_lines:
        parts.append("<h2>Map</h2>")
        parts.append(f'<pre class="map">{escape(chr(10).join(map_lines))}</pre>')

    if result.cuts:
        parts.append("<h2>Cuts</h2>")
        parts.extend(_cut_table(result.cuts, sequence, result.circular))
        parts.append(
            '<p class="note">A nick written <em>n</em> falls between base '
            "<em>n</em> and base <em>n</em>&nbsp;+&nbsp;1 on the forward strand. "
            "The overhang is the forward-strand text between the two nicks, which "
            "is the form in which two ends can be compared: ends anneal when both "
            "the text and the end type agree.</p>"
        )

    parts.append("<h2>Fragments</h2>")
    parts.extend(_fragment_table(result.fragments, longest))
    parts.append(
        f'<p class="note">Lengths are measured between top-strand nicks and sum '
        f"to {thousands(result.total_length)} bp, the length of the molecule. On a "
        f"circular molecule <em>n</em> cuts give <em>n</em> fragments; on a linear "
        f"one they give <em>n</em>&nbsp;+&nbsp;1.</p>"
    )

    if bands:
        parts.append("<h2>By size</h2>")
        parts.extend(_gel_table(bands, resolution))
        parts.append(
            '<p class="note">Sizes are placed on a logarithmic axis. This is an '
            "ordering by size with the bands that will not separate marked, not a "
            "simulation of electrophoresis: no agarose percentage, buffer or run "
            "time is modelled.</p>"
        )

    if systems:
        parts.extend(_methylation_section(overlaps, systems))

    parts.extend(_relationships(result.enzymes))

    parts.extend(
        [
            "<footer>",
            "<p>Coordinates are 1-based and inclusive on the forward strand. A "
            "fragment whose end is smaller than its start spans the origin of a "
            "circular molecule.</p>",
            f"<p>restrict-cut {escape(version)}. Enzyme specificities from REBASE "
            "(Roberts <em>et al.</em>, Nucleic Acids Research).</p>",
            "</footer>",
            "</body>",
            "</html>",
            "",
        ]
    )
    return "\n".join(parts)
