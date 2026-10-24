"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import IO, Sequence

from restrict_cut import __version__
from restrict_cut import enzymes as catalogue
from restrict_cut.digest import Digest, digest
from restrict_cut.enzymes import END_TYPES, SUPPLIERS, Enzyme
from restrict_cut.gel import (
    DEFAULT_RESOLUTION,
    Band,
    lanes,
    render as render_gel,
    summary as gel_summary,
)
from restrict_cut.iupac import CODES, reverse_complement
from restrict_cut.mapview import render as render_map
from restrict_cut.methylation import (
    DEFAULT_SYSTEMS,
    SYSTEMS,
    exposure,
    resolve_systems,
    site_overlaps,
)
from restrict_cut.reader import SequenceError, read_one, read_records
from restrict_cut.selection import choose
from restrict_cut.text import base_pairs, overhang_display, plural, thousands

STDOUT_FORMATS = ("table", "tsv", "json", "none")

_FRAGMENT_COLUMNS = ("sequence", "length", "start", "end", "five_prime", "three_prime")
_CUT_COLUMNS = (
    "sequence",
    "enzyme",
    "top_nick",
    "bottom_nick",
    "site_start",
    "strand",
    "end_type",
    "overhang",
)
_SELECT_COLUMNS = (
    "enzyme",
    "site",
    "cuts",
    "positions",
    "end_type",
    "overhang",
    "site_length",
    "suppliers",
)
_ENZYME_COLUMNS = (
    "enzyme",
    "site",
    "notation",
    "site_length",
    "specific_bases",
    "end_type",
    "overhang",
    "overhang_depends_on_orientation",
    "cuts_outside_site",
    "expected_frequency",
    "suppliers",
)


def _non_negative(text: str) -> int:
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from None
    if value < 0:
        raise argparse.ArgumentTypeError(f"must not be negative: {value}")
    return value


def _positive(text: str) -> int:
    value = _non_negative(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"must be 1 or greater: {value}")
    return value


def _percentage(text: str) -> float:
    """Read a non-negative, finite percentage.

    ``float`` accepts ``nan`` and ``inf``, and a negative test lets both
    through: every comparison against ``nan`` is false, and ``inf`` exceeds any
    difference, so either one silently marks every pair of bands as
    unresolved - including 100 bp beside 10,000 bp.
    """
    try:
        value = float(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {text!r}") from None
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError(
            f"must be a finite number, not {text!r}"
        )
    if value < 0:
        raise argparse.ArgumentTypeError(f"must not be negative: {value}")
    return value


def _enzyme_list(text: str) -> list[str]:
    names = [part.strip() for part in text.replace(" ", ",").split(",")]
    chosen = [name for name in names if name]
    if not chosen:
        raise argparse.ArgumentTypeError("no enzyme names given")
    return chosen


def _site(text: str) -> str:
    pattern = text.strip().upper()
    if not pattern:
        raise argparse.ArgumentTypeError("empty recognition site")
    unknown = sorted(set(pattern) - set(CODES))
    if unknown:
        shown = ", ".join(repr(code) for code in unknown)
        raise argparse.ArgumentTypeError(f"not IUPAC nucleotide codes: {shown}")
    return pattern


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        prog="restrict-cut",
        description=(
            "Map restriction sites, simulate digests and choose cloning enzymes, "
            f"from the {len(catalogue.all_enzymes())} REBASE specificities whose "
            "cut position is known."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"restrict-cut {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    _add_digest_parser(subparsers)
    _add_select_parser(subparsers)
    _add_enzymes_parser(subparsers)
    return parser


def _add_catalogue_filters(parser: argparse.ArgumentParser) -> None:
    """Flags shared by the commands that search the whole catalogue."""
    parser.add_argument(
        "--available",
        action="store_true",
        help="consider only enzymes a company sells",
    )
    parser.add_argument(
        "--supplier",
        metavar="CODE",
        help=(
            "consider only enzymes sold by this REBASE supplier code "
            f"({', '.join(sorted(SUPPLIERS))})"
        ),
    )
    parser.add_argument(
        "--min-site",
        type=_non_negative,
        default=0,
        metavar="N",
        help="consider only sites at least N bases long (default: 0, no bound)",
    )
    parser.add_argument(
        "--max-site",
        type=_non_negative,
        default=0,
        metavar="N",
        help="consider only sites at most N bases long (default: 0, no bound)",
    )
    parser.add_argument(
        "--end-type",
        choices=END_TYPES,
        help="consider only enzymes leaving this kind of end",
    )


def _add_digest_parser(subparsers: argparse._SubParsersAction) -> None:
    digest_parser = subparsers.add_parser(
        "digest",
        help="map and cut one sequence with named enzymes",
        description=(
            "Cut one sequence with the enzymes named and report the sites, the "
            "cuts, the fragments they produce and how those fragments would "
            "order by size. Give --circular for a plasmid: on a circular "
            "molecule a site may span the origin and n cuts give n fragments, "
            "not n + 1."
        ),
    )
    digest_parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="PATH",
        help="FASTA file holding the sequence, or '-' for standard input",
    )
    digest_parser.add_argument(
        "-e",
        "--enzyme",
        type=_enzyme_list,
        action="append",
        metavar="NAMES",
        help="enzyme to apply, comma-separated or repeated; omit for --all",
    )
    digest_parser.add_argument(
        "--all",
        action="store_true",
        help=(
            "apply every enzyme in the catalogue, subject to the filters below, "
            "and report those that cut or have a site"
        ),
    )
    digest_parser.add_argument(
        "-c",
        "--circular",
        action="store_true",
        help="treat the sequence as a closed circle, as a plasmid is",
    )
    digest_parser.add_argument(
        "--methylation",
        metavar="NAMES",
        help=(
            "methyltransferase activities to check sites against: "
            f"{', '.join(sorted(SYSTEMS))}, or 'none' "
            f"(default: {','.join(DEFAULT_SYSTEMS)})"
        ),
    )
    digest_parser.add_argument(
        "--resolution",
        type=_percentage,
        default=DEFAULT_RESOLUTION,
        metavar="PERCENT",
        help=(
            "call two fragments unresolved when their sizes differ by less than "
            f"this percentage (default: {DEFAULT_RESOLUTION:g})"
        ),
    )
    digest_parser.add_argument(
        "--ladder",
        metavar="NAME|SIZES",
        help=(
            "comma-separated fragment sizes to show as a ladder beside the "
            "sample, as --ladder 1000,500,250"
        ),
    )
    digest_parser.add_argument(
        "--map-width",
        type=_positive,
        default=64,
        metavar="N",
        help="width of the text map in characters (default: 64)",
    )
    digest_parser.add_argument(
        "--no-map", action="store_true", help="omit the map from the summary"
    )
    digest_parser.add_argument(
        "--no-gel", action="store_true", help="omit the size ordering from the summary"
    )
    digest_parser.add_argument(
        "--report", metavar="PATH", help="write a self-contained HTML report to PATH"
    )
    digest_parser.add_argument(
        "-o",
        "--output",
        metavar="PATH",
        help="write one TSV row per fragment to PATH",
    )
    digest_parser.add_argument(
        "--cut-table",
        metavar="PATH",
        help="write one TSV row per cut to PATH, with nick positions and overhang",
    )
    digest_parser.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        dest="stdout_format",
        help="what to print on standard output (default: table)",
    )
    _add_catalogue_filters(digest_parser)


def _add_select_parser(subparsers: argparse._SubParsersAction) -> None:
    select_parser = subparsers.add_parser(
        "select",
        help="find the enzymes that cut a sequence a given number of times",
        description=(
            "Search the catalogue for enzymes that cut a target a set number of "
            "times, optionally leaving other sequences untouched. Asking for one "
            "cut finds the enzymes that linearise a vector; asking for none finds "
            "the sites that can safely be introduced into it."
        ),
    )
    select_parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="PATH",
        help="FASTA file holding the target sequence, or '-' for standard input",
    )
    select_parser.add_argument(
        "-c",
        "--circular",
        action="store_true",
        help="treat the target as a closed circle",
    )
    counts = select_parser.add_mutually_exclusive_group()
    counts.add_argument(
        "-n",
        "--cuts",
        type=_non_negative,
        default=None,
        metavar="N",
        help="required number of cuts in the target (default: 1)",
    )
    counts.add_argument(
        "--max-cuts",
        type=_non_negative,
        metavar="N",
        help="accept any number of cuts up to N, instead of an exact count",
    )
    select_parser.add_argument(
        "-p",
        "--protect",
        action="append",
        metavar="PATH",
        help=(
            "FASTA file whose sequence must not be cut at all; repeat for "
            "several. An insert belongs here: an enzyme that cuts it is no use "
            "for moving it, however cleanly it cuts the vector"
        ),
    )
    select_parser.add_argument(
        "--protect-circular",
        action="store_true",
        help="treat the protected sequences as circular too",
    )
    select_parser.add_argument(
        "--limit",
        type=_non_negative,
        default=0,
        metavar="N",
        help="show only the first N enzymes (default: 0, show all)",
    )
    select_parser.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        dest="stdout_format",
        help="what to print on standard output (default: table)",
    )
    _add_catalogue_filters(select_parser)


def _add_enzymes_parser(subparsers: argparse._SubParsersAction) -> None:
    enzymes_parser = subparsers.add_parser(
        "enzymes",
        help="browse the enzyme catalogue",
        description=(
            "List the catalogue, or describe one enzyme in full: how it cuts, "
            "what end it leaves, which enzymes it can be swapped for, and whether "
            "its site can carry a methylated base."
        ),
    )
    enzymes_parser.add_argument(
        "name",
        nargs="?",
        help="describe this enzyme in detail instead of listing the catalogue",
    )
    enzymes_parser.add_argument(
        "--site",
        type=_site,
        metavar="SEQUENCE",
        help="list only enzymes recognising exactly this sequence",
    )
    enzymes_parser.add_argument(
        "--compatible",
        action="store_true",
        help="group enzymes by the end they leave, showing which can be ligated together",
    )
    enzymes_parser.add_argument(
        "-f",
        "--format",
        choices=STDOUT_FORMATS,
        default="table",
        dest="stdout_format",
        help="what to print on standard output (default: table)",
    )
    _add_catalogue_filters(enzymes_parser)


#: Output flags of the digest command, in the order they are written.
_OUTPUT_FLAGS = (("--output", "output"), ("--cut-table", "cut_table"), ("--report", "report"))


def _identity(path: str) -> tuple[object, ...]:
    """Identify the file a path names, not the path itself.

    Comparing resolved paths is not enough. ``Path.resolve`` follows symbolic
    links, but a *hard* link is a second genuine name for one inode: the two
    paths differ, resolve differently, and still refer to the same bytes on
    disk. An output aimed at a hard link of the input would pass a path
    comparison and then overwrite the input anyway.

    An existing file is therefore identified by its device and inode, which is
    what the filesystem itself considers sameness. A file that does not exist
    yet has no inode, so it falls back to the resolved path; that is sound,
    because a path with nothing at it cannot already be an input.
    """
    expanded = Path(path).expanduser()
    try:
        status = expanded.stat()
    except OSError:
        return ("path", expanded.resolve())
    return ("inode", status.st_dev, status.st_ino)


def _check_output_paths(args: argparse.Namespace) -> None:
    """Refuse to overwrite the input, or to send two outputs to one file.

    Each output is written independently and the sequence has already been read
    by then, so without this check ``--output`` aimed at the input file
    destroys the source and still exits successfully. Two outputs sharing a
    file would quietly leave only whichever was written last.
    """
    inputs: dict[tuple[object, ...], str] = {}
    if args.input != "-":
        inputs[_identity(args.input)] = args.input

    claimed: dict[tuple[object, ...], str] = {}
    for flag, attribute in _OUTPUT_FLAGS:
        value = getattr(args, attribute, None)
        if value is None:
            continue
        identity = _identity(value)
        if identity in inputs:
            named = inputs[identity]
            also = "" if Path(named) == Path(value) else f" (the same file as {named})"
            raise ValueError(f"{flag} would overwrite an input file: {value}{also}")
        if identity in claimed:
            raise ValueError(
                f"{claimed[identity]} and {flag} would both write to {value}"
            )
        claimed[identity] = flag


def _filtered_catalogue(args: argparse.Namespace) -> list[Enzyme]:
    """Apply the shared catalogue filters."""
    return catalogue.select(
        available_only=args.available,
        min_length=args.min_site,
        max_length=args.max_site,
        end_type=args.end_type,
        supplier=args.supplier,
    )


def _chosen_enzymes(args: argparse.Namespace) -> list[Enzyme]:
    """Resolve the enzymes a digest should apply."""
    names: list[str] = []
    for group in args.enzyme or []:
        names.extend(group)
    if names and args.all:
        raise ValueError("give either --enzyme or --all, not both")
    if names:
        return catalogue.resolve(names)
    if not args.all:
        raise ValueError(
            "name the enzymes with --enzyme, or pass --all to apply every enzyme "
            "that cuts"
        )
    return _filtered_catalogue(args)


def _resolve_ladder(text: str | None) -> tuple[list[int], str]:
    """Read the ladder sizes given on the command line.

    Only explicit sizes are accepted. No supplier's marker is built in: the
    sizes of a commercial ladder are a product specification, and a ladder
    derived from a reference sequence can be produced with a digest of its own,
    which keeps it checkable.
    """
    if not text:
        return [], ""
    sizes = []
    for part in text.replace(" ", ",").split(","):
        if not part:
            continue
        try:
            size = int(part)
        except ValueError:
            raise ValueError(f"not a fragment size: {part!r}") from None
        if size < 1:
            raise ValueError(f"fragment size must be 1 or greater: {size}")
        sizes.append(size)
    if not sizes:
        raise ValueError(f"no ladder sizes in {text!r}")
    return sizes, "given sizes"


def _methylation_choice(text: str | None) -> list[str] | None:
    if text is None:
        return None
    if text.strip().lower() in {"none", ""}:
        return []
    return [part for part in text.replace(" ", ",").split(",") if part]


def format_digest_table(
    result: Digest,
    sequence: str,
    bands: list[Band],
    ladder_bands: list[Band],
    ladder_name: str,
    overlaps: list,
    systems: list,
    args: argparse.Namespace,
) -> str:
    """Render a readable summary of one digest."""
    lines: list[str] = [
        f"> {result.name}  {thousands(result.length)} bp  {result.topology}"
    ]
    counts = result.cut_counts
    for enzyme in result.enzymes:
        cuts = result.cuts_for(enzyme.name)
        site_count = len(result.sites_for(enzyme.name))
        if not cuts:
            detail = "does not cut"
            if site_count:
                detail = f"{plural(site_count, 'site')} but no cut within the molecule"
            lines.append(f"  {enzyme.name:<12} {enzyme.notation:<24} {detail}")
            continue
        overhang = overhang_display(
            cuts[0].overhang(sequence, result.circular), cuts[0].end_type
        )
        lines.append(
            f"  {enzyme.name:<12} {enzyme.notation:<24} "
            f"{plural(len(cuts), 'cut'):<10} {overhang}"
        )
    if len(result.enzymes) > 12:
        cutting = sum(1 for count in counts.values() if count)
        lines = lines[:1]
        lines.append(
            f"  {plural(len(result.enzymes), 'enzyme')} applied, "
            f"{cutting} of them cutting"
        )

    if result.uncut_sites:
        listed = ", ".join(
            f"{name} ({plural(count, 'site')})"
            for name, count in sorted(result.uncut_sites.items())
        )
        lines.append(f"  note       sites present but cut falls off the end: {listed}")
    if result.partial_sites:
        listed = ", ".join(
            f"{name} ({plural(count, 'site')})"
            for name, count in sorted(result.partial_sites.items())
        )
        lines.append(
            f"  note       sites cut on one side only, too near the end to "
            f"excise: {listed}"
        )
    if result.uncertain:
        listed = ", ".join(
            f"{name} ({plural(count, 'position')})"
            for name, count in sorted(result.uncertain.items())
        )
        lines.append(f"  note       uncalled, sequence is ambiguous there: {listed}")

    if result.is_uncut:
        lines.append("  no enzyme cuts this sequence")
        return "\n".join(lines) + "\n"

    if not args.no_map:
        lines.append("")
        lines.extend(
            ("  " + line).rstrip()
            for line in render_map(
                [(cut.top, cut.enzyme) for cut in result.cuts],
                result.length,
                width=args.map_width,
                circular=result.circular,
            )
        )

    lines.append("")
    lines.append(
        f"  fragments  {plural(len(result.fragments), 'fragment')}, "
        f"{thousands(result.total_length)} bp total"
    )
    shared = result.shared_boundaries
    if shared:
        listed = ", ".join(
            f"{'+'.join(names)} at {thousands(position)}"
            for position, names in list(shared.items())[:4]
        )
        if len(shared) > 4:
            listed += f", and {len(shared) - 4} more"
        lines.append(
            f"             {plural(len(result.cuts), 'cut')} fall on "
            f"{result.boundaries} boundaries: {listed} break the same bond"
        )
    for fragment in result.fragments[:20]:
        note = " spans origin" if fragment.spans_origin else ""
        if fragment.circular:
            note = " uncut circle"
        lines.append(
            f"    {base_pairs(fragment.length):>12}  "
            f"{thousands(fragment.start)}..{thousands(fragment.end)}  "
            f"{fragment.left} / {fragment.right}{note}"
        )
    if len(result.fragments) > 20:
        lines.append(f"    ... and {len(result.fragments) - 20} more")

    if bands and not args.no_gel:
        lines.append("")
        lines.append(f"  by size    {gel_summary(bands, args.resolution)}")
        lines.extend(
            ("  " + line).rstrip()
            for line in render_gel(
                bands, ladder=ladder_bands or None, ladder_name=ladder_name
            )
        )

    if systems:
        lines.append("")
        checked = ", ".join(system.name for system in systems)
        if overlaps:
            affected = len({(item.enzyme, item.site_start) for item in overlaps})
            lines.append(
                f"  methylation  {plural(affected, 'site')} carrying "
                f"{plural(len(overlaps), 'methylated base')} ({checked})"
            )
            for overlap in overlaps[:10]:
                lines.append(
                    f"    {overlap.enzyme} site at {thousands(overlap.site_start)}: "
                    f"{overlap.system} modifies base {thousands(overlap.position)}"
                )
            if len(overlaps) > 10:
                lines.append(f"    ... and {len(overlaps) - 10} more")
            lines.append(
                "    overlap is not the same as blocked; check the enzyme against "
                "its supplier's table"
            )
        else:
            lines.append(f"  methylation  no site carries a methylated base ({checked})")

    return "\n".join(lines) + "\n"


def format_digest_tsv(result: Digest) -> str:
    """One TSV row per fragment."""
    rows = [list(_FRAGMENT_COLUMNS)]
    for fragment in result.fragments:
        rows.append(
            [
                result.name,
                str(fragment.length),
                str(fragment.start),
                str(fragment.end),
                fragment.left,
                fragment.right,
            ]
        )
    return "".join("\t".join(row) + "\n" for row in rows)


def format_cuts_tsv(result: Digest, sequence: str) -> str:
    """One TSV row per cut."""
    rows = [list(_CUT_COLUMNS)]
    for cut in result.cuts:
        rows.append(
            [
                result.name,
                cut.enzyme,
                str(cut.top),
                str(cut.bottom),
                str(cut.site_start),
                cut.strand,
                cut.end_type,
                cut.overhang(sequence, result.circular),
            ]
        )
    return "".join("\t".join(row) + "\n" for row in rows)


def format_digest_json(
    result: Digest,
    sequence: str,
    bands: list[Band],
    overlaps: list,
    systems: list,
    parameters: dict[str, str],
) -> str:
    """Every cut, fragment and band as JSON."""
    payload = {
        "tool": "restrict-cut",
        "version": __version__,
        "rebase_release": catalogue.REBASE_RELEASE,
        "parameters": parameters,
        "sequence": {
            "name": result.name,
            "length": result.length,
            "topology": result.topology,
        },
        "enzymes_applied": result.considered or len(result.enzymes),
        "enzymes_reported": len(result.enzymes),
        "enzymes": [
            {
                "name": enzyme.name,
                "site": enzyme.pattern,
                "notation": enzyme.notation,
                "end_type": enzyme.end_type(),
                "cuts_outside_site": enzyme.cuts_outside_site,
                "sites": len(result.sites_for(enzyme.name)),
                "cuts": result.cut_counts[enzyme.name],
                "suppliers": list(enzyme.supplier_names),
            }
            for enzyme in result.enzymes
        ],
        "sites": [
            {
                "enzyme": site.enzyme,
                "start": site.start,
                "end": site.end,
                "strand": site.strand,
                "spans_origin": site.spans_origin,
            }
            for site in result.sites
        ],
        "cuts": [
            {
                "enzyme": cut.enzyme,
                "top_nick": cut.top,
                "bottom_nick": cut.bottom,
                "site_start": cut.site_start,
                "strand": cut.strand,
                "end_type": cut.end_type,
                "overhang": cut.overhang(sequence, result.circular),
            }
            for cut in result.cuts
        ],
        "fragments": [
            {
                "length": fragment.length,
                "start": fragment.start,
                "end": fragment.end,
                "five_prime": fragment.left,
                "three_prime": fragment.right,
                "spans_origin": fragment.spans_origin,
                "circular": fragment.circular,
            }
            for fragment in result.fragments
        ],
        "total_fragment_length": result.total_length,
        "bands": [
            {
                "size": band.size,
                "fragments": band.count,
                "resolved": band.is_resolved,
                "unresolved_with": list(band.unresolved_with),
            }
            for band in bands
        ],
        "uncertain_positions": result.uncertain,
        "uncut_sites": result.uncut_sites,
        "partial_sites": result.partial_sites,
        "methylation": {
            "systems": [
                {"name": system.name, "target": system.label, "base": system.base}
                for system in systems
            ],
            "overlaps": [
                {
                    "enzyme": overlap.enzyme,
                    "site_start": overlap.site_start,
                    "system": overlap.system,
                    "modified_base": overlap.position,
                }
                for overlap in overlaps
            ],
        },
    }
    return json.dumps(payload, indent=2) + "\n"


def _digest_parameters(args: argparse.Namespace, result: Digest) -> dict[str, str]:
    return {
        "Input": args.input,
        "Sequence": f"{thousands(result.length)} bp, {result.topology}",
        "Enzymes": f"{plural(len(result.enzymes), 'enzyme')} applied",
        "Resolution": f"{args.resolution:g}%",
        "Methylation": args.methylation or ",".join(DEFAULT_SYSTEMS),
        "REBASE release": catalogue.REBASE_RELEASE,
    }


def run_digest(args: argparse.Namespace, stdout: IO[str]) -> int:
    """Execute the digest subcommand. Returns a process exit code."""
    _check_output_paths(args)
    enzymes = _chosen_enzymes(args)
    record = read_one(args.input)
    if not record.sequence:
        raise SequenceError(f"{record.identifier} holds no sequence")

    result = digest(
        record.sequence, enzymes, circular=args.circular, name=record.identifier
    )
    if args.all:
        # Naming every enzyme in the catalogue and then reporting the 712 with
        # nothing to say about this sequence is noise, not information. An
        # enzyme is kept if it cut, or if it has a site it could not cut, or if
        # the sequence is ambiguous where it might have cut.
        result.considered = len(enzymes)
        speaks = {cut.enzyme for cut in result.cuts}
        speaks |= {site.enzyme for site in result.sites}
        speaks |= set(result.uncertain)
        result.enzymes = [item for item in result.enzymes if item.name in speaks]

    systems = resolve_systems(_methylation_choice(args.methylation))
    overlaps = []
    for enzyme in result.enzymes:
        overlaps.extend(
            site_overlaps(
                record.sequence,
                enzyme,
                result.sites_for(enzyme.name),
                systems,
                result.circular,
            )
        )
    overlaps.sort(key=lambda item: (item.site_start, item.enzyme, item.system))

    sizes = [fragment.length for fragment in result.fragments]
    labels: dict[int, tuple[str, ...]] = {}
    for fragment in result.fragments:
        labels.setdefault(fragment.length, ())
    bands = lanes(sizes, labels, args.resolution)

    ladder_sizes, ladder_name = _resolve_ladder(args.ladder)
    ladder_bands = lanes(ladder_sizes, None, args.resolution) if ladder_sizes else []

    map_lines = (
        []
        if args.no_map
        else render_map(
            [(cut.top, cut.enzyme) for cut in result.cuts],
            result.length,
            width=args.map_width,
            circular=result.circular,
        )
    )
    parameters = _digest_parameters(args, result)

    if args.output:
        Path(args.output).write_text(format_digest_tsv(result), encoding="utf-8")
    if args.cut_table:
        Path(args.cut_table).write_text(
            format_cuts_tsv(result, record.sequence), encoding="utf-8"
        )
    if args.report:
        from restrict_cut.report import render_html

        Path(args.report).write_text(
            render_html(
                result,
                record.sequence,
                parameters=parameters,
                version=__version__,
                bands=bands,
                resolution=args.resolution,
                overlaps=overlaps,
                systems=systems,
                map_lines=map_lines,
            ),
            encoding="utf-8",
        )

    if args.stdout_format == "table":
        stdout.write(
            format_digest_table(
                result,
                record.sequence,
                bands,
                ladder_bands,
                ladder_name,
                overlaps,
                systems,
                args,
            )
        )
    elif args.stdout_format == "tsv":
        stdout.write(format_digest_tsv(result))
    elif args.stdout_format == "json":
        stdout.write(
            format_digest_json(
                result, record.sequence, bands, overlaps, systems, parameters
            )
        )
    return 0


def _excluded_note(selection) -> list[str]:
    """Say which enzymes a protected sequence's ambiguity ruled out."""
    if not selection.excluded_uncertain:
        return []
    shown = ", ".join(selection.excluded_uncertain[:8])
    if len(selection.excluded_uncertain) > 8:
        shown += f", and {len(selection.excluded_uncertain) - 8} more"
    return [
        f"excluded: {shown} - a protected sequence carries an ambiguity code "
        "where these might cut, so none can be called safe"
    ]


def format_select_table(candidates: list, args: argparse.Namespace, selection) -> str:
    """Render the enzymes meeting a selection."""
    total = selection.considered
    excluded = _excluded_note(selection)
    if not candidates:
        wanted = (
            f"at most {args.max_cuts}"
            if args.max_cuts is not None
            else str(1 if args.cuts is None else args.cuts)
        )
        lines = [
            f"no enzyme of the {total} considered cuts {wanted} time(s)"
            + (" without cutting the protected sequences" if args.protect else "")
        ]
        # The exclusions matter most in exactly this case: without them an empty
        # result reads as "nothing cuts here" when the truth is "nothing can be
        # called safe here".
        lines.extend(excluded)
        return "\n".join(lines) + "\n"
    lines = [
        f"{len(candidates)} of {total} enzymes considered meet the selection",
        f"{'ENZYME':<12}{'SITE':<22}{'CUTS':>5}  {'END':<13}{'POSITIONS':<24}SUPPLIERS",
    ]
    for candidate in candidates:
        enzyme = candidate.enzyme
        positions = ", ".join(thousands(position) for position in candidate.positions[:3])
        if len(candidate.positions) > 3:
            positions += ", ..."
        if candidate.overhangs:
            # What the sequence actually yields, which for a non-palindromic
            # site differs between a forward and a reverse occurrence.
            end = " ".join(
                overhang_display(overhang, end_type)
                for end_type, overhang in candidate.overhangs
            )
        else:
            end = enzyme.end_type()
        suppliers = "".join(enzyme.suppliers) or "-"
        warn = " ?" if not candidate.is_clean else ""
        lines.append(
            f"{enzyme.name:<12}{enzyme.notation:<22}{candidate.cuts:>5}  "
            f"{end:<13}{positions:<24}{suppliers}{warn}"
        )
    if any(not candidate.is_clean for candidate in candidates):
        lines.append(
            "? the target carries an ambiguity code where this enzyme might also cut"
        )
    lines.extend(excluded)
    return "\n".join(lines) + "\n"


def format_select_tsv(candidates: list) -> str:
    """One TSV row per selected enzyme."""
    rows = [list(_SELECT_COLUMNS)]
    for candidate in candidates:
        enzyme = candidate.enzyme
        rows.append(
            [
                enzyme.name,
                enzyme.pattern,
                str(candidate.cuts),
                ";".join(str(position) for position in candidate.positions),
                ";".join(end for end, _ in candidate.overhangs) or enzyme.end_type(),
                ";".join(overhang for _, overhang in candidate.overhangs),
                str(enzyme.length),
                "".join(enzyme.suppliers),
            ]
        )
    return "".join("\t".join(row) + "\n" for row in rows)


def run_select(args: argparse.Namespace, stdout: IO[str]) -> int:
    """Execute the select subcommand. Returns a process exit code."""
    enzymes = _filtered_catalogue(args)
    record = read_one(args.input)
    if not record.sequence:
        raise SequenceError(f"{record.identifier} holds no sequence")

    protect = []
    for path in args.protect or []:
        records = read_records(path)
        for guarded in records:
            if not guarded.sequence:
                # Skipping it would leave the protection the user asked for
                # simply unapplied, and the run would look like a success.
                raise SequenceError(
                    f"{path}: record {guarded.identifier!r} holds no sequence, so "
                    "nothing can be protected from being cut; remove the record or "
                    "give it a sequence"
                )
            protect.append(
                (guarded.identifier, guarded.sequence, args.protect_circular)
            )

    selection = choose(
        record.sequence,
        enzymes,
        cuts=None if args.max_cuts is not None else (1 if args.cuts is None else args.cuts),
        max_cuts=args.max_cuts,
        circular=args.circular,
        name=record.identifier,
        protect=protect,
    )
    candidates = selection.candidates
    if args.limit:
        candidates = candidates[: args.limit]

    if args.stdout_format == "table":
        stdout.write(format_select_table(candidates, args, selection))
    elif args.stdout_format == "tsv":
        stdout.write(format_select_tsv(candidates))
    elif args.stdout_format == "json":
        stdout.write(
            json.dumps(
                {
                    "tool": "restrict-cut",
                    "version": __version__,
                    "rebase_release": catalogue.REBASE_RELEASE,
                    "target": {
                        "name": record.identifier,
                        "length": len(record.sequence),
                        "topology": "circular" if args.circular else "linear",
                    },
                    "protected": list(selection.protected),
                    "considered": selection.considered,
                    "excluded_protected_uncertain": list(selection.excluded_uncertain),
                    "enzymes": [
                        {
                            "name": candidate.enzyme.name,
                            "site": candidate.enzyme.pattern,
                            "notation": candidate.enzyme.notation,
                            "cuts": candidate.cuts,
                            "positions": list(candidate.positions),
                            "end_type": candidate.enzyme.end_type(),
                            "ends_observed": [
                                {"end_type": end, "overhang": overhang}
                                for end, overhang in candidate.overhangs
                            ],
                            "uncertain_positions": candidate.uncertain,
                            "suppliers": list(candidate.enzyme.supplier_names),
                        }
                        for candidate in candidates
                    ],
                },
                indent=2,
            )
            + "\n"
        )
    return 0


def format_enzyme_detail(enzyme: Enzyme) -> str:
    """Describe one enzyme in full."""
    lines = [f"{enzyme.name}  {enzyme.notation}"]
    rows: list[tuple[str, str]] = [
        ("recognition site", f"{enzyme.pattern} ({plural(enzyme.length, 'base')})"),
        (
            "specificity",
            f"{enzyme.specific_bases} fixed bases"
            + (
                f" of {enzyme.length}, the rest ambiguity codes or spacers"
                if enzyme.is_ambiguous
                else ""
            )
            + f", expected every {thousands(round(enzyme.expected_frequency))} bp"
            " of random DNA",
        ),
        ("end left", enzyme.end_type()),
    ]
    overhang = enzyme.fixed_overhang
    if overhang is None:
        rows.append(
            (
                "overhang",
                "set by the DNA, because the enzyme cuts outside its site",
            )
        )
    elif enzyme.overhang_depends_on_orientation:
        mirrored = reverse_complement(overhang)
        rows.append(
            (
                "overhang",
                f"{overhang_display(overhang, enzyme.end_type())} at a forward site, "
                f"{overhang_display(mirrored, enzyme.end_type())} at a reverse one; "
                "the site is not a palindrome",
            )
        )
    else:
        rows.append(("overhang", overhang_display(overhang, enzyme.end_type())))
    rows.append(
        (
            "cut offsets",
            ", ".join(f"top {top}, bottom {bottom}" for top, bottom in enzyme.cuts),
        )
    )
    if enzyme.excises_site:
        rows.append(("note", "cuts on both sides, excising its own recognition site"))
    if not enzyme.searches_one_strand:
        rows.append(("strands", "searched on both strands"))
    rows.append(
        ("suppliers", ", ".join(enzyme.supplier_names) if enzyme.is_available else "none")
    )

    same = catalogue.isoschizomers(enzyme)
    other = catalogue.neoschizomers(enzyme)
    if same:
        rows.append(("isoschizomers", ", ".join(same)))
    if other:
        rows.append(("same site, other cut", ", ".join(other)))

    for system in SYSTEMS.values():
        rows.append((f"{system.name} overlap", exposure(enzyme, system)))

    width = max(len(name) for name, _ in rows)
    for name, value in rows:
        lines.append(f"  {name:<{width}}  {value}")
    lines.append(
        "  methylation overlap says whether a modified base can fall inside the"
    )
    lines.append(
        "  site, not whether this enzyme is blocked by it; check the supplier"
    )
    return "\n".join(lines) + "\n"


def format_enzymes_table(enzymes: list[Enzyme]) -> str:
    """List the catalogue."""
    lines = [
        f"{len(enzymes)} enzymes",
        f"{'ENZYME':<14}{'NOTATION':<24}{'LEN':>4}{'END':>14}  {'OVERHANG':<10}SUPPLIERS",
    ]
    orientation = False
    for enzyme in enzymes:
        overhang = enzyme.fixed_overhang
        if overhang is None:
            shown = "variable"
        elif enzyme.overhang_depends_on_orientation:
            # The enzyme leaves this overhang at a forward site and its reverse
            # complement at a reverse one, so printing one of the two as though
            # it were the answer would mislead.
            shown = f"{overhang}/{reverse_complement(overhang)}"
            orientation = True
        else:
            shown = overhang or "blunt"
        lines.append(
            f"{enzyme.name:<14}{enzyme.notation:<24}{enzyme.length:>4}"
            f"{enzyme.end_type():>14}  {shown:<10}{''.join(enzyme.suppliers) or '-'}"
        )
    if orientation:
        lines.append(
            "a pair separated by / is one overhang at a forward site and the other "
            "at a reverse one"
        )
    lines.append("'variable' means the enzyme cuts outside its site, so the DNA sets the overhang")
    return "\n".join(lines) + "\n"


def format_enzymes_tsv(enzymes: list[Enzyme]) -> str:
    """One TSV row per enzyme."""
    rows = [list(_ENZYME_COLUMNS)]
    for enzyme in enzymes:
        overhang = enzyme.fixed_overhang
        rows.append(
            [
                enzyme.name,
                enzyme.pattern,
                enzyme.notation,
                str(enzyme.length),
                str(enzyme.specific_bases),
                enzyme.end_type(),
                "" if overhang is None else overhang,
                "yes" if enzyme.overhang_depends_on_orientation else "no",
                "yes" if enzyme.cuts_outside_site else "no",
                f"{enzyme.expected_frequency:.0f}",
                "".join(enzyme.suppliers),
            ]
        )
    return "".join("\t".join(row) + "\n" for row in rows)


def format_enzymes_json(enzymes: list[Enzyme]) -> str:
    """Render enzyme specificities as JSON."""
    return (
        json.dumps(
            {
                "tool": "restrict-cut",
                "version": __version__,
                "rebase_release": catalogue.REBASE_RELEASE,
                "suppliers": SUPPLIERS,
                "enzymes": [
                    {
                        "name": enzyme.name,
                        "site": enzyme.pattern,
                        "notation": enzyme.notation,
                        "cuts": [list(pair) for pair in enzyme.cuts],
                        "site_length": enzyme.length,
                        "specific_bases": enzyme.specific_bases,
                        "ambiguous_site": enzyme.is_ambiguous,
                        "end_type": enzyme.end_type(),
                        "overhang": enzyme.fixed_overhang,
                        "overhang_depends_on_orientation": (
                            enzyme.overhang_depends_on_orientation
                        ),
                        "matching_sequences": enzyme.matching_sequences,
                        "cuts_outside_site": enzyme.cuts_outside_site,
                        "excises_site": enzyme.excises_site,
                        "searches_one_strand": enzyme.searches_one_strand,
                        "expected_frequency": enzyme.expected_frequency,
                        "isoschizomers": list(catalogue.isoschizomers(enzyme)),
                        "neoschizomers": list(catalogue.neoschizomers(enzyme)),
                        "methylation_overlap": {
                            system.name: exposure(enzyme, system)
                            for system in SYSTEMS.values()
                        },
                        "suppliers": list(enzyme.supplier_names),
                    }
                    for enzyme in enzymes
                ],
            },
            indent=2,
        )
        + "\n"
    )


def format_compatible(enzymes: list[Enzyme], stdout_format: str = "table") -> str:
    """Group enzymes by the end they leave.

    Only ends shared by more than one enzyme are interchangeable, so the table
    lists those and says how many ends were unique to a single enzyme rather
    than counting them in a total it does not print.
    """
    groups = catalogue.compatible_groups(enzymes)
    shared = {end: names for end, names in groups.items() if len(names) > 1}
    unique = len(groups) - len(shared)

    if stdout_format == "none":
        return ""
    if stdout_format == "json":
        return (
            json.dumps(
                {
                    "tool": "restrict-cut",
                    "version": __version__,
                    "rebase_release": catalogue.REBASE_RELEASE,
                    "considered": len(enzymes),
                    "ends_unique_to_one_enzyme": unique,
                    "shared_ends": [
                        {"end": end, "enzymes": list(names)}
                        for end, names in shared.items()
                    ],
                },
                indent=2,
            )
            + "\n"
        )
    if stdout_format == "tsv":
        rows = [["end", "enzymes"]]
        rows.extend([end, ";".join(names)] for end, names in shared.items())
        return "".join("\t".join(row) + "\n" for row in rows)

    lines = [
        f"{len(shared)} ends are shared by more than one of the "
        f"{len(enzymes)} enzymes considered",
        "enzymes sharing a line leave ends that will ligate to one another",
    ]
    for end, names in shared.items():
        lines.append(f"  {end:<18} {', '.join(names)}")
    if unique:
        lines.append(
            f"a further {unique} end{'' if unique == 1 else 's'} are left by "
            f"only one enzyme each, so have nothing to pair with here"
        )
    lines.append(
        "enzymes cutting outside their site are omitted: their overhang comes "
        "from the DNA"
    )
    return "\n".join(lines) + "\n"


#: Flags that narrow the catalogue, and the value each has when unused. Naming
#: a single enzyme bypasses the catalogue entirely, so combining the two is a
#: contradiction rather than a preference, and is refused instead of ignored.
_FILTER_DEFAULTS = {
    "available": False,
    "supplier": None,
    "min_site": 0,
    "max_site": 0,
    "end_type": None,
    "site": None,
}


def _active_filters(args: argparse.Namespace) -> list[str]:
    """Which catalogue filters the caller actually set."""
    return [
        "--" + attribute.replace("_", "-")
        for attribute, unused in _FILTER_DEFAULTS.items()
        if getattr(args, attribute, unused) != unused
    ]


def run_enzymes(args: argparse.Namespace, stdout: IO[str]) -> int:
    """Execute the enzymes subcommand. Returns a process exit code."""
    if args.name:
        narrowing = _active_filters(args)
        if narrowing:
            raise ValueError(
                f"{', '.join(narrowing)} narrow the catalogue, which naming a "
                f"single enzyme does not search; drop the name to list the "
                f"catalogue, or the filters to describe {args.name}"
            )
        enzyme = catalogue.get(args.name)
        if args.stdout_format == "table":
            stdout.write(format_enzyme_detail(enzyme))
        elif args.stdout_format == "tsv":
            stdout.write(format_enzymes_tsv([enzyme]))
        elif args.stdout_format == "json":
            stdout.write(format_enzymes_json([enzyme]))
        return 0

    enzymes = _filtered_catalogue(args)
    if args.site:
        enzymes = [enzyme for enzyme in enzymes if enzyme.pattern == args.site]

    if args.compatible:
        stdout.write(format_compatible(enzymes, args.stdout_format))
        return 0

    if args.stdout_format == "table":
        stdout.write(format_enzymes_table(enzymes))
    elif args.stdout_format == "tsv":
        stdout.write(format_enzymes_tsv(enzymes))
    elif args.stdout_format == "json":
        stdout.write(format_enzymes_json(enzymes))
    return 0


def main(argv: Sequence[str] | None = None, stdout: IO[str] | None = None) -> int:
    """Entry point. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    out = sys.stdout if stdout is None else stdout
    try:
        if args.command == "digest":
            return run_digest(args, out)
        if args.command == "select":
            return run_select(args, out)
        if args.command == "enzymes":
            return run_enzymes(args, out)
    except (SequenceError, KeyError, OSError, ValueError) as error:
        message = error.args[0] if error.args else str(error)
        parser.exit(2, f"restrict-cut: error: {message}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
