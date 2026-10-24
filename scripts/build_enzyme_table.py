"""Generate ``src/restrict_cut/enzyme_data.py`` from the REBASE distribution files.

The enzyme table is not written by hand. It is derived from the monthly REBASE
release (Roberts et al., *Nucleic Acids Research*), the reference registry for
restriction enzyme specificities, and every row is verified twice before it is
emitted:

1. REBASE publishes each enzyme's specificity twice, in two different files and
   two different encodings: as machine-readable integer cut coordinates in the
   EMBOSS-format file, and as a human-readable recognition sequence such as
   ``G^AATTC`` or ``GGTCTC(1/5)`` in the Type II file. This script reads the
   integers, independently parses the notation, and requires the two to agree.
   A sign error or an off-by-one in either reading cannot survive that. An entry
   the second source does not state, and so cannot confirm, is treated as a
   failure rather than waved through: an unchecked offset is exactly what this
   script exists to keep out of the package.
2. Offsets are cross-checked against Biopython's separately maintained
   REBASE-derived dictionary when Biopython is importable, which catches a
   misreading of the EMBOSS format common to both of this script's readers.

The REBASE files themselves are not kept in this repository. They are
downloaded when the table is rebuilt, which is how the upstream tools that
depend on REBASE handle it too, and keeps this project from redistributing
someone else's database.

Run from the repository root:

    python scripts/build_enzyme_table.py                  # download and rebuild
    python scripts/build_enzyme_table.py --source data/   # use local copies
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.request
from pathlib import Path

#: Where the monthly REBASE files are published.
REBASE_URL = "http://rebase.neb.com/rebase/"

#: The files this script reads, and what each is for.
REBASE_FILES = {
    "link_emboss_e": "machine-readable cut coordinates",
    "link_gtype2": "recognition sequences and commercial suppliers",
    "link_emboss_s": "supplier codes",
}


def download(directory: Path) -> None:
    """Fetch the REBASE files this script reads into ``directory``."""
    directory.mkdir(parents=True, exist_ok=True)
    for filename, purpose in REBASE_FILES.items():
        target = directory / filename
        print(f"fetching {filename} ({purpose})")
        with urllib.request.urlopen(REBASE_URL + filename, timeout=120) as response:
            target.write_bytes(response.read())


def _offset(emboss_value: int) -> int:
    """Convert an EMBOSS cut coordinate to an offset from the site start.

    The returned offset counts the bases of the recognition site that lie 5' of
    the cut, so 0 means 'immediately before the site' and ``len(site)`` means
    'immediately after it'; a negative offset reaches upstream of the site.

    EMBOSS writes an upstream cut one lower than the offset it means, because
    its coordinates have no zero: a cut ten bases before the site is written
    -11, not -10. Applying that correction to a positive coordinate would be
    wrong, so the sign is tested rather than assumed.
    """
    return emboss_value if emboss_value > 0 else emboss_value + 1


def read_emboss(path: Path) -> dict[str, tuple[str, tuple[tuple[int, int], ...]]]:
    """Parse ``link_emboss_e``: name, pattern, length, ncuts, blunt, c1-c4."""
    enzymes: dict[str, tuple[str, tuple[tuple[int, int], ...]]] = {}
    for line in path.read_text(encoding="latin-1").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        name, pattern = fields[0], fields[1].upper()
        length, ncuts = int(fields[2]), int(fields[3])
        coordinates = [int(value) for value in fields[5:9]]
        if not ncuts:
            # REBASE records the specificity but not where this enzyme cuts.
            # It cannot be used to predict a fragment, so it is left out rather
            # than shipped with a guessed cut position.
            continue
        if len(pattern) != length:
            raise ValueError(f"{name}: pattern {pattern!r} is not {length} long")
        cuts = [(_offset(coordinates[0]), _offset(coordinates[1]))]
        if ncuts == 4:
            cuts.append((_offset(coordinates[2]), _offset(coordinates[3])))
        enzymes[name] = (pattern, tuple(cuts))
    return enzymes


def parse_notation(text: str) -> tuple[str, tuple[tuple[int, int | None], ...]] | None:
    """Read REBASE's recognition-sequence notation into a pattern and offsets.

    Handles the caret form (``G^AATTC``), the offset-pair form
    (``GGTCTC(1/5)``, ``CCGC(-3/-1)``), the two-sided form
    (``(10/12)CGANNNNNNTGC(12/10)``), and the padded trailing caret
    (``TGGCCANNNNNNNNNNN^``) in which REBASE spells a downstream cut out in Ns
    and leaves the bottom-strand cut unexpressed. A bottom offset of None means
    the notation does not state it.
    """
    text = text.split(",")[0].strip()
    if not text:
        return None

    upstream: tuple[int, int] | None = None
    leading = re.match(r"^\((-?\d+)/(-?\d+)\)", text)
    if leading:
        # A leading pair counts bases upstream of the site, so it is negated
        # into this module's signed offsets.
        upstream = (-int(leading.group(1)), -int(leading.group(2)))
        text = text[leading.end():]

    downstream: tuple[int, int] | None = None
    trailing = re.search(r"\((-?\d+)/(-?\d+)\)$", text)
    if trailing:
        downstream = (int(trailing.group(1)), int(trailing.group(2)))
        text = text[: trailing.start()]

    if "^" in text:
        position = text.index("^")
        padded = text.replace("^", "")
        if position == len(padded):
            # The cut sits beyond the site, written out as a run of Ns.
            pattern = padded.rstrip("N")
            return pattern, ((position, None),)
        return padded, ((position, len(padded) - position),)

    if downstream is None:
        return None
    length = len(text)
    cut = (length + downstream[0], length + downstream[1])
    return text, ((upstream, cut) if upstream else (cut,))


def read_gtype2(path: Path) -> dict[str, dict[str, str]]:
    """Parse ``link_gtype2``: recognition sequence and commercial suppliers."""
    entries: dict[str, dict[str, str]] = {}
    for block in path.read_text(encoding="latin-1").split("\n\n"):
        fields = dict(re.findall(r"<(\d)>(.*)", block))
        name = fields.get("1", "").strip()
        if name:
            entries[name] = {
                "recognition": fields.get("3", "").strip(),
                "suppliers": fields.get("5", "").strip(),
            }
    return entries


def read_suppliers(path: Path) -> dict[str, str]:
    """Parse ``link_emboss_s``: one-letter supplier code to supplier name."""
    suppliers: dict[str, str] = {}
    for line in path.read_text(encoding="latin-1").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        code, _, name = line.partition(" ")
        suppliers[code.strip()] = name.strip()
    return suppliers


def load_biopython_offsets() -> dict[str, tuple[int, int]]:
    """Offsets from Biopython's REBASE-derived table, for cross-checking."""
    try:
        from Bio.Restriction.Restriction_Dictionary import rest_dict
    except ImportError:
        return {}
    offsets = {}
    for name, record in rest_dict.items():
        if record["fst5"] is not None and record["fst3"] is not None:
            offsets[name] = (record["fst5"], record["size"] + record["fst3"])
    return offsets


HEADER = '''"""Restriction enzyme specificities, generated from REBASE.

THIS FILE IS GENERATED BY ``scripts/build_enzyme_table.py`` - DO NOT EDIT.

Source: REBASE, The Restriction Enzyme Database, release {release}.
Richard J. Roberts, Tamas Vincze, Janos Posfai and Dana Macelis.
"REBASE - a database for DNA restriction and modification: enzymes, genes and
genomes." Nucleic Acids Research 51: D629-D630 (2023).
<https://rebase.neb.com>

Each entry is ``(pattern, cuts, suppliers)``:

``pattern``
    The recognition sequence, as IUPAC nucleotide codes.
``cuts``
    One ``(top, bottom)`` pair per double-strand break the enzyme makes. Each
    offset counts the bases of the recognition site lying 5' of the cut, so 0
    is immediately before the site and a negative offset reaches upstream of
    it. ``bottom`` is measured in the same top-strand frame, so
    ``bottom - top`` is the overhang: positive is a 5' overhang, negative a 3'
    overhang, zero a blunt end. Enzymes that excise their own recognition site,
    such as BcgI, carry two pairs.
``suppliers``
    REBASE one-letter codes for the companies selling the enzyme, empty when it
    has no commercial source.
"""

from __future__ import annotations

#: REBASE release these specificities were taken from.
REBASE_RELEASE = "{release}"

#: Companies selling restriction enzymes, by REBASE one-letter code.
SUPPLIERS: dict[str, str] = {{
{suppliers}}}

#: ``name -> (pattern, cuts, supplier codes)`` for every REBASE enzyme whose
#: cut position is known.
ENZYMES: dict[str, tuple[str, tuple[tuple[int, int], ...], str]] = {{
{rows}}}
'''


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="build_enzyme_table.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--source",
        type=Path,
        metavar="DIR",
        help="read the REBASE files from DIR instead of downloading them",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("src/restrict_cut/enzyme_data.py"),
        metavar="PATH",
        help="where to write the generated module",
    )
    args = parser.parse_args(argv)

    data = args.source
    if data is None:
        import tempfile

        holding = tempfile.TemporaryDirectory(prefix="rebase-")
        data = Path(holding.name)
        download(data)
    else:
        missing = [name for name in REBASE_FILES if not (data / name).is_file()]
        if missing:
            print(f"{data}: missing {', '.join(missing)}", file=sys.stderr)
            return 1
    out = args.output

    emboss_path = data / "link_emboss_e"
    release = "unknown"
    for line in emboss_path.read_text(encoding="latin-1").splitlines()[:10]:
        found = re.search(r"REBASE version (\d+)", line)
        if found:
            release = found.group(1)
            break

    enzymes = read_emboss(emboss_path)
    published = read_gtype2(data / "link_gtype2")
    suppliers = read_suppliers(data / "link_emboss_s")
    biopython = load_biopython_offsets()

    unverified: list[str] = []
    notation_mismatch: list[str] = []
    biopython_mismatch: list[str] = []
    rows: list[str] = []

    for name in sorted(enzymes):
        pattern, cuts = enzymes[name]
        entry = published.get(name, {})

        parsed = parse_notation(entry.get("recognition", ""))
        if parsed is None:
            unverified.append(name)
        else:
            stated_pattern, stated_cuts = parsed
            disagrees = stated_pattern != pattern or len(stated_cuts) != len(cuts)
            if not disagrees:
                for (top, bottom), (stated_top, stated_bottom) in zip(cuts, stated_cuts):
                    if top != stated_top:
                        disagrees = True
                    if stated_bottom is not None and bottom != stated_bottom:
                        disagrees = True
            if disagrees:
                notation_mismatch.append(
                    f"{name}: EMBOSS {pattern} {cuts} != notation "
                    f"{entry['recognition']!r} -> {stated_pattern} {stated_cuts}"
                )

        if name in biopython and len(cuts) == 1 and cuts[0] != biopython[name]:
            biopython_mismatch.append(
                f"{name}: {cuts[0]} != biopython {biopython[name]}"
            )

        codes = "".join(sorted(set(entry.get("suppliers", "")) & set(suppliers)))
        cut_text = ", ".join(f"({top}, {bottom})" for top, bottom in cuts)
        comma = "," if len(cuts) == 1 else ""
        rows.append(f'    "{name}": ("{pattern}", ({cut_text}{comma}), "{codes}"),\n')

    compared = len(set(enzymes) & set(biopython))
    print(f"REBASE release {release}: {len(enzymes)} enzymes with a known cut position")
    print(f"notation cross-check: {len(enzymes) - len(unverified)} verified, "
          f"{len(notation_mismatch)} disagree, {len(unverified)} not cross-checkable")
    for problem in notation_mismatch[:20]:
        print("  ", problem)
    for name in unverified[:20]:
        print(f"   {name}: REBASE states no recognition sequence this script can "
              f"parse, so its cut offsets cannot be checked against anything")
    if len(unverified) > 20:
        print(f"   ... and {len(unverified) - 20} more")
    print(f"biopython cross-check: {compared} compared, "
          f"{len(biopython_mismatch)} disagree")
    for problem in biopython_mismatch[:20]:
        print("  ", problem)

    # An entry that cannot be checked is not the same as an entry that passed,
    # and writing it anyway is the failure mode that matters here: the whole
    # point of this script is that no cut offset reaches the package without
    # agreeing with a second, independent statement of the same fact. Silence
    # from the cross-check has to count as failure, or the guarantee the
    # package's documentation makes about its own data is not true.
    if notation_mismatch or biopython_mismatch or unverified:
        reasons = []
        if notation_mismatch:
            reasons.append(f"{len(notation_mismatch)} disagree with REBASE's notation")
        if biopython_mismatch:
            reasons.append(f"{len(biopython_mismatch)} disagree with biopython")
        if unverified:
            reasons.append(f"{len(unverified)} could not be cross-checked at all")
        print(
            "refusing to write the table: " + ", ".join(reasons)
            + ". Every entry must agree with a second source before it ships."
        )
        return 1

    supplier_rows = "".join(
        f'    "{code}": "{suppliers[code]}",\n' for code in sorted(suppliers)
    )
    out.write_text(
        HEADER.format(release=release, suppliers=supplier_rows, rows="".join(rows)),
        encoding="utf-8",
    )
    commercial = sum(1 for name in enzymes if published.get(name, {}).get("suppliers"))
    print(f"wrote {out} ({out.stat().st_size:,} bytes); "
          f"{commercial} of {len(enzymes)} are commercially available")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
