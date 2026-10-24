"""Re-fetch the example sequences from NCBI.

Every example in this directory is a real database record, downloaded
unmodified. None of it is simulated, and nothing is sliced or edited, so a
result computed from one of these files can be checked against the public
record.

    python examples/fetch_examples.py

Set a contact address in the ``NCBI_EMAIL`` environment variable if you want
NCBI to be able to reach you about the requests; E-utilities works without one.
"""

from __future__ import annotations

import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

#: ``filename -> (accession, what it is, expected length)``. The lengths are
#: asserted after download, so a silently truncated or superseded record is
#: caught rather than quietly replacing a known example.
RECORDS: dict[str, tuple[str, str, int]] = {
    "puc19.fasta": ("L09137.2", "cloning vector pUC19", 2_686),
    "pbr322.fasta": ("J01749.1", "cloning vector pBR322", 4_361),
    "lambda.fasta": ("NC_001416.1", "Enterobacteria phage lambda", 48_502),
    "aav2.fasta": ("NC_001401.2", "adeno-associated virus 2", 4_679),
}

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fetch(accession: str, email: str | None) -> str:
    """Download one record as FASTA."""
    query = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": "restrict-cut-examples",
    }
    if email:
        query["email"] = email
    with urllib.request.urlopen(f"{EFETCH}?{urllib.parse.urlencode(query)}", timeout=60) as response:
        return response.read().decode("utf-8")


def main() -> int:
    email = os.environ.get("NCBI_EMAIL") or None
    directory = Path(__file__).parent
    failures = 0

    for filename, (accession, description, expected) in RECORDS.items():
        try:
            text = fetch(accession, email)
        except OSError as error:
            print(f"{filename}: download failed: {error}", file=sys.stderr)
            failures += 1
            continue
        sequence = "".join(line.strip() for line in text.splitlines()[1:])
        if len(sequence) != expected:
            print(
                f"{filename}: {accession} is {len(sequence):,} bases, expected "
                f"{expected:,}; the record may have been revised",
                file=sys.stderr,
            )
            failures += 1
            continue
        (directory / filename).write_text(text, encoding="utf-8")
        print(f"{filename}: {accession}, {len(sequence):,} bp, {description}")
        time.sleep(0.4)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
