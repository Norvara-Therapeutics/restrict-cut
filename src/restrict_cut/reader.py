"""FASTA input.

A restriction map needs the whole molecule at once - an enzyme that cuts
outside its recognition site can reach across the end of a circular sequence -
so records are read whole rather than streamed. Plain text, gzip, bzip2 and xz
are accepted from a path or from standard input, and the compression is
recognised from the file's magic bytes rather than from its name.
"""

from __future__ import annotations

import bz2
import gzip
import io
import lzma
import sys
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Callable, Iterator

from restrict_cut.iupac import CODES


class SequenceError(ValueError):
    """Raised when input cannot be parsed as FASTA, or is not DNA."""


#: Characters permitted in a sequence: the IUPAC nucleotide codes, including
#: the ambiguity codes.
VALID_RESIDUES: frozenset[str] = frozenset("".join(CODES) + "".join(CODES).lower())

#: Deleting every valid residue leaves only the offending characters, which is
#: a C-speed check and so affordable on every line of a large file.
_RESIDUE_SIEVE = str.maketrans("", "", "".join(sorted(VALID_RESIDUES)))

#: Compression recognised from the first bytes of the stream, with the opener
#: that wraps an already-open binary stream rather than a file name.
_MAGIC: tuple[tuple[bytes, Callable[[IO[bytes]], IO[bytes]]], ...] = (
    (b"\x1f\x8b", lambda stream: gzip.GzipFile(fileobj=stream)),
    (b"BZh", bz2.BZ2File),
    (b"\xfd7zXZ\x00", lzma.LZMAFile),
)


@dataclass(frozen=True)
class Record:
    """One FASTA record.

    Attributes:
        identifier: text up to the first whitespace of the header line.
        description: the remainder of the header line, empty when absent.
        sequence: residues, uppercased, with all whitespace removed.
    """

    identifier: str
    description: str
    sequence: str

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def header(self) -> str:
        """The header line without its '>' prefix."""
        return f"{self.identifier} {self.description}".strip()

    @property
    def ambiguous_bases(self) -> int:
        """How many residues are not A, C, G or T."""
        return sum(1 for base in self.sequence if base not in "ACGT")


def _check_residues(residues: str, lineno: int) -> None:
    """Reject sequence data that is not nucleotide text.

    Without this, a file that is not a sequence at all - prose, a saved error
    page, a column of numbers - is absorbed as residues and reported with cut
    positions, which is worse than an error.
    """
    invalid = sorted(set(residues.translate(_RESIDUE_SIEVE)))
    if not invalid:
        return
    shown = ", ".join(repr(character) for character in invalid[:5])
    if len(invalid) > 5:
        shown += ", ..."
    raise SequenceError(
        f"line {lineno}: not nucleotide sequence; unexpected character(s): {shown}"
    )


def _decompressed(binary: IO[bytes]) -> IO[str]:
    """Wrap a binary stream in text mode, transparently decompressing it."""
    buffered = binary if hasattr(binary, "peek") else io.BufferedReader(binary)  # type: ignore[arg-type]
    head = bytes(buffered.peek(6)[:6])  # type: ignore[attr-defined]
    stream: IO[bytes] = buffered
    for magic, opener in _MAGIC:
        if head.startswith(magic):
            stream = opener(buffered)
            break
    return io.TextIOWrapper(stream, encoding="utf-8-sig", newline=None)


def open_text(source: str | Path | IO[str]) -> IO[str]:
    """Open a path, '-' for standard input, or pass through a text handle."""
    if hasattr(source, "read"):
        return source  # type: ignore[return-value]
    if source == "-":
        return _decompressed(sys.stdin.buffer)
    path = Path(source)  # type: ignore[arg-type]
    if path.is_dir():
        raise IsADirectoryError(f"not a file: {source}")
    if not path.is_file():
        raise FileNotFoundError(f"no such file: {source}")
    return _decompressed(path.open("rb"))


def decoded(handle: IO[str], source: object = "input") -> Iterator[str]:
    """Yield lines, reporting a damaged file as an input error.

    Decompression and decoding both happen lazily, as lines are pulled, so a
    truncated gzip or a file saved as UTF-16 does not fail when it is opened -
    it fails in the middle of reading, as an EOFError or a UnicodeDecodeError
    that has nothing to do with sequence and reaches the user as a traceback.

    The whole stream is wrapped in one ``try`` rather than guarding each line,
    because ``yield from`` hands the iteration straight to the underlying file
    object and so costs nothing per line.
    """
    label = "standard input" if source == "-" else source
    try:
        yield from handle
    except EOFError:
        raise SequenceError(
            f"{label}: compressed input ended before its end-of-stream marker, "
            "so the file is truncated"
        ) from None
    except UnicodeDecodeError as error:
        raise SequenceError(
            f"{label}: not UTF-8 text (byte {error.start} is not valid), so it is "
            "not a FASTA file; a sequence saved as UTF-16 or a binary file reads "
            "like this"
        ) from None
    except (zlib.error, lzma.LZMAError, OSError) as error:
        raise SequenceError(f"{label}: cannot be read: {error}") from None


def parse_fasta(lines: object) -> Iterator[Record]:
    """Yield records from FASTA text."""
    identifier = ""
    description = ""
    chunks: list[str] = []
    started = False

    for lineno, raw in enumerate(lines, start=1):  # type: ignore[call-overload]
        if lineno == 1:
            # A file exported on Windows may open with a byte-order mark, which
            # would otherwise make the first header unrecognisable.
            raw = raw.lstrip("\ufeff")
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if started:
                yield Record(identifier, description, "".join(chunks))
            head = line[1:].strip()
            if not head:
                raise SequenceError(f"line {lineno}: header has no identifier")
            fields = head.split(None, 1)
            identifier = fields[0]
            description = fields[1] if len(fields) > 1 else ""
            chunks = []
            started = True
            continue
        if not started:
            raise SequenceError(f"line {lineno}: sequence data before any '>' header")
        residues = "".join(line.split())
        _check_residues(residues, lineno)
        chunks.append(residues.upper())

    if not started:
        raise SequenceError("input contains no FASTA records")
    yield Record(identifier, description, "".join(chunks))


def read_records(source: str | Path | IO[str], limit: int = 0) -> list[Record]:
    """Read every FASTA record from a path, handle or '-'.

    Args:
        source: path, open text handle, or '-' for standard input.
        limit: stop after this many records; 0 reads everything.
    """
    handle = open_text(source)
    records = []
    for record in parse_fasta(decoded(handle, source)):
        records.append(record)
        if limit and len(records) >= limit:
            break
    return records


def read_one(source: str | Path | IO[str]) -> Record:
    """Read a single FASTA record, refusing a file that holds several.

    Raises:
        SequenceError: if the input holds no records or more than one.
    """
    records = read_records(source, limit=2)
    if len(records) > 1:
        raise SequenceError(
            f"expected one sequence but {source} holds several; "
            "name the record you mean with a single-record file"
        )
    return records[0]
