"""Restriction enzyme mapping and digest simulation.

Coordinates are 1-based and inclusive on the forward strand. A cut is held as
the pair of nicks it leaves, each named by the base immediately 5' of it, so a
nick written ``n`` falls between base ``n`` and base ``n + 1``; keeping both is
what lets the overhang, and so the end type, be reported.
"""

from restrict_cut.digest import Cut, Digest, Fragment, digest
from restrict_cut.enzymes import (
    REBASE_RELEASE,
    Enzyme,
    all_enzymes,
    get,
    isoschizomers,
    neoschizomers,
    resolve,
)
from restrict_cut.gel import Band, lanes, unresolved_pairs
from restrict_cut.iupac import matches, reverse_complement
from restrict_cut.methylation import SYSTEMS, Overlap, exposure, site_overlaps
from restrict_cut.reader import Record, SequenceError, read_one, read_records
from restrict_cut.search import Site, find_sites
from restrict_cut.selection import Candidate, Selection, Survey, choose, survey

__version__ = "1.0.0rc1"

__all__ = [
    "Band",
    "Candidate",
    "Cut",
    "Digest",
    "Enzyme",
    "Fragment",
    "Overlap",
    "REBASE_RELEASE",
    "Record",
    "SYSTEMS",
    "Selection",
    "SequenceError",
    "Site",
    "Survey",
    "__version__",
    "all_enzymes",
    "choose",
    "digest",
    "exposure",
    "find_sites",
    "get",
    "isoschizomers",
    "lanes",
    "matches",
    "neoschizomers",
    "read_one",
    "read_records",
    "resolve",
    "reverse_complement",
    "site_overlaps",
    "survey",
    "unresolved_pairs",
]
