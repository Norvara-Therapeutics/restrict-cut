# restrict-cut

Restriction mapping from the command line: where the enzymes cut, what fragments
that leaves, which enzyme will linearise a vector without touching the insert,
and whether the sites you are relying on are methylated before you ever open the
tube. Pure Python with no third-party dependencies, and a self-contained HTML
report that opens from a local disk without a network connection.

Specificities come from [REBASE](https://rebase.neb.com), the reference registry
of restriction enzymes: 720 enzymes whose cut position is known, 581 of them
commercially available.

## Install

Install the prerelease from PyPI:

```bash
python -m pip install --pre norvara-restrict-cut
```

The Python package and command keep the shorter project name:

```python
import restrict_cut
```

```bash
restrict-cut --version
```

## Quickstart

```bash
restrict-cut digest -i examples/puc19.fasta -e EcoRI,PstI,PvuII --circular --no-gel
```

```
> L09137.2  2,686 bp  circular
  EcoRI        G^AATTC                  1 cut      5'-AATT
  PstI         CTGCA^G                  1 cut      3'-TGCA
  PvuII        CAG^CTG                  2 cuts     blunt

         | ||   |
  o==============================================================o
  1                                                          2,686
    | PvuII at 308
    | EcoRI at 396
    | PstI at 439
    | PvuII at 630

  fragments  4 fragments, 2,686 bp total
        2,364 bp  631..308  PvuII / PvuII spans origin
          191 bp  440..630  PstI / PvuII
           88 bp  309..396  PvuII / EcoRI
           43 bp  397..439  EcoRI / PstI

  methylation  no site carries a methylated base (dam, dcm)
```

The other question a map does not answer is which enzyme to use. `select`
searches the whole catalogue instead of a named set, and `--protect` names
sequence that must survive untouched:

```bash
restrict-cut select -i examples/puc19.fasta --circular --cuts 1 \
    --protect examples/aav2.fasta --available --min-site 6 --limit 8
```

```
8 of 421 enzymes considered meet the selection
ENZYME      SITE                   CUTS  END          POSITIONS               SUPPLIERS
SbfI        CCTGCA^GG                 1  3'-TGCA      439                     INSV
SdaI        CCTGCA^GG                 1  3'-TGCA      439                     BV
XbaI        T^CTAGA                   1  5'-CTAG      423                     BIJMNQRSVX
SspI        AAT^ATT                   1  blunt        2,503                   BIJNQRSVX
NdeI        CA^TATG                   1  5'-TA        184                     BJNQRSVX
SphI        GCATG^C                   1  3'-CATG      445                     IJNQRSVX
```

Every one of those linearises pUC19 and leaves the AAV2 genome intact. Rarer
sites come first, and within a tier the enzymes more companies sell.

Protection is treated as a guarantee, so it is deliberately stricter than the
target side: an enzyme is excluded if the protected sequence carries a confirmed
site **or** an ambiguity code where a site cannot be ruled out. A guarded
sequence reading `GANTTC` has no confirmed EcoRI site and cannot be called free
of one either, so EcoRI is named as excluded rather than returned as safe. The
consequence is worth knowing: protecting a sequence full of Ns will exclude a
great many enzymes, which is the honest answer rather than a fault.

## Commands

| Command | What it does |
| --- | --- |
| `digest` | cut one sequence with named enzymes; sites, cuts, fragments, map, sizes |
| `select` | search the catalogue for enzymes cutting a target a set number of times |
| `enzymes` | browse the catalogue, or describe one enzyme in full |

```
restrict-cut enzymes BsaI
```

```
BsaI  GGTCTC(1/5)
  recognition site  GGTCTC (6 bases)
  specificity       6 fixed bases, expected every 2,048 bp of random DNA
  end left          5' overhang
  overhang          set by the DNA, because the enzyme cuts outside its site
  cut offsets       top 7, bottom 11
  strands           searched on both strands
  isoschizomers     Bso31I, Eco31I
  dam overlap       never
  dcm overlap       possible
  cpg overlap       possible
```

## Options

```
digest
-i, --input PATH         FASTA file holding the sequence, or '-' for stdin
-e, --enzyme NAMES       enzymes to apply, comma-separated or repeated
    --all                apply every enzyme in the catalogue, subject to the
                         filters, and report those with a cut or a site
-c, --circular           treat the sequence as a closed circle, as a plasmid is
    --methylation NAMES  activities to check: dam, dcm, cpg, or none
                         (default: dam,dcm)
    --resolution PERCENT call two fragments unresolved below this difference
                         (default: 5)
    --ladder SIZES       comma-separated sizes to show as a ladder
    --map-width N        width of the text map (default: 64)
    --no-map, --no-gel   drop a section from the summary

select
-n, --cuts N             required number of cuts in the target (default: 1)
    --max-cuts N         accept any count up to N instead of an exact one
                         (alternative to --cuts, not combinable with it)
-p, --protect PATH       sequence that must not be cut at all; repeatable
                         (an enzyme that *might* cut it is excluded too)
    --limit N            show only the first N enzymes

enzymes [NAME]
    --site SEQUENCE      only enzymes recognising exactly this sequence
    --compatible         group enzymes by the end they leave

shared by all three
    --available          only enzymes a company sells
    --supplier CODE      only enzymes sold by this REBASE supplier
    --min-site N         only sites at least N bases long
    --max-site N         only sites at most N bases long
    --end-type TYPE      only enzymes leaving a 5' overhang, blunt or 3' overhang
```

## Outputs

| Flag | Output |
| --- | --- |
| `--format table` (default) | the readable summary on stdout |
| `--format tsv` | one row per fragment (`digest`) or per enzyme (`select`, `enzymes`) |
| `--format json` \| `none` | every site, cut, fragment and band as JSON, or nothing |
| | `digest --all` records `enzymes_applied` beside `enzymes_reported` |
| `--output PATH` | one TSV row per fragment |
| `--cut-table PATH` | one TSV row per cut, with nick positions and overhang |
| `--report PATH` | the self-contained HTML report |

Any combination may be requested in one pass. TSV and JSON carry exact values
and full overhang sequences; only the terminal summary abbreviates.

An output path that is also an input, or that another output already claims, is
refused before anything is written. The sequence is read into memory before the
outputs are produced, so without that check `--output` aimed at the input file
would destroy the source and still report success. Sameness is decided by device
and inode rather than by comparing paths, so a hard link or a symbolic link to an
input is caught as well: those are different paths naming one file, and a path
comparison lets them through.

## Reading the output

- Coordinates are 1-based and inclusive on the forward strand. A cut is held as
  the two nicks it leaves, each named by the base 5' of it: a nick written `n`
  falls between base `n` and base `n + 1`.
- **On a circular molecule, `n` cuts give `n` fragments**; on a linear one they
  give `n + 1`. A plasmid cut once yields one full-length linear fragment, not
  two pieces. A fragment whose end is smaller than its start spans the origin.
- The count that sets the number of fragments is the number of *distinct
  top-strand nicks*, not the number of cuts. Two enzymes can break the same bond
  while nicking the bottom strand in different places — AccII and Alw26I both cut
  after base 3 of pUC19, one bluntly and one leaving a four-base overhang — which
  is one boundary made by two cuts. Both enzymes are named on the fragments they
  bound, joined by `+`, and the summary says so rather than leaving a cut count
  above the fragment count looking like a contradiction.
- Fragment lengths are measured between top-strand nicks. For a sticky end the
  two strands differ in length, so some convention is needed; this one makes the
  fragment lengths sum to the length of the molecule, which is printed so you can
  check it.
- An overhang is written `5'-GATC` or `3'-TGCA`, never bare. Two ends anneal when
  both the bases and the strand agree, and a 5' GATC will not ligate to a 3' GATC.
- `enzymes --compatible` groups the enzymes whose ends will ligate to each other,
  so BamHI, BglII and BclI appear together. It lists only ends shared by more
  than one enzyme, since an end unique to a single enzyme has nothing to pair
  with, and says how many such ends it left out. Two kinds of enzyme are omitted
  entirely, for the same reason: the end is not a property of the enzyme. One
  that cuts outside its site, such as BsaI, takes its overhang from whatever DNA
  sits there. One with a sticky end and a non-palindromic site, such as BauI,
  leaves `5'-ACGA` where its site reads forward and `5'-TCGT` where it reads
  backward — and those do not anneal to each other, so pairing it with BsiI on
  the strength of the forward reading would promise something half the sites
  will not deliver. The list view writes such an overhang as `ACGA/TCGT`, and
  `select` reports the end read from the sequence in hand.
- Site rarity counts the sequences an enzyme actually recognises, over both
  strands. Doubling the one-strand count is wrong when a pattern overlaps its
  own reverse complement: Hin4I's `GAYNNNNNVTC` does, so doubling would call it
  one site every 341 bp when counting distinct sequences gives 512 bp.
- Site length is not specificity. AhdI's `GACNNNNNGTC` is eleven bases but only
  six are fixed, so it is no rarer than EcoRI's `GAATTC`. `select` ranks on
  expected frequency, which accounts for both ambiguity codes and spacer Ns.
- An ambiguity code in the *sequence* is never reported as a site. Those
  positions are counted separately as uncalled, because sending someone to the
  bench with a digest that cannot work is worse than reporting nothing. In a
  `--protect` sequence the same uncertainty excludes the enzyme outright, since
  there the question is whether the insert is safe, not whether the site is real.
- A Type IIS site near the end of a linear molecule is reported as present but
  not cutting, because the cut would fall off the end. Make the molecule
  `--circular` and it cuts. An enzyme that cuts on *both* sides of its site, such
  as BcgI, can land one break and lose the other; that site is reported
  separately as cut but not excised, rather than as uncut.
- A ladder and the sample share one size axis, so a band level with a ladder
  band really is that size. Each band gets its own row even where several sizes
  round to the same place, because a lane that silently drops a band is worse
  than a crowded one.

## Methylation

Plasmid grown in an ordinary *E. coli* strain is methylated before you see it.
With `--methylation dam,dcm` (the default) every recognition site is checked
against the sequence and its context, and a site carrying a methylated base is
reported.

What this reports is **overlap**, not blocking, and the difference matters:

- BamHI's `GGATCC` always contains a Dam `GATC`, and BamHI cuts it anyway.
- MboI's site *is* `GATC`, and Dam methylation stops it outright.
- DpnI needs the methylation before it will cut at all.

Whether an enzyme is blocked is an empirical property of that enzyme, and no
sequence can settle it. What a sequence can settle is whether *this* site, in
*this* construct, carries a modified base — which is the part you cannot look
up. `enzymes NAME` reports the static picture as `always`, `possible` or
`never`, derived from the recognition sequence rather than from a hand-kept
list; check the enzyme against its supplier's sensitivity table before acting.

## By size

`digest` orders the fragments on a logarithmic size axis and marks the bands
that will not come apart, because a digest that looks informative as a table of
lengths can be unreadable on a gel: two fragments of equal length run as one
band, and a 40 bp difference is plain at 200 bp and invisible at 10 kb.

This is an ordering by size, **not a simulation of electrophoresis**. No agarose
percentage, buffer, voltage or run time is modelled, and no migration distance is
claimed. `--resolution` sets the percentage below which two bands are called
unresolved, and is a rule of thumb rather than a measurement.

## Limitations

- Input is FASTA, so topology comes from `--circular` rather than from the file.
  Plain text, gzip, bzip2 and xz are accepted and recognised from the file's
  magic bytes rather than its name. A truncated or corrupt archive, and a file
  saved as UTF-16 rather than UTF-8, are reported as input errors naming the
  file rather than as a decompression traceback.
- Only complete digests are modelled. Partial digests, star activity, enzyme
  buffer compatibility, double-digest buffer choice and incubation conditions are
  out of scope: those belong to the supplier's tables, not to a sequence.
- The catalogue holds the 720 REBASE enzymes whose cut position is known. The
  roughly 330 further entries whose specificity is recorded but whose cut site is
  not are left out rather than shipped with a guessed position, so an enzyme
  missing from `enzymes` is missing deliberately.
- Isoschizomers are reported from within that catalogue, so an enzyme whose only
  relatives have no known cut position lists none.
- Methylation covers Dam, Dcm and CpG. Other methyltransferases, and the
  host-specificity systems, are not modelled.
- Nothing here designs primers, picks a buffer, or orders anything.

## Examples

Every file in `examples/` is an unmodified database record, not simulated
sequence: pUC19 (`L09137.2`), pBR322 (`J01749.1`), phage lambda
(`NC_001416.1`) and AAV2 (`NC_001401.2`). `examples/fetch_examples.py` re-downloads
them and checks each length, so a superseded record is caught rather than
quietly replacing a known example.

One worked check is worth stating, because it is a difference rather than a
match. Digesting the deposited lambda sequence with HindIII gives seven
fragments — 23130, 9416, 6682, 4361, 2322, 2027 and 564 — and not the eight of
the commercial lambda size standard, which includes a 6557 and a 125. Those two
sum to 6682, and reaching them needs a seventh HindIII site near position 37584;
the deposited sequence reads `AAGCTC` there, not `AAGCTT`. The tool reports what
is in the file it was given.

## Where the enzyme table comes from

`src/restrict_cut/enzyme_data.py` is generated, never edited by hand. REBASE
publishes each specificity twice, in two files and two encodings: as integer cut
coordinates in the EMBOSS-format file, and as a recognition sequence such as
`G^AATTC` or `GGTCTC(1/5)` in the Type II file. `scripts/build_enzyme_table.py`
reads the integers, independently parses the notation, and refuses to write the
table unless the two agree for every enzyme — which they do, for all 720. An
entry the second source does not state, and so cannot confirm, counts as a
failure too: silence is not agreement, and an unchecked cut offset is exactly
what this script exists to keep out of the package. Where Biopython is installed
it is cross-checked against that separately maintained REBASE-derived table as
well.

The REBASE files are downloaded when the table is rebuilt rather than kept in
this repository:

```bash
python scripts/build_enzyme_table.py
```

If you use this tool in published work, cite REBASE: Roberts, Vincze, Posfai and
Macelis, "REBASE — a database for DNA restriction and modification: enzymes,
genes and genomes", *Nucleic Acids Research* 51: D629–D630 (2023).

## Development

```bash
pip install -e ".[dev]"
pytest -q
python -m pyflakes src/restrict_cut tests scripts
```

The suite includes a brute-force differential test: a deliberately naive
reference implementation, sharing no code with the library, is run over the
whole catalogue against real plasmids, random sequence and sequence full of
ambiguity codes, on both topologies, and every cut position must agree. It
found the bug it was written to find — recognition patterns were compiled into
character classes of *bases*, so an `N` in the sequence failed to match an `N`
in the site. `tests/test_circular.py` rotates a plasmid through eight origins
and requires every fragment length and overhang to come back unchanged, which is
what caught the overhang of a cut spanning the origin being computed from
wrapped coordinates.

## License

Licensed under the
[PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0).
See [LICENSE](LICENSE) for the complete terms. Enzyme specificities are derived
from REBASE, copyright Richard J. Roberts, and are redistributed here in
generated form with attribution.
