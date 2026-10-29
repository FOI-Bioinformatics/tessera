"""Typed genotype/lineage names for references, sourced from typing schemes.

The recombination caller competes references by accession label. This module attaches
a *typed* genotype name to each label -- e.g. ``GII.P16-GII.4``, ``CRF01_AE``,
``B.1.1.7`` -- so the report can name a parent the way virologists do, rather than by
bare accession. Names are sourced, in priority order, from (1) a user-supplied map,
(2) the lineage note carried in a fetched genome's FASTA header (NCBI datasets
``isolate.lineage`` / ``pangolinClassification``), and (3) a genotype token mined from
the GenBank title in the header.

Everything here is pure text / small-file IO with no network use, so the recomb step
stays self-contained: a standalone ``tessera recomb -i panel.msa.fasta`` can load the
sidecar map and render typed names without importing the discover/NCBI layer.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from collections.abc import Container, Iterable, Mapping
from pathlib import Path

from ..core.io import strip_sequence_extension

# A label-to-genotype map keyed by the MSA leaf label (a normalized accession).
LineageMap = dict[str, str]

# The sidecar file written next to a built panel and read back at report time.
LINEAGES_TSV = "lineages.tsv"
_TSV_HEADER = "label\tgenotype\tsource"

# A lineage designation in a title: a token carrying a digit, anchored by a letter, e.g.
# "CRF01_AE", "GII.P16-GII.1", "BA.2.10.1". The leading-letter requirement drops bare
# years ("2012") and most accessions' numeric tails while keeping real designations.
_LINEAGE_TOKEN = re.compile(r"[A-Za-z][\w.\-]*\d[\w.\-]*")
_MIN_TOKEN_LEN = 4

# Recombinant lineage names (HIV CRF/URF, Pango X-lineages): a clean parental panel
# excludes them, since they carry both parents' segments and mask the true parents.
_RECOMBINANT_LINEAGE = re.compile(r"^(CRF|URF)\d|^X[A-Z]|recombinant", re.IGNORECASE)


def is_recombinant_lineage(name: str) -> bool:
    """True when a lineage name denotes a recombinant (CRF/URF, a Pango X-lineage,
    or the literal 'recombinant')."""
    return bool(_RECOMBINANT_LINEAGE.search(name))


def _base_accession(label: str) -> str:
    """Drop the version suffix so 'U54771' and 'U54771.1' compare equal.

    Only a trailing ``.<digits>`` is a version. Splitting on the first dot instead
    matched any two labels sharing a prefix, so a strain-style name was typed by an
    unrelated one: ``strain.B.1`` took ``strain.A.1``'s lineage because both reduced
    to ``strain``.
    """
    return _VERSION_SUFFIX.sub("", label.strip())


# A trailing ".<digits>" on an accession (U54771.1 -> U54771); nothing else is one.
_VERSION_SUFFIX = re.compile(r"\.\d+$")


def organism_from_title(title: str) -> str:
    """The organism/species portion of a title (drop isolate/strain tails)."""
    return title.split(",")[0].split(" isolate ")[0].split(" strain ")[0].strip()


def dominant_lineage_token(titles: list[str], min_frac: float = 0.5) -> str | None:
    """The lineage designation shared by most titles (e.g. a query's own lineage).

    ``titles`` are descriptions WITHOUT a leading accession (as in a BLAST hit title).
    """
    counts: Counter[str] = Counter()
    for title in titles:
        for token in {t.strip("._-/") for t in _LINEAGE_TOKEN.findall(title)}:
            if len(token) >= _MIN_TOKEN_LEN:
                counts[token] += 1
    if not counts:
        return None
    token, seen = counts.most_common(1)[0]
    return token if seen >= max(2, int(min_frac * len(titles))) else None


def genotype_from_title(title: str, organism: str | None = None) -> str | None:
    """Mine a single genotype token from one reference's FASTA defline.

    ``title`` is a full defline whose first whitespace token is the accession (which is
    dropped). Returns the first designation-like token (carries a digit, anchored by a
    letter, at least ``_MIN_TOKEN_LEN`` long) in the description's leading clause, or
    ``None`` -- in which case the caller falls back to the bare accession. Tokens that
    are part of ``organism`` (e.g. an influenza subtype in the species name) are skipped.
    """
    parts = title.split(None, 1)
    desc = parts[1] if len(parts) == 2 else ""
    if not desc:
        return None
    skip: set[str] = set()
    if organism:
        skip = {t.strip("._-/") for t in _LINEAGE_TOKEN.findall(organism)}
    head = desc.split(",")[0].split(" isolate ")[0].split(" strain ")[0]
    candidates = [t.strip("._-/") for t in _LINEAGE_TOKEN.findall(head)]
    candidates = [t for t in candidates if len(t) >= _MIN_TOKEN_LEN and t not in skip]
    # Prefer the longest designation: "CRF01_AE" over the organism abbreviation "HIV-1",
    # "GII.P16-GII.4" over a stray date token. Ties keep the first (left-most).
    return max(candidates, key=len) if candidates else None


# What the pool builders write in the label slot when a genome carries no clade: an
# empty note, Nextclade's "NA", the marker given to a dataset's example sequences, and
# the value Nextclade writes for a tip it places in no clade ("unassigned"; the mumps
# dataset uses "unclassified"). These are not clades, so a genome carrying one must not
# be given a lineage of its own: on the mpox dataset "unassigned" is the reconstructed
# ancestral sequence, which is a maximum-likelihood reconstruction rather than an
# observed genome.
_NO_POOL_LABEL = frozenset({"", "na", "example", "unassigned", "unclassified"})


def label_from_pool_header(title: str) -> str | None:
    """The clade/lineage label of a pool genome's defline, verbatim.

    Tessera writes pool genomes itself -- Nextclade tree tips and NCBI Virus genomes --
    with the defline ``{accession} {label}``, so the label is everything after the
    accession and needs no mining. :func:`genotype_from_title` is for free-text GenBank
    titles: it wants a digit-bearing token of four characters or more, which drops
    ``B.1``, ``XBB``, ``D8``, ``IIb``, the HIV pure subtypes and any multi-word clade
    (``clade 2 wild-type``). Returns ``None`` for a genome that carries no label.
    """
    parts = title.split(None, 1)
    label = parts[1].strip() if len(parts) == 2 else ""
    return None if label.lower() in _NO_POOL_LABEL else label


def _is_composite_clade(clade: str, clades: Container[str]) -> bool:
    """True when ``clade`` is a ``/``-joined pair of clades that ``clades`` also holds.

    A tip labelled ``Ib/IIb`` beside standalone ``Ib`` and ``IIb`` spans both and is no
    lineage of its own; given one it wins a panel slot and competes as a donor, which on
    a clonal query reads as a mosaic. The parts must themselves be clades of the same
    pool, because a slash is also part of ordinary clade names: enterovirus D68 names a
    clade ``A2/D`` while naming neither ``A2`` nor ``D``, and Ebola names
    ``Ebov-2013/r2021``.
    """
    parts = [part.strip() for part in clade.split("/")]
    return len(parts) > 1 and all(part in clades for part in parts)


def pool_labels(files: Iterable[Path]) -> dict[str, str]:
    """``label -> clade`` for pool genome files, read from their deflines.

    A clade that is a composite of two other clades in the same pool is left out
    (:func:`_is_composite_clade`), so those genomes type as untyped and dereplicate by
    ANI rather than seeding a lineage.
    """
    out: dict[str, str] = {}
    for path in files:
        label = label_from_pool_header(first_header(path))
        if label:
            out[strip_sequence_extension(path.name)] = label
    clades = set(out.values())
    return {label: clade for label, clade in out.items()
            if not _is_composite_clade(clade, clades)}


def first_header(path: Path) -> str:
    """Return the first FASTA defline (without the leading '>'), or '' if none.

    Handles a gzip-compressed FASTA so the query's own header can be mined too.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fo:
            for line in fo:
                if line.startswith(">"):
                    return line[1:].rstrip("\n")
    except OSError:
        return ""
    return ""


def titles_from_collection(files: Iterable[Path]) -> dict[str, str]:
    """Map ``label -> full FASTA defline`` for each genome file in a collection.

    The label is ``strip_sequence_extension(file.name)`` -- the same string the aligner
    uses as the MSA leaf, so the resulting map keys match the caller's labels.
    """
    out: dict[str, str] = {}
    for path in files:
        header = first_header(path)
        if header:
            out[strip_sequence_extension(path.name)] = header
    return out


def read_lineage_rows(path: Path) -> list[tuple[str, str, str]]:
    """Read ``(label, genotype, source)`` rows from a lineages TSV (best-effort)."""
    rows: list[tuple[str, str, str]] = []
    if not Path(path).exists():
        return rows
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2 or not fields[0].strip() or not fields[1].strip():
            continue  # skip the header and any garbled line
        if fields[0] == "label" and fields[1] == "genotype":
            continue
        source = fields[2].strip() if len(fields) >= 3 else ""
        rows.append((fields[0].strip(), fields[1].strip(), source))
    return rows


def write_lineage_map(path: Path, rows: Iterable[tuple[str, str, str]]) -> None:
    """Write ``(label, genotype, source)`` rows to a lineages TSV."""
    lines = [_TSV_HEADER]
    lines += [f"{label}\t{genotype}\t{source}" for label, genotype, source in rows]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def lineage_map_from_rows(rows: Iterable[tuple[str, str, str]]) -> LineageMap:
    """Reduce ``(label, genotype, source)`` rows to a ``label -> genotype`` map."""
    return {strip_sequence_extension(label): genotype for label, genotype, _ in rows}


def load_lineage_map(path: Path | None) -> LineageMap:
    """Load a ``label -> genotype`` map from a lineages TSV (empty if absent)."""
    if path is None:
        return {}
    return lineage_map_from_rows(read_lineage_rows(Path(path)))


def build_lineage_map(
    *,
    user_tsv: Path | None = None,
    datasets_rows: Iterable[tuple[str, str]] | None = None,
    title_by_label: Mapping[str, str] | None = None,
    organism: str | None = None,
    pool_labels: Mapping[str, str] | None = None,
) -> list[tuple[str, str, str]]:
    """Merge typed names from all sources by priority into sorted sidecar rows.

    Priority (highest wins): user map > NCBI datasets lineage > pool label > mined
    title token. ``pool_labels`` are the structured labels of genomes that came from a
    pool Tessera built (see :func:`label_from_pool_header`); a genome that has one is
    never typed from a mined token. Keys are normalized with
    ``strip_sequence_extension`` so they match MSA leaves.
    """
    merged: dict[str, tuple[str, str]] = {}  # label -> (genotype, source)
    # Source 4 (lowest): a token mined from each reference's title/header note.
    for label, title in (title_by_label or {}).items():
        genotype = genotype_from_title(title, organism)
        if genotype:
            merged[strip_sequence_extension(label)] = (genotype, "title")
    # Source 3: the label a pool builder wrote into the header, taken verbatim.
    for label, clade in (pool_labels or {}).items():
        if clade:
            merged[strip_sequence_extension(label)] = (clade, "pool-header")
    # Source 2: the structured NCBI datasets lineage / pangolin classification.
    for label, lineage in datasets_rows or []:
        if lineage:
            merged[strip_sequence_extension(label)] = (lineage, "ncbi-datasets")
    # Source 1 (highest): the user-supplied override map.
    for label, genotype, _ in read_lineage_rows(user_tsv) if user_tsv else []:
        merged[strip_sequence_extension(label)] = (genotype, "user")
    return [(label, g, src) for label, (g, src) in sorted(merged.items())]


def lineage_of(label: str, lineage_map: LineageMap | None) -> str | None:
    """The genotype for ``label``, or ``None`` if untyped.

    Falls back to a version-insensitive accession match so a versioned label
    (``MK573073.1``) still resolves against an unversioned map key (``MK573073``).
    """
    if not lineage_map:
        return None
    genotype = lineage_map.get(label)
    if genotype is None:
        base = _base_accession(label)
        # Sorted, so two versions of one accession carrying different lineages resolve
        # the same way on every run rather than by dict order.
        for key in sorted(lineage_map):
            if _base_accession(key) == base:
                return lineage_map[key]
    return genotype


def typed(label: str, lineage_map: LineageMap | None) -> str:
    """Render ``label`` as ``"GENOTYPE (label)"`` when typed, else the bare ``label``."""
    genotype = lineage_of(label, lineage_map)
    return f"{genotype} ({label})" if genotype else label
