"""Genome staging.

progressiveMauve and the other aligners cannot read gzip, and need plain FASTA
files in a working directory. :func:`stage_genomes` populates a directory by
decompressing ``.gz`` inputs (via the :mod:`gzip` module, no shell) and
symlinking plain inputs, normalizing every file to a ``.fasta`` extension so the
leaf name in the resulting MSA is a clean genome label.
"""

from __future__ import annotations

import gzip
import logging
import re
import shutil
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import BinaryIO, TextIO

from .errors import UserInputError

# Extensions stripped to derive a genome's label (its leaf name in the MSA).
SEQUENCE_EXTENSIONS = (".fna", ".fasta", ".fa")


def _open_text(path: str | Path) -> TextIO:
    """Open a plain or gzip-compressed text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA file into a list of ``(header-first-token, sequence)``.

    A ``.gz`` input is read through :mod:`gzip`: staging already accepts compressed
    genomes, so the readers that look inside a query must accept them too.

    Whitespace inside a sequence line is dropped: a trailing space or tab is not a
    base, and kept it lengthens the sequence (the backbone row then no longer matches
    the rows aligned to it). A header with no name -- ``>`` or ``> `` -- reads as an
    unnamed record (``""``).
    """
    records: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    with _open_text(path) as fo:
        for line in fo:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                tokens = line[1:].split()
                name = tokens[0] if tokens else ""
                seq = []
            else:
                seq.append("".join(line.split()))
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def write_fasta_record(handle: TextIO, name: str, seq: str, width: int = 80) -> None:
    """Write one ``>name`` record to ``handle``, wrapping the sequence at ``width``."""
    handle.write(f">{name}\n")
    for pos in range(0, len(seq), width):
        handle.write(seq[pos : pos + width] + "\n")


def normalize_reference(
    genomes: Sequence[Path],
    reference: Path | None,
    *,
    tool: str,
    min_genomes: int = 1,
    ensure_member: bool = False,
) -> tuple[list[Path], Path]:
    """Resolve the backbone reference and validate the genome set for an aligner.

    The reference defaults to the first genome. With ``ensure_member`` an explicit
    reference not already in ``genomes`` is prepended (the pairwise backends align
    every query against it, so it must be present). Raises :class:`UserInputError`
    when fewer than ``min_genomes`` genomes remain.
    """
    genomes = list(genomes)
    if not genomes:
        raise UserInputError(f"{tool} alignment needs at least {min_genomes} genomes.")
    if reference is None:
        reference = genomes[0]
    elif ensure_member and reference not in genomes:
        genomes = [reference, *genomes]
    if len(genomes) < min_genomes:
        raise UserInputError(
            f"{tool} alignment needs at least {min_genomes} genomes."
        )
    return genomes, reference


def strip_sequence_extension(name: str) -> str:
    """Return ``name`` without a trailing ``.gz`` and/or sequence extension."""
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    for ext in SEQUENCE_EXTENSIONS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def safe_filename_stem(name: str, fallback: str = "sequence") -> str:
    """A file or directory name derived from an untrusted label (a FASTA header).

    Path separators and other punctuation become ``_``, and leading/trailing dots are
    removed so the result can never be ``.`` or ``..`` or climb out of the directory it
    is joined to. ``fallback`` is returned when nothing usable is left.
    """
    return re.sub(r"[^\w.-]+", "_", name).strip(".") or fallback


def collection_genomes(directory: Path) -> list[Path]:
    """The genome files of a collection directory, sorted by name.

    Hidden files are skipped. A collection is "every file in this directory", and a
    directory a person has opened in a file browser holds files they never put there
    (``.DS_Store`` on macOS): staged as a genome, one becomes the backbone or crashes
    the aligner. Anything else that is not FASTA is rejected by name at staging.
    """
    return sorted(
        p for p in Path(directory).iterdir() if p.is_file() and not p.name.startswith(".")
    )


# Written beside a working collection when Tessera creates it, so a later run can tell
# its own scratch directory from one the user happened to name the same.
_WORKING_COPY_MARKER = ".tessera-working-collection"
# Files a run leaves beside its working copy. Output directories written before the
# marker existed are recognised by these, so re-running into one still works.
_RUN_ARTEFACTS = ("round1.msa.fasta", "panel.msa.fasta", "fill_summary.tsv", "panel_lineages.tsv")


def _clear_working_copy(dest: Path) -> None:
    """Remove a previous run's working collection at ``dest``.

    The working copy is cleared at the start of every run. That is only safe for a
    directory Tessera made: ``<output>/collection`` is also a natural place for a person
    to keep their own genomes, and clearing that would destroy them. An existing,
    non-empty directory is removed only when a previous run left its marker (or, for
    output written by an older release, its other files) beside it; otherwise refuse.
    """
    if not dest.exists():
        return
    parent = dest.parent
    ours = (parent / _WORKING_COPY_MARKER).exists() or any(
        (parent / name).exists() for name in _RUN_ARTEFACTS
    )
    if not ours and (not dest.is_dir() or any(dest.iterdir())):
        raise UserInputError(
            f"{dest} already exists and was not created by Tessera. A run clears that "
            "directory to hold its working copy of the references, which would delete "
            "what is in it. Move it elsewhere or choose a different output directory."
        )
    shutil.rmtree(dest)


def _mark_working_copy(dest: Path) -> None:
    (dest.parent / _WORKING_COPY_MARKER).write_text(
        "The collection/ directory here is Tessera's working copy; it is cleared and "
        "rebuilt at the start of every run.\n"
    )


def new_working_collection(dest: Path) -> None:
    """Start an empty working collection at ``dest``, clearing a previous run's."""
    dest = Path(dest)
    _clear_working_copy(dest)
    dest.mkdir(parents=True)
    _mark_working_copy(dest)


def copy_collection(source: Path, dest: Path) -> None:
    """Replace ``dest`` with a fresh copy of the collection at ``source``.

    The working copy lives at ``<output>/collection`` and is rebuilt on every run so
    the user's input is never modified. That guarantee inverts when the input *is*
    that directory -- the natural way to continue from a previous run's output -- as
    clearing the destination first would delete the collection before copying it.
    Refuse, rather than destroy the input.
    """
    source, dest = Path(source), Path(dest)
    src, dst = source.resolve(), dest.resolve()
    if src == dst or dst in src.parents or src in dst.parents:
        raise UserInputError(
            f"The collection {source} is, contains, or lies inside this run's working "
            f"copy ({dest}), which is cleared at the start of every run. Choose a "
            "different output directory, or point --collection at a copy elsewhere."
        )
    _clear_working_copy(dest)
    shutil.copytree(source, dest)
    _mark_working_copy(dest)


def _require_fasta(source: Path) -> None:
    """Reject a staged input whose first non-blank character is not ``>``."""
    try:
        if source.name.endswith(".gz"):
            with gzip.open(source, "rb") as fo:
                head = fo.read(4096)
        else:
            with open(source, "rb") as fo:
                head = fo.read(4096)
    except (OSError, EOFError) as exc:
        raise UserInputError(f"Cannot read {source} as a (gzip) FASTA file: {exc}") from exc
    if not head.lstrip().startswith(b">"):
        raise UserInputError(
            f"{source} does not look like a FASTA file (it does not start with '>'). "
            "Every file in the collection directory is treated as a genome; move "
            "other files out of it."
        )


_WHITESPACE = re.compile(rb"\s")


def _has_sequence_whitespace(source: Path) -> bool:
    """True when a plain FASTA carries whitespace inside a sequence line.

    A trailing space or tab, a space within the line, or a carriage return (CRLF line
    endings). Blank lines and the header's own spaces do not count.
    """
    with open(source, "rb") as fo:
        for line in fo:
            if line.startswith(b">"):
                continue
            body = line[:-1] if line.endswith(b"\n") else line
            if body and _WHITESPACE.search(body):
                return True
    return False


def _write_clean(src: Iterable[bytes], dst: BinaryIO) -> None:
    """Copy a FASTA, dropping whitespace from sequence lines (and blank lines)."""
    for line in src:
        if line.startswith(b">"):
            dst.write(line.rstrip(b"\r\n") + b"\n")
        else:
            body = b"".join(line.split())
            if body:
                dst.write(body + b"\n")


def _stage_one(source: Path, target_dir: Path, logger: logging.Logger) -> Path:
    """Place one genome into ``target_dir`` as ``<label>.fasta``; return the path.

    The staged file is what the aligner reads, while Tessera reads the same genome
    through :func:`read_fasta`, which ignores whitespace in sequence lines. The two
    must agree on the genome's length, and not every aligner ignores a trailing space
    (minimap2 counts it as a base), so a genome that carries such whitespace is staged
    as a cleaned copy rather than a link to the original.
    """
    label = strip_sequence_extension(source.name)
    target = target_dir / f"{label}.fasta"
    if source.name.endswith(".gz"):
        logger.debug("Decompressing %s -> %s", source, target)
        with gzip.open(source, "rb") as src, open(target, "wb") as dst:
            _write_clean(src, dst)
    elif _has_sequence_whitespace(source):
        logger.warning(
            "%s has whitespace inside its sequence lines (trailing spaces, tabs or "
            "Windows line endings); aligning a cleaned copy.", source,
        )
        with open(source, "rb") as src, open(target, "wb") as dst:
            _write_clean(src, dst)
    else:
        logger.debug("Linking %s -> %s", source, target)
        target.symlink_to(source.resolve())
    return target


def stage_genomes(
    query: Path,
    collection_dir: Path,
    target_dir: Path,
    logger: logging.Logger,
) -> tuple[list[Path], Path]:
    """Stage the query and every collection file into ``target_dir``.

    Returns ``(all_genome_paths, query_path)``. Raises :class:`UserInputError`
    when the query or collection is missing or the collection is empty.
    """
    query = Path(query)
    collection_dir = Path(collection_dir)
    if not query.exists():
        raise UserInputError(f"Query file not found: {query}")
    if not collection_dir.is_dir():
        raise UserInputError(f"Collection directory not found: {collection_dir}")

    collection_files = collection_genomes(collection_dir)
    if not collection_files:
        raise UserInputError(f"Collection directory is empty: {collection_dir}")
    for source in [*collection_files, query]:
        _require_fasta(source)

    _reject_colliding_labels(collection_files, query)

    target_dir.mkdir(parents=True, exist_ok=True)
    staged: list[Path] = []
    for source in collection_files:
        staged.append(_stage_one(source, target_dir, logger))
    query_staged = _stage_one(query, target_dir, logger)
    staged.append(query_staged)
    return staged, query_staged


def _reject_colliding_labels(collection_files: list[Path], query: Path) -> None:
    """Refuse two inputs that would stage to the same label.

    Genomes are staged as ``<label>.fasta`` with the label taken from the filename, so
    ``X.fa`` and ``X.fasta`` -- or a query whose stem matches a collection member --
    collide. Staging would then fail on the symlink with a bare ``FileExistsError``,
    reported to the user as "Unexpected error". The label is also the sequence's name
    in the alignment, so a collision is genuinely ambiguous rather than something to
    resolve silently.
    """
    seen: dict[str, Path] = {}
    for source in [*collection_files, Path(query)]:
        label = strip_sequence_extension(source.name)
        if label in seen:
            raise UserInputError(
                f"Two inputs share the label '{label}': {seen[label]} and {source}. "
                "Genomes are identified in the alignment by filename without its "
                "extension, so these cannot be told apart -- rename one."
            )
        seen[label] = source


def select_reference(
    genomes: list[Path],
    query: Path,
    query_as_backbone: bool,
    explicit_reference: str | None = None,
) -> Path:
    """Choose the backbone genome for a reference-anchored alignment.

    Preference: an explicitly named reference (by label or filename), else the
    query when ``query_as_backbone`` is set, else the first non-query genome.
    Raises :class:`UserInputError` when an explicit reference cannot be matched.
    """
    if explicit_reference is not None:
        wanted = strip_sequence_extension(Path(explicit_reference).name)
        for genome in genomes:
            if genome.stem == wanted:
                return genome
        raise UserInputError(
            f"Reference '{explicit_reference}' not found among the staged genomes."
        )

    if query_as_backbone:
        return query

    for genome in genomes:
        if genome != query:
            return genome
    raise UserInputError("No reference candidate found (collection has only the query).")
