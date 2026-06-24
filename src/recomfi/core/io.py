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
import shutil
from pathlib import Path

from .errors import UserInputError

# Extensions stripped to derive a genome's label (its leaf name in the MSA).
SEQUENCE_EXTENSIONS = (".fna", ".fasta", ".fa")


def strip_sequence_extension(name: str) -> str:
    """Return ``name`` without a trailing ``.gz`` and/or sequence extension."""
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    for ext in SEQUENCE_EXTENSIONS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _stage_one(source: Path, target_dir: Path, logger: logging.Logger) -> Path:
    """Place one genome into ``target_dir`` as ``<label>.fasta``; return the path."""
    label = strip_sequence_extension(source.name)
    target = target_dir / f"{label}.fasta"
    if source.name.endswith(".gz"):
        logger.debug("Decompressing %s -> %s", source, target)
        with gzip.open(source, "rb") as src, open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
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

    collection_files = sorted(p for p in collection_dir.iterdir() if p.is_file())
    if not collection_files:
        raise UserInputError(f"Collection directory is empty: {collection_dir}")

    target_dir.mkdir(parents=True, exist_ok=True)
    staged: list[Path] = []
    for source in collection_files:
        staged.append(_stage_one(source, target_dir, logger))
    query_staged = _stage_one(query, target_dir, logger)
    staged.append(query_staged)
    return staged, query_staged


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
