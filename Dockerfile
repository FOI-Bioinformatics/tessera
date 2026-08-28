# Tessera with its aligner backends, so a result can be reproduced rather than only
# described.
#
# The Python package is pure and pip-installable; the aligners are not, which is the
# whole reason this image exists. `run_provenance.json` records which aligner and which
# version produced an alignment -- pinning them here is what makes that record
# actionable for a reader who was not there.
#
# Build:
#   docker build -t tessera:1.1.0 .
# Run (mount your data, write results back out):
#   docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" tessera:1.1.0 \
#     recomb --msa /data/panel.msa.fasta --query query --output /data/out
#
# `--user` is not optional on Linux. The image runs as a non-root user, and a bind
# mount keeps its ownership from the host, so without it the container cannot write
# its results into your directory. Passing your own uid/gid also means the outputs
# come back owned by you rather than by root. (Docker Desktop on macOS maps uids
# behind the scenes and hides this, which is exactly how it reaches Linux users
# unnoticed.) Everything the image needs at runtime is world-readable so any uid
# works.
#
# Excluded on purpose: Cactus (a pipeline in its own right, with Toil and its own
# containers) and progressiveMauve (no conda build for linux-aarch64). Both remain
# available through environment.yml on a host that can install them.

FROM mambaorg/micromamba:1.5.8-jammy

LABEL org.opencontainers.image.title="Tessera" \
      org.opencontainers.image.description="Detect recombination in a query genome against a reference panel" \
      org.opencontainers.image.source="https://github.com/FOI-Bioinformatics/tessera" \
      org.opencontainers.image.licenses="MIT"

USER root
WORKDIR /opt/tessera

# Pinned: an unpinned image is not a reproducible environment, which is the point of
# having one. Bump deliberately, and record the new versions in the changelog.
# All of these resolve on both linux-64 and linux-aarch64, so the image builds
# natively on x86 and Apple silicon alike.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.13 \
        sibeliaz=1.2.7 \
        mafft=7.526 \
        minimap2=2.31 \
        skani=0.3.2 \
        skder=1.3.7 \
        entrez-direct=26.0 \
        ncbi-datasets-cli=18.36.0 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

# Install the package itself from the source tree. Copy the metadata first so a code
# change does not invalidate the dependency layer.
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN pip install --no-cache-dir .

# Cache location inside the container; mount over it to persist recruited panels
# between runs (otherwise every run refetches).
ENV TESSERA_CACHE=/cache
# An arbitrary --user uid has no home directory, so anything that caches under $HOME
# either warns or falls back on every run: matplotlib wants MPLCONFIGDIR, and
# fontconfig (pulled in by the plotting stack) wants XDG_CACHE_HOME. Point both at a
# writable path so a normal run is quiet.
ENV MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/xdg-cache \
    HOME=/tmp
RUN mkdir -p /cache /data /tmp/matplotlib /tmp/xdg-cache \
    && chmod 777 /cache /data /tmp/matplotlib /tmp/xdg-cache \
    && chmod -R a+rX /opt/conda

# Identify the caller to NCBI. Supply your own at run time:
#   docker run -e NCBI_API_KEY=... -e NCBI_EMAIL=you@example.org ...
ENV NCBI_EMAIL=""

USER $MAMBA_USER
WORKDIR /data

ENTRYPOINT ["tessera"]
CMD ["--help"]
