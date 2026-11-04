# Fourth audit: implementation plan index

**Spec:** `docs/superpowers/specs/2026-10-12-audit-swot-gap.md`. **Written:** 2026-10-12, against
commit `4b35e39` on branch `fix-audit-caller-behaviour`.

The audit spans five subsystems, so it is implemented as four plans, each shippable on its
own, split by risk in the way plans A, B and C were. Plans D, E and F do not change which
regions the scan calls; plan G does, and every behaviour-changing task in it is gated on the
harnesses.

| Plan | File | Tasks | Changes region calls? | Branch |
|---|---|---|---|---|
| D -- alignment integrity and adapter resources | `2026-10-12-audit-d-alignment-integrity.md` | 12 | no; changes the alignment that reaches the scan | `fix-audit4-alignment` |
| E -- panel discovery loop | `2026-10-12-audit-e-discovery-loop.md` | 11 | no; changes which genomes a recruited panel holds | `fix-audit4-discovery` |
| F -- reporting faithfulness, CLI and reassortment | `2026-10-12-audit-f-reporting-cli.md` | 18 | no; checked by an identical specificity run before and after | `fix-audit4-reporting` |
| G -- caller statistics (harness-gated) | `2026-10-12-audit-g-caller-statistics.md` | 17 | yes; every such task carries a gate and a revert rule | `fix-audit4-callers` |

Every task's tests were run while the plans were written, in a scratch worktree of
`4b35e39`: each new test fails on the unchanged code for the reason its Step 2 records, and
passes with the candidate implementation. With each plan's candidates applied on its own,
the fast suite stays green, with ruff and mypy clean:

| Plan | Fast suite with candidates | Other checks run while planning |
|---|---|---|
| D | 794 passed (757 + 37) | `requires_binary` 8 passed, 1 skipped (mafft, minimap2, sibeliaz installed; progressiveMauve absent); specificity output identical |
| E | 793 passed (757 + 36) | none (no network, no harness reaches the changed branches) |
| F | 844 passed (757 + 87) | `run_specificity.py --reps 3 --min-methods 1` identical before and after |
| G | 796 passed (757 + 39) before Task 16 was added | specificity at both gates and the regime harnesses at 3 replicates; redundant-lineage regimes at 10 |

The four candidate sets were not applied together. Their combined test count is not
measured, and the overlaps listed under "Landing order" are where a merge needs care.

## Before starting

1. **Merge `fix-audit-caller-behaviour` first, or branch from it.** All four plans are
   written against its tip. It is not yet pushed; its final gate (`run_hybrids.py` at
   `--min-methods 1`) is recorded in its own plan.
2. **Plan C Task 6 is held** in `stash@{0}` and in
   `validation/data/audit-c/task6/task6.patch`, awaiting the maintainer's decision. It edits
   the same lines of `_type_panel` and `_select_from` as plan E Task 5. Plan E Task 5 states
   how to apply it in either order. Do not pop or drop that stash; when moving work between
   worktrees, use `git diff > file` and `git apply`, never `git stash`, because the stash is
   shared by every worktree of the repository.
3. **Append the aligner environment to `PATH`, never prepend it.** Prepending replaces the
   Python interpreter (plan F Task 17 corrects the documentation that says otherwise).

## Landing order

**D, then E, then F, then G.** The reasons, in order of weight:

- **G last.** It is the only plan that changes region calls, its gates compare against the
  state the other three leave, and its Task 16 consumes plan F's `not_run` mechanism.
- **F before G.** Plan F Task 5 and plan G Tasks 4, 13 and 16 edit the same help strings,
  `RecombParams` comments and documentation paragraph, and G changes what they must say.
  Both sides carry an "if the other landed first" note with the replacement text.
- **D first.** It changes what alignment reaches the scan. Landing it before G means G's
  gates are measured on corrected alignments.
- **E before F.** They share three command files and `docs/reference-panels.md` (table
  below). The edits do not conflict in meaning; landing E first means F's validation and help
  changes are written against E's seed-source resolution, which moves into `FillParams`.
  E is independent of G. Its one cross-plan interaction outside F is plan C Task 6, above.

Files more than one plan edits, beyond `CHANGELOG.md`:

| File | Plans | Note |
|---|---|---|
| `src/tessera/recomb/run.py` | D (Task 6), F (Tasks 2-8), G (Tasks 5, 11, 13, 16) | edits are anchored on quoted text; rebase on whichever lands first |
| `src/tessera/cli/cmd_recomb.py` | F (Task 5), G (Tasks 13, 16) | see the F-before-G overlap notes above |
| `docs/detection-methods.md` | F (Tasks 2-6, 8), G (Tasks 4, 5, 10-13, 15, 16) | as above |
| `src/tessera/cli/cmd_fill_references.py` | E (Task 1), F (Tasks 10-12) | E1 range-checks `--min-improvement` and F11 the coverage fractions and counts, both after the same validation line; keep both checks |
| `src/tessera/cli/cmd_build_panel.py` | E (Task 11), F (Tasks 11, 12) | E11 removes `_seed_source` and moves seed-source resolution into `FillParams`; F12 adds `--pool-consensus` to `build-panel`, which E11's check then governs. F12's test passes `--nextclade` so it holds in either order |
| `src/tessera/cli/cmd_detect.py` | E (Task 11), F (Task 11) | different lines |
| `docs/reference-panels.md` | E (Tasks 1, 2, 5-8, 11), F (Tasks 12-15) | different sections |

`src/tessera/core/io.py`, `src/tessera/discover/nextclade.py` and `docs/aligners.md` are
each edited by plan D alone.

## Gates that were not run while planning

These are executor gates. A plan's pull request states which ran and what they showed; a
gate that could not be run is written as not run, with the reason.

- **`validation/run_hybrids.py`** (network): plan G, before and after, at both
  `--min-methods` settings, as `CONTRIBUTING.md` requires for any caller change.
- **`validation/run_validation.py`** (aligners and fetched data): plan G, and plan D as a
  check that the alignments it changes leave the curated verdicts unchanged.
- **One paired real-data `detect` run** before plan E Task 2 merges, since no offline
  harness reaches the download or curation path.
- **The orthopox example with `--aligner sibeliaz` on macOS** after plan D Task 3, to test
  whether the "sibeliaz runs out of memory on orthopox" note is the BSD wrapper
  incompatibility rather than memory.

## Traceability: audit item to task

### Callers (spec section 1.1)

| Item | Plan and task | Note |
|---|---|---|
| A1 cluster-size weighting | G10 | mean within a cluster; gated on G2's redundant-backbone regime |
| A2 redundant donor splits an event | G11 | shared clusters for merge, Bootscan and HMM; a residual case under informative-site windowing is recorded, not fixed |
| A3 descent floor | G12 | 3SEQ tests every donor with a faster exact DP; MaxChi keeps its floor for runtime, counting shallow donors at p = 1 |
| A4 `--min-region` vs trimmed span | F5 (help), G13 (behaviour) | |
| A5 Bootscan permutation | G6 | the spec's example numbers did not reproduce; the defect is confirmed another way |
| A6 BH family pruned | G12 | measured: did not change a call on the simulated panels |
| A7 MaxChi int64 overflow (suspected) | G7 | widening with a guard test |
| A8 merged `n_windows` units | G5 | |
| A9 sign test underflows to zero | G3 | log-space tail |
| A10 HMM half-width | G4 | median step; the literal half-window reading changes significance and goes to the decision record |
| A11 merged identity | G5 | |
| A12 3SEQ and MaxChi not independent | G15 | documentation correction and decision-record addendum |
| A13 memory sites | G8 | |
| A14 barcode vote | G15 (record only) | maintainer decision |
| A14 terminal breakpoint interval | G9 | also covers `maxchi.py:170`, which the spec missed |
| trailing unscanned columns | G14 | |
| incompatible-pair list | G8 | |

### Orchestration and reporting (spec section 1.2)

| Item | Plan and task |
|---|---|
| B1 `--no-phi` reported as untestable | F2 (also covers too few sites and not testable) |
| B2 `--min-region` scope | F5 (help), G13 (behaviour) |
| B3 low-information gap headline | F4 |
| B4 `n_windows`, `posterior` after trimming | F5 (docs), G4 (code) |
| B5 provenance window geometry | F6 |
| B6 `--margin`, `--merge-gap` scope | F5 (help), G16 (`--merge-gap` reaches the coverage diagnostic) |
| B7 one-state HMM | F3 (reported as not run), G16 (no longer counts toward the gate) |
| B8 pair-plot colours | F9 |
| B9 headline backbone | F7 (confirmed while planning) |
| smaller reporting items | F4, F6, F8 |

### Core, aligners and converters (spec section 1.3)

| Item | Plan and task |
|---|---|
| C1 progressiveMauve backbone row | D2 |
| C2 duplicate sequence ID in one genome | D4 |
| C3 minimap2 thread budget | D5 |
| C4 Nextclade pool cache install | D8 |
| C5 unplaced-genome warning | D6 |
| C6 truncated MAF row | D7 |
| C7 non-directory at the collection path | D9 |
| standard input, log handlers | D10 |
| cactus sample names, XMFA flank path | D11 |
| sibeliaz on macOS (found while planning) | D3, separable |

### Panel discovery (spec section 1.4)

| Item | Plan and task |
|---|---|
| D1 download fan-out | E2 (with excluding curated-away genomes; approve together) |
| D2 stop test | E1 (coverage shortfall; `--min-improvement` becomes a fraction) |
| D3 lineage sidecar | E5 |
| D4 false saturation | E3 |
| D5 zero RefSeq | E4 |
| D6 curation never ran | E6 |
| D7 all-sibling pool | E7 |
| D8 `--fetch-limit` cache key | E9 |
| D9 short-query windows | E8 |
| D10 non-clade markers | E10 |

### CLI and reassortment (spec section 1.5)

| Item | Plan and task |
|---|---|
| E1 output path is a file | F10 (seven commands, not five) |
| E2 `reassort --margin` bounds | F13 |
| E3 reassortant with one group | F14 (verdict becomes undetermined) |
| E4 pairwise evidence discarded | F14 (`segment_pairs.tsv`) |
| E5 seed-source help, `find-references -q`, `msa` conflict and overwrite, `build-panel` options | F12 |
| E5 coverage fractions, unbounded counts | F11 |
| E5 `--ani-floor` fraction, two `--margin` units | F13 |
| E5 `reassort --threads`, clade tie | F15 |
| E5 `--candidate-pool` vs `--nextclade`, inert `--pool-consensus`, pool reader extensions | E11 |
| E14 empty candidate list | F14 |

### Repository and documentation (spec section 1.6)

| Item | Plan and task |
|---|---|
| F1 harness argument parsing | F16 (also `run_regimes.py`, `characterise_false_regions.py`) |
| F2 `PATH` lines | F17 (also `CONTRIBUTING.md`; three lines in `validation/README.md`) |
| F3 mypy text, F4 lint line, F5 README false-region rate | F17 |

## Not planned, and why

**Maintainer decisions.** These change what an option or a gate means, and the plans record
the options and evidence instead of choosing:

- The `--min-methods` default, still open in
  `docs/superpowers/specs/2026-10-01-min-methods-default-decision.md`. Plan G Task 15 adds
  the A12 dependence numbers and the A14 options to it, and plan G records that the literal
  half-window reading of A10 cut false regions from 7/12 to 1/12 at the default gate while
  the HMM lost its near-identical calls (3/3 to 0/3); that belongs with option B.
- Whether the barcode caller, which has no significance test, counts toward `--min-methods`.
- Whether `--metric`, which has one legal value, is removed.
- Whether genomes in subdirectories of a collection are read, warned about or refused.
- Whether seed options are refused when `--collection` is given (they are inert there).
- Whether to curate a freshly seeded panel before round 1 (an existing test pins the current
  behaviour as deliberate; plan E Task 6 reports instead).

**Suspected items not reproduced.** The negative-lineage seeder searches the first 1500
bases, not the whole query, as its docstring claims; and reference typing resolves its
Nextclade dataset from whichever untyped genome sorts first. Both are confirmed in code, with
no consequence reproduced.

**Found while planning, outside the spec, not planned.**

- A donor-absent region's `parent_free_support` is never computed, so it always reads "no".
  This is the not-run-as-negative pattern again; fixing it changes a reported value.
- Under informative-site windowing, a donor copy 0.03 % away is not pooled by lineage
  clustering, so that redundant-donor event is still reported as two rows. Closing it means
  choosing plan C's informative-column window or floor, gated on `run_hybrids.py`.
- GENECONV has the same pre-test floor and family pruning as 3SEQ and MaxChi.
- `validation/run_benchmark.py`, `run_coalescent_benchmark.py` and `run_caller_benchmark.py`
  hand-parse their arguments, as `run_specificity.py` did.
- `reassort/scan.py:154` discards a computed value that nothing uses; an existing test pins
  the three-value return.

## Corrections to the spec found while planning

The spec is an audit record and is not rewritten; these corrections are listed here and in
its closing section.

- **Aligner binaries.** The spec states that no aligner binary is installed. The `recomfi-aln`
  environment holds mafft, minimap2 and sibeliaz 1.2.7; only progressiveMauve and cactus are
  absent, and CI installs mafft alone. The audit passes ran without that environment on
  `PATH`.
- **The default backend fails on macOS.** sibeliaz 1.2.7's wrapper uses two GNU-only
  commands, so on macOS it wrote an empty alignment for every input tried, and the adapter
  blamed memory (plan D Task 3).
- **D1 magnitudes.** At the default `--max-hits 5`, a gap can add up to 20 genomes per round
  for 1 kb and 145 for 6 kb, not 60 and 435.
- **D9.** Hepatitis B already gets two windows; the affected range is 1.5 to 3 kb, for
  example hepatitis delta.
- **E1** affects seven commands; `find-references` and `curate-panel` were missed.
- **B1** also covers PHI runs with too few informative sites and those not testable at their
  window. **B9** reproduces, in the case with no called region.
- **C7** is narrower than described: the raw error needs a previous run's marker file, or a
  symlink at the path.
- **A3's remedy as written** makes a 16-genome scan about ten times slower; plan G takes a
  cheaper route.
- **A5's example numbers** do not reproduce as stated.
- **Line numbers** for B2 and B6 are 16 lines off; each plan records the verified lines.
