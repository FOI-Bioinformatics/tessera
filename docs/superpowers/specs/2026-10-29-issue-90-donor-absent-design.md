# Issue #90: naming a donor when the donor's lineage is absent

Date: 2026-10-29. A design note, not a plan to build. No code is proposed for merge
here, and nothing in it has been measured on the harness: the one measurement it asks
for needs network access the audit environment does not have. Written because the issue
is the last one open, its premise has already been corrected once, and four candidate
fixes have already been measured and eliminated on it. The purpose is to record which
kinds of evidence remain, what would falsify each, and what the decision actually costs
-- so the next attempt does not re-run the four that failed.

## The failure, as isolated on the issue

With the donor's whole sub-clade (`A.D.1`) removed from an RSV-A panel, the scan reports:

| | |
|---|---|
| true donor lineage | `A.D.1`, entirely absent |
| named donor | `EPI_ISL_15120672`, clade `A.D.5` |
| identity over the tract | 0.9851 |
| callers agreeing | four (hmm, 3seq, maxchi, bootscan) |
| `donor_absent` / `donor_undercovered` | no / no |

So a donor from the wrong sub-clade is named with four-method agreement and no caveat.
The detection is right and the attribution is wrong, and nothing in the output says so.

This matters more than its one harness case: an unsampled true parent is the ordinary
situation in real use, and the honest output is "something recombined here and I cannot
name the donor".

## What is already ruled out, and why

From the issue's own measurements across 37 scorable hybrid cases. Repeating any of
these is wasted work.

| # | signal | why it fails |
|---|---|---|
| 1 | absolute identity threshold | stand-in sits at 0.9851-0.9854 against a 0.950 bar; no absolute bar catches it without catching most real donors |
| 2 | donor identity relative to the backbone's own flank identity | donor-present deficits span -0.0152 to +0.0367; this case sits at +0.0107 to +0.0131, inside the range, with seven true positives equal or worse |
| 3 | donor-lineage agreement among overlapping regions | this case names 2 distinct sub-clades, the same count as `rsv_a` with the donor present |
| 4 | the tract's best match against the panel's within-sub-clade identity scale | the scales overlap heavily, and only n=3 pairs are available |

**They share one cause.** RSV's `A.D.*` sub-clades are nearly equidistant from the
tract, so no statistic computed from pairwise identities distinguishes "the donor's own
lineage" from "a neighbouring one". The stand-in's identity, 0.9851, is within 0.0003 of
the same-sub-clade stand-in's 0.9854 from the original premise. Identity carries
essentially no information here, and every signal above is a function of identity.

The conclusion to carry forward: **a fix must use evidence of a different kind, not
another threshold on identity.** That is the issue's own conclusion and this note does
not improve on it.

## The three kinds of evidence that remain

### Option A: placement, not best match

The quantity that separates "the donor is here" from "something near the donor is here"
is not how well the best reference matches, but *how the tract attaches to the panel's
structure*. A tract from a panel member attaches by a short pendant branch. A tract from
an unsampled lineage attaches by a long one, because the sites that distinguish its true
parent from every panel member are all on that branch.

Concretely, and within the dependency-free constraint (neighbour joining is a few lines
of numpy; no tree library, no UShER, no external binary):

1. restrict the alignment to the region's columns;
2. build a distance matrix over the panel plus the query's tract;
3. neighbour-join it, and read the length of the pendant branch by which the tract
   attaches;
4. compare that length against the pendant branches of the panel's *own* members on the
   same columns -- the panel calibrates itself, which is what makes the statistic
   scale-free across pathogens and divergences.

**The test statistic** is the tract's pendant length as a quantile of the panel's own
pendant lengths over the same columns. Donor present: the tract should sit inside the
panel's distribution. Donor absent: above it.

**Why this is not discriminator 2 in another form.** A deficit in identity is a function
of the single best match. A pendant length is a function of the sites the whole panel
fails to explain, jointly. The two are correlated but not the same, and the audit's
measurements constrain only the first.

**What would falsify it**, and this must be measured before any of it is built: compute
the statistic on all 37 scorable hybrid cases with the donor present, and on
`donorabsent_rsv`. If the donor-absent value falls inside the donor-present
distribution, as the identity deficits did, Option A is dead and should be recorded here
as a fifth eliminated discriminator. **The prior is not favourable**: these sub-clades
are nearly equidistant, which is exactly the geometry that makes pendant lengths similar
too. This is a candidate to test, not a fix to schedule.

Cost if it works: one distance matrix and one NJ per called region, over a region's
columns only. Acceptable. Cost if it is built without the measurement: a threshold
calibrated on a single case, which is the thing the issue explicitly declines to do.

### Option B: say the resolution is insufficient, rather than that the donor is absent

`donor_absent` and `donor_undercovered` are binary claims about the panel. Neither is
true here: a close relative *is* present, and the named lineage is wrong by one
sub-clade, not absent. A third, weaker statement is available and is true:

> the donor was resolved to `A.D.5`, but `A.D.3` and `A.D.5` are within the panel's own
> noise of each other over this region, so the sub-clade is not identifiable.

This needs no new evidence. It is computed from quantities already on the region: the
margin between the best donor and the runner-up from a *different* lineage, against the
panel's within-lineage spread over the same columns. Where the between-lineage margin
does not exceed the within-lineage spread, the sub-clade call is not resolved and the
report should say which lineages are tied rather than name one.

**This does not fix the issue's scored case.** The harness asks for `donor_absent` or
`donor_undercovered`, and this would set neither. It is a smaller, honest improvement to
what a reader sees, and it would stop the *confident* presentation that the issue
objects to -- four callers agreeing on a name that is not identifiable. It is also the
only option here that is clearly implementable today, and discriminator 4's measurement
("the scales overlap heavily") is weak evidence *for* it: overlapping scales are
precisely the condition it would report.

### Option C: accept it as a limit of an identity-based screen

The project has a `known_limit` mechanism, reported but not gating, used for `neg_hiv1`.
The case qualifies on the stated grounds: in a panel whose sub-clades are equidistant
from the tract, the information needed may genuinely not be present in the alignment.

The README already says Tessera is an indicative screen rather than a full phylogenetic
test, and that strong candidates should be confirmed with a dedicated method. A limit
that follows directly from that framing is defensible to document.

This is **the maintainer's call and not the agent's**, for the reason the issue gives:
marking one's own unsolved problem as a known limitation is a convenient move. It is
recorded here as a real option, not recommended from here.

## Recommendation

Ordered, and the first item is the only one that should start now.

1. **Measure Option A before building any of it.** One script over the hybrid cases,
   reporting the tract's pendant-length quantile for every scorable case plus
   `donorabsent_rsv`. Publish the distribution in `validation/` next to the attribution
   results, whichever way it comes out. This is the project's own evaluate-first
   precedent from `docs/gap-analysis.md` G1, where the barcode caller was measured
   before promotion and found inert.
2. **Option B is worth doing regardless of Option A's outcome**, because it is true,
   cheap, and addresses the part of the complaint that is about presentation: a name
   offered with four-caller agreement when the sub-clade is not identifiable. It should
   be its own change with its own measurement, not folded into Option A.
3. **Option C, if Option A's measurement comes out as the prior suggests.** That is a
   decision about what the tool claims, and it belongs to the maintainer.

## What blocks item 1 today

`run_hybrids.py` needs `data.clades.nextstrain.org` and the NCBI hosts to build its
pools, and the audit environment's network policy denies both. The measurement needs a
machine with those hosts allowed, and a warm Nextclade cache. Nothing in this note
should be built before it runs.
