# Elliptic rank 30: certificate-first research pipeline

This repository studies the open constructive target

\[
\operatorname{rank} E(\mathbb Q) \ge 30
\]

by optimizing for a final proof certificate, not for a large analytic-rank score.

## Current truth status

**new search method**

No curve with 30 certified independent rational points has been found here.
The package does independently close the baseline reproduction gate for the
current public rank-29 record and preserves all evidence needed to audit that claim.

## Baseline result

For

```text
y^2 + x*y = x^3 - A*x + B
A = 27006183241630922218434652145297453784768054621836357954737385
B = 55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497
```

this repository verifies:

- all 29 complete rational points in `points.json` exactly;
- nonzero discriminant and a global minimal model;
- trivial rational torsion;
- an exact rank-29 local mod-2 independence certificate;
- canonical-height matrices in SageMath 10.9 and Magma V2.29-9;
- nonzero regulator and positive-definite numerical height matrix;
- exact conductor, local reduction data, and global root number `-1`;
- p-saturation of the listed subgroup for every prime below `4096`, independently in Sage and Magma.

Therefore, unconditionally,

\[
\operatorname{rank} E(\mathbb Q) \ge 29.
\]

The package does **not** promote the published conditional upper bound to an
unconditional exact-rank statement.

## Repository layout

```text
curve.json, points.json    exact baseline inputs
candidate_record.json      canonical candidate record
certificates/              compact machine-readable proof summaries
evidence/                  raw Sage/Magma/CI evidence and preserved failed runs
paper/paper.tex            single continuously updated research paper source
paper/paper.pdf            compiled paper
src/                       candidate, audit, sieve, and workstream utilities
search/                    explicit-family workers and missing-data request
tests/                     package integrity tests
verify_exact.py            dependency-free exact lower-bound verifier
verify_sage.py             independent Sage verifier
verify_magma.m             independent Magma verifier
```

The paper is updated in place. Do not create versioned `paper_v2`, `paper_v3`,
etc. Historical states belong in Git history.

## Reproduce the exact certificate

```bash
python3 verify_exact.py
python3 rank_packet_obstruction_tests.py
python3 -m unittest discover -s tests -v
sha256sum -c MANIFEST.sha256
```

## Reproduce in SageMath

```bash
sage -python verify_sage.py
sage -python verify_sage.py --saturate
```

The second command directly tests p-saturation for all primes below 4096 and is
substantially more expensive.

## Reproduce in Magma

```bash
magma verify_magma.m
```

The raw computations used for this package are archived under `evidence/` so
verification does not depend on GitHub Actions retention.

## A long and useful failed route

Treating 13 extra points over a rank-17 fibration as 13 independent quadratic
splitting characters forces a multiquadratic auxiliary curve of genus at least
`20481`. Low genus is lost already at the fourth independent character. This is
a restricted obstruction to that construction mechanism, not a universal rank bound.

## What the failure reveals

One character need not carry only one point. The package verifies an explicit
family in which one nontrivial quadratic character carries three independent
Mordell-Weil directions. Product-twist channels can also contribute sections
that were not used to construct the original covers.

## The decisive change of setting

Search over complete **Galois-character height packets**. For every low-genus
branch code, compute all twist Mordell-Weil lattices, same-character rank
multiplicities, product twists, successive minima, and global-solubility data.
Only then sieve specializations.

## Current search obstruction

The rank-29 announcement identifies the successful parameter of the second
rank-17 K3 fibration, but the public primary material examined does not provide
a complete executable package containing the exact `Q(t)` model, all 17 labelled
generic sections, the bad-parameter locus, and the coordinate map to the record
fiber. `search/REQUEST_FOR_FIBRATION_DATA.md` states the precise missing input. A
dated public-source audit is in `search/source_search_2026-08-06.md`, and
`search/ingest_family.py` provides a fail-closed exact-family ingestion gate.
No blind coefficient search is substituted.

## Certificate policy

A future rank-30 promotion requires a specified Weierstrass equation, 30 full
rational points, exact substitutions, torsion, a canonical-height matrix,
positive determinant, an independent exact independence certificate, saturation,
software versions, provenance, hashes, and verification in a second system.

See `STATUS.json`, `candidate_record.json`, and `paper/paper.pdf`.
