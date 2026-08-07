# Elliptic Rank 30

**Certificate-first research toward an elliptic curve over** \(\mathbf Q\) **with at least 30 independent rational points.**

Author and project lead: **Francesco Giannicola**

## Truth status

> **new intermediate theorem**

No rank-30 curve is claimed here yet. The repository independently certifies
the current rank-at-least-29 baseline and develops exact construction and
obstruction results aimed at a final Level-5 certificate.

The historical slowdown of rank records is not treated as evidence of a
universal upper bound. A universal claim would require a universal theorem.

## Reproduced baseline

For the published record model

```text
y^2 + x*y = x^3 - A*x + B
A = 27006183241630922218434652145297453784768054621836357954737385
B = 55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497
```

this repository verifies:

- all 29 complete rational points exactly;
- nonsingularity and a global minimal model;
- trivial rational torsion;
- an exact local mod-2 independence certificate of rank 29;
- canonical-height matrices in SageMath and Magma;
- nonzero regulator and positive-definite numerical height matrix;
- conductor, local reduction data, and global root number;
- prime saturation of the displayed subgroup for every prime below 4096.

Therefore, unconditionally,

\[
\operatorname{rank}E(\mathbf Q)\ge 29.
\]

The published conditional upper bound is not promoted to an unconditional
exact-rank statement.

## Current mathematical frontier

### 1. Galois-character packets

For a multiquadratic packet with character dimension \(k\) and total branch
support \(b\),

\[
k\le b-1,
\qquad
g=1+2^{k-2}(b-4).
\]

Thirteen independent quadratic characters force genus at least 20481.
However, one character can carry several Mordell--Weil directions, so searches
must compute complete twist packets, including hidden product-twist channels.

### 2. One-function visibility barrier

If one rational function has zeros \(P_1,\ldots,P_N\) and pole divisor \(NO\),
then

\[
P_1+\cdots+P_N=O.
\]

A split degree-30 norm identity therefore spans rank at most 29. The corrected
one-function target uses 31 zeros.

### 3. Certified split \(E_8\) surface

The repository certifies a rational elliptic surface with twelve \(I_1\)
fibres and split Mordell--Weil lattice \(E_8\). This is the maximum-capacity
rational base: rank 8 with no reducible-fibre cost.

### 4. Cyclic cubic trace codes

A degree-three **genus-one** pullback is the first degree-three low-genus
mechanism whose Hodge capacity reaches rank 30. For a cyclic trisection orbit:

- the three sections have projected Gram matrix
  \(\begin{psmallmatrix}4&-2\\-2&4\end{psmallmatrix}\), independent of the
  trisection's intersection with the zero section;
- the trace vectors define an exact ternary code in \(E_8/3E_8\);
- the norm-six shell has 6720 classes;
- there are 2240 maximal totally isotropic four-spaces;
- each contains 240 minimal trisection classes.

The decisive finite search is now: find **eleven Hermitian-independent cyclic
trisection orbits defining one positive-rank genus-one cover**. Together with
the invariant \(E_8(3)\), they would force generic rank 30.

See:

- [`research/cyclic_cubic_trace_codes.md`](research/cyclic_cubic_trace_codes.md)
- [`certificates/cyclic_trisection_trace_code.json`](certificates/cyclic_trisection_trace_code.json)
- [`paper/paper.pdf`](paper/paper.pdf)

## Repository layout

```text
curve.json, points.json          exact rank-29 baseline inputs
candidate_record.json            canonical candidate record
certificates/                    compact machine-readable certificates
evidence/                        preserved Sage/Magma/CI evidence
paper/paper.tex                  living paper entry point
paper/paper.pdf                  compiled paper
paper/source/                    modular LaTeX sources
research/                        theorems, proofs, and exact certificates
search/                          discovery workers and exact search outputs
src/                             candidate, sieve, and audit utilities
tests/                           integrity and theorem tests
verify_exact.py                  dependency-free baseline verifier
verify_sage.py                   independent SageMath verifier
verify_magma.m                   independent Magma verifier
```

## Reproduce the exact checks

```bash
python3 verify_exact.py
python3 research/cyclic_trisection_trace_code_certificate.py
python3 -m unittest discover -s tests -v
sha256sum -c MANIFEST.sha256
```

SageMath:

```bash
sage -python verify_sage.py
sage -python verify_sage.py --saturate
```

Magma:

```bash
magma verify_magma.m
```

## Candidate policy

Only a Level-5 package is a record. It must contain one exact curve, thirty
complete rational points, exact substitutions, torsion, a rigorous height
matrix, nonzero determinant, an independent exact independence proof,
saturation, a second-CAS rerun, hashes, versions, and provenance.

Analytic scores, Selmer dimensions, finite-characteristic ranks, and machine-
learning rankings are discovery evidence only.

## Contributing

Contributions are most useful when they provide one of:

- an explicit family with exact generic sections and a Shioda Gram matrix;
- a cyclic trisection or packet with exact branch and height data;
- a globally soluble covering producing a genuinely new quotient direction;
- an independent verifier or adversarial audit;
- executable data for the second rank-17 K3 fibration used in the rank-29
  construction.

Every mathematical claim must state whether it is proved, conditional,
computational evidence, or heuristic.
