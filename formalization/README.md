# Lean formalization boundary

The repository contains a Lean 4 / mathlib formalization of the exact
algebraic and finite arithmetic layer used by the research certificates.

Pinned versions:

```text
Lean 4.32.2
mathlib 4.32.2
```

Build it with:

```bash
lake update
lake exe cache get
lake build
```

The imported root is `EllipticRank30.lean`.

## Formally checked in Lean

- the birational identity taking `u^3 = c*t*(t-1)` to
  `Y^2 = X^3 + 16*c^4`;
- positivity of `mu^2-mu+1` and the resulting real obstruction to the
  minimal constant equation;
- the repeated-factor and quadratic-twist identities at `mu=2`;
- the generic marked coefficient factorization and the two symmetric root
  identities;
- both normalized solutions used by the minimal-scale certificate;
- all numerical Euler, root-rank, rank-capacity, trace-code, and shell-count
  arithmetic appearing in the six-channel reduction.

## Not yet formalized

Lean does not currently certify the geometric inputs that turn those
identities into Mordell--Weil theorems. In particular, the repository still
uses separately checked mathematical arguments and CAS certificates for:

- Kodaira fibre classification and minimal-model analysis;
- Shioda--Tate and Hodge/Picard bounds;
- the quadratic and cubic character decomposition of Mordell--Weil groups;
- irreducibility and field-of-definition claims not represented in the Lean
  files;
- finite-field point orders, function-field ranks, canonical heights,
  saturation, and rational-point independence.

The paper marks this boundary explicitly. A green Lean workflow means that
all declarations in `EllipticRank30` compile without `sorry` or `admit`; it is
not by itself a proof of a rank-30 curve.
