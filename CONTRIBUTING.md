# Contributing

This repository is a certificate-first mathematical research project.  The
priority is truth status, not apparent importance.

## Claims

Every contribution must separate:

- exact algebraic statements;
- proved computational statements with software/proof flags;
- conditional statements;
- finite-field discovery evidence;
- heuristics;
- failed routes.

A curve is not a rank-30 candidate because of a high Mestre--Nagao score,
Selmer dimension, analytic-looking data, or 30 unverified points.

## Reproducibility

Use exact integers and rationals.  Preserve complete inputs, software versions,
commands, hashes, transformation histories, and raw evidence.  Verification
scripts must work from the equation and points alone and must not depend on the
search database.

## Paper

There is one living paper:

```text
paper/paper.tex
paper/paper.pdf
```

Modify it in place; Git history is the version archive.  The author line is
`Francesco Giannicola`.

## Final promotion

A certified result requires an explicit nonsingular elliptic curve over `Q`,
30 complete rational points, exact equation checks, rigorous independence,
canonical-height evidence, necessary saturation analysis, an adversarial
audit, and independent verification in a second computer algebra system.
