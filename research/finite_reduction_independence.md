# Finite-reduction independence certificates

For an elliptic curve `E/Q` with good reduction at odd primes `p_1,...,p_s`, reduction followed by quotienting by doubles gives a homomorphism

\[
E(\mathbf Q)/2E(\mathbf Q) \longrightarrow \prod_j E(\mathbf F_{p_j})/2E(\mathbf F_{p_j}).
\]

If the images of rational points `P_1,...,P_n` have row rank `n` over `F_2`, their classes in `E(Q)/2E(Q)` are independent. If `E(Q)_tors=0`, every integral relation has even coefficients; division by two and infinite descent prove that the points are `Z`-independent.

For good reduction, rational torsion injects into the finite groups away from the residue characteristic. A collection of exact group orders with gcd one therefore proves trivial rational torsion.

The implementation in `research/finite_reduction_independence.py` is deliberately independent of the cubic-field Kummer certificate. It uses exact rational arithmetic, direct finite-group enumeration, canonical coset labels, binary Gaussian elimination, stored torsion witnesses, and a canonical SHA-256 hash. It does not use analytic rank, BSD, or floating-point heights.

Reproduction:

```bash
python3 research/finite_reduction_independence.py \
  --verify baseline/rank29_finite_reduction_certificate.json
```

This is a lower-bound certificate, not a proof that the full Mordell–Weil rank is exactly 29 and not a saturation certificate.
