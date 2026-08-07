# Incremental independence

For every new point `P`:

1. verify `P` exactly on the curve;
2. reject duplicates, negatives, torsion, and known combinations;
3. test modular dependence;
4. reduce against the known Mordell--Weil lattice;
5. update the height Gram matrix and exact local certificate;
6. retain every coordinate transformation and provenance record.

The full basis is not recomputed from scratch after every point.
