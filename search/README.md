# Search workstreams

Every program in this directory is either an exact construction/certification
step or clearly labelled discovery evidence.  A large analytic score is never
a candidate level by itself.

## Active exact tracks

- `sage_solve_j0_e8_sections.py` — complete minimal-section scheme for marked
  cyclic `j=0` rational elliptic surfaces.
- `verify_auxiliary_cm_curve.sage.py` — proved Sage rank/torsion computation for
  the auxiliary CM descent curve.
- `sage_finite_field_function_rank.py` — finite-field function-field rank upper
  bound experiment.
- `sage_search_split_e8_pencils.py` — general-position cubic pencils; promotion
  requires a square-free degree-12 discriminant and split rank eight.
- `sage_scan_split_e8_packets_finite_field.py` — exact characteristic-p packet
  ranks for all three branch channels; discovery evidence only until lifted.
- `sage_norm31_seed_search.py` — exact finite-field seeds for the 31-zero norm
  equation, plus explicitly heuristic centered integer lifts.

## Promotion discipline

No file in `search/results/` proves a rank statement over `Q` unless its
machine-readable record includes explicit characteristic-zero equations and
sections and passes the independent verification package.  Failed and partial
runs remain useful for bottleneck analysis but are not silently promoted.
