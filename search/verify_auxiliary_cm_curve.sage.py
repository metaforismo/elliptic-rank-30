#!/usr/bin/env sage-python
from sage.all import EllipticCurve, proof
import json
from pathlib import Path

proof.all(True)
E = EllipticCurve([0, 0, 0, 0, -432])
T = E.torsion_subgroup()
rank = E.rank(proof=True)
P = E(12, 36)
assert 3 * P == E(0)
assert rank == 0
assert T.order() == 3
out = {
    "status": "pass",
    "curve": "y^2=x^3-432",
    "rank": int(rank),
    "torsion_order": int(T.order()),
    "torsion_invariants": [int(x) for x in T.invariants()],
    "point_12_36_order": int(P.order()),
    "proof_flags": "proof.all(True)",
}
Path("certificates").mkdir(exist_ok=True)
Path("certificates/auxiliary_cm_curve_rank0.json").write_text(
    json.dumps(out, indent=2, sort_keys=True) + "\n"
)
print(json.dumps(out, sort_keys=True))
