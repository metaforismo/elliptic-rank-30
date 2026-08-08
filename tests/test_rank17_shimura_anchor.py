from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "rank17_shimura_anchor.py"
CERTIFICATE_PATH = ROOT / "certificates" / "rank17_shimura_anchor.json"

SPEC = importlib.util.spec_from_file_location(
    "rank17_shimura_anchor", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class Rank17ShimuraAnchorTests(unittest.TestCase):
    def test_certificate_matches_committed_json(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)

    def test_exact_anchor_claims(self) -> None:
        result = MODULE.compute_certificate()
        self.assertEqual(result["shimura_metadata"]["N"], 474)
        self.assertEqual(result["shimura_metadata"]["factorization"], [6, 79])
        self.assertEqual(result["curve"]["genus"], 2)
        self.assertTrue(result["curve"]["squarefree"])
        self.assertEqual(len(result["rational_affine_points"]), 8)
        self.assertEqual(
            len(
                result["curve"][
                    "rational_points_at_infinity_weighted_projective_T_Z_U"
                ]
            ),
            2,
        )
        self.assertEqual(
            result["bielliptic_involution"]["quotient_genus"], 1
        )
        self.assertTrue(
            result["bielliptic_involution"]["quotient_nonsingular"]
        )
        self.assertEqual(result["conditional_assumptions"], [])

    def test_unproved_moduli_claims_remain_explicit(self) -> None:
        result = MODULE.compute_certificate()
        claims = result["source_claims_not_proved_here"]
        self.assertEqual(len(claims), 3)
        self.assertTrue(any("non-CM" in claim for claim in claims))
        self.assertTrue(any("rank 17" in claim for claim in claims))
        self.assertIn("missing_bridge", result)


if __name__ == "__main__":
    unittest.main()
