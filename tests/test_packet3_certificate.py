from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "packet3_certificate.py"
CERTIFICATE_PATH = ROOT / "certificates" / "packet3_certificate.json"
SPEC = importlib.util.spec_from_file_location("packet3_certificate", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class Packet3CertificateTests(unittest.TestCase):
    def test_exact_certificate_matches_committed_json(self) -> None:
        computed = MODULE.compute_certificate()
        committed = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
        self.assertEqual(computed, committed)

    def test_cover_and_rank_conclusion(self) -> None:
        result = MODULE.compute_certificate()
        self.assertEqual(result["cover"]["degree"], 4)
        self.assertEqual(result["cover"]["genus"], 0)
        self.assertEqual(result["squareclasses"]["rank_over_F2"], 2)
        self.assertEqual(len(result["sections"]), 3)
        self.assertEqual(len(result["non_torsion_specialization"]["twists"]), 3)
        self.assertEqual(result["conditional_assumptions"], [])
        self.assertEqual(
            result["claim"],
            "For the displayed E/Q(t) and biquadratic genus-zero cover L/Q(t), "
            "rank E(L) >= rank E(Q(t)) + 3.",
        )


if __name__ == "__main__":
    unittest.main()
