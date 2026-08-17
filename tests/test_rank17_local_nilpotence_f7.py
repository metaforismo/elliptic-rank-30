from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "rank17_local_nilpotence_f7.py"
CERTIFICATE = ROOT / "certificates" / "rank17_local_nilpotence_f7.json"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_local_nilpotence_f7", SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ImportError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Rank17LocalNilpotenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.fresh = cls.module.compute_certificate()
        cls.committed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    def test_committed_certificate_matches_recomputation(self):
        self.assertEqual(self.committed, self.fresh)

    def test_initial_ideal_is_maximal_primary(self):
        consequence = self.fresh["mathematical_consequence"]
        self.assertEqual(
            consequence["initial_ideal_contains"],
            ["b^3", "w^4", "t^12"],
        )
        self.assertEqual(consequence["local_krull_dimension"], 0)
        self.assertEqual(consequence["local_length_upper_bound"], 144)
        self.assertFalse(consequence["reduced_local_branch_through_seed"])

    def test_membership_targets_are_exact(self):
        targets = {
            record["target"]: tuple(record["target_monomial"])
            for record in self.fresh["membership_certificates"]
        }
        self.assertEqual(
            targets,
            {"b^3": (3, 0, 0), "w^4": (0, 4, 0), "t^12": (0, 0, 12)},
        )

    def test_truth_boundary_is_explicit(self):
        status = self.fresh["truth_status"]
        self.assertIn("over F_7", status)
        self.assertIn("no p-adic", status)
        self.assertTrue(
            any("vertical Z_7 lift" in item for item in self.fresh["limitations"])
        )
        self.assertTrue(
            any("rank-30" in item for item in self.fresh["limitations"])
        )


if __name__ == "__main__":
    unittest.main()
