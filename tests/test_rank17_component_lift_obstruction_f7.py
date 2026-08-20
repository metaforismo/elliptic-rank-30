from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "rank17_component_lift_obstruction_f7.py"
CERTIFICATE = ROOT / "certificates" / "rank17_component_lift_obstruction_f7.json"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_component_lift_obstruction_f7", SCRIPT
    )
    if spec is None or spec.loader is None:
        raise ImportError(SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Rank17ComponentLiftObstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()
        cls.fresh = cls.module.compute_certificate()
        cls.committed = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    def test_committed_certificate_matches_recomputation(self):
        self.assertEqual(self.committed, self.fresh)

    def test_full_radicalized_jacobian_has_column_rank_17(self):
        linear = self.fresh["linearized_lift_problem"]
        self.assertEqual(linear["jacobian_shape"], [25, 17])
        self.assertEqual(linear["jacobian_rank_mod_7"], 17)
        self.assertEqual(linear["pivot_columns"], list(range(17)))

    def test_left_kernel_witness_obstructs_mod_49(self):
        linear = self.fresh["linearized_lift_problem"]
        self.assertEqual(linear["witness_times_jacobian_mod_7"], [0] * 17)
        self.assertEqual(linear["witness_dot_rhs_mod_7"], 3)
        consequence = self.fresh["mathematical_consequence"]
        self.assertFalse(consequence["mod_49_lift_exists"])
        self.assertFalse(consequence["Z_7_lift_reducing_to_seed_exists"])

    def test_component_chart_is_explicit(self):
        seed = self.fresh["seed"]
        self.assertEqual(seed["component_tangent_ratios_mod_7"], [1, 1])
        self.assertEqual(len(self.fresh["radical_component_equations"]), 6)

    def test_truth_boundary_is_explicit(self):
        self.assertIn("no global rank-30 conclusion", self.fresh["truth_status"])
        limitations = self.fresh["limitations"]
        self.assertTrue(any("Other F_7" in item for item in limitations))
        self.assertTrue(any("IV" in item for item in limitations))
        self.assertTrue(any("rank at least 29" in item for item in limitations))


if __name__ == "__main__":
    unittest.main()
