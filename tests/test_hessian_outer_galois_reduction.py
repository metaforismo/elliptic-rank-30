import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "hessian_outer_galois_certificate.py"
CERTIFICATE = ROOT / "certificates" / "hessian_outer_galois_reduction.json"


class HessianOuterGaloisReductionTest(unittest.TestCase):
    def test_symbolic_certificate(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertTrue(all(data["symbolic_checks"].values()))
        self.assertEqual(
            data["generic_cyclic_hyperplanes_per_trace_class_upper_bound"], 8
        )
        self.assertEqual(
            data["generic_candidates_per_240_class_trace_code_upper_bound"], 1920
        )
        self.assertEqual(data["cyclicity_condition"], "3*a*Q - L^2 = 0")


if __name__ == "__main__":
    unittest.main()
