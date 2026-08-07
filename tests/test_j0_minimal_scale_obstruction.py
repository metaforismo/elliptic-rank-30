import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "certificates" / "j0_minimal_scale_obstruction.json"
SCRIPT = ROOT / "research" / "j0_minimal_scale_obstruction_certificate.py"


class J0MinimalScaleObstructionTest(unittest.TestCase):
    def test_exact_certificate(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERT.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertEqual(data["scale_cube_classes_in_Qstar_mod_cubes"], ["1", "2"])
        self.assertEqual(data["K_rational_rho_values"], ["-1/2"])
        self.assertEqual(data["paired_base_rank_certificate"]["cube_class_1_exact_rank"], 0)
        self.assertEqual(data["paired_base_rank_certificate"]["cube_class_2_exact_rank"], 0)
        self.assertIn("restricted-family theorem", data["truth_status"])


if __name__ == "__main__":
    unittest.main()
