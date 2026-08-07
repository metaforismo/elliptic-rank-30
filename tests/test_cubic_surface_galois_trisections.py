import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "cubic_surface_galois_trisection_certificate.py"
CERTIFICATE = ROOT / "certificates" / "cubic_surface_galois_trisections.json"


class CubicSurfaceGaloisTrisectionTest(unittest.TestCase):
    def test_regeneration_and_counts(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertEqual(data["e8_root_count"], 240)
        self.assertEqual(data["inner_product_one_neighbors_per_root"], 56)
        self.assertEqual(data["unordered_inner_product_one_root_pairs"], 6720)
        self.assertEqual(data["e8_norm_six_vector_count"], 6720)
        self.assertTrue(data["unique_unordered_root_pair_for_every_norm_six_vector"])
        self.assertEqual(data["intersection_checks"]["S_alpha_dot_S_beta"], 0)
        self.assertEqual(data["intersection_checks"]["K_X_square_after_three_contractions"], 3)


if __name__ == "__main__":
    unittest.main()
