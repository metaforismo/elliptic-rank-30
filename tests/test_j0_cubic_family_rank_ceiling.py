import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "certificates" / "j0_cubic_family_rank_ceiling.json"
SCRIPT = ROOT / "research" / "j0_cubic_family_rank_ceiling_certificate.py"


class J0CubicFamilyRankCeilingTest(unittest.TestCase):
    def test_rank_ceiling(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERT.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertEqual(data["invariant_surface"]["geometric_mordell_weil_rank_upper_bound"], 8)
        self.assertEqual(data["quadratic_twist_surface"]["geometric_mordell_weil_rank_upper_bound"], 4)
        self.assertEqual(data["combined_geometric_rank_upper_bound"], 12)
        self.assertEqual(data["specialization_jump_needed_for_rank_30"], 18)


if __name__ == "__main__":
    unittest.main()
