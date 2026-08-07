import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "cyclic_cubic_hodge_certificate.py"
CERTIFICATE = ROOT / "certificates" / "cyclic_cubic_hodge_representation.json"


class CyclicCubicHodgeRepresentationTest(unittest.TestCase):
    def test_exact_invariants(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertEqual(data["surface_invariants"]["b2"], 38)
        self.assertEqual(data["surface_invariants"]["h11"], 32)
        self.assertEqual(
            data["H2_character_dimensions"],
            {"trivial": 10, "zeta3": 14, "zeta3_squared": 14},
        )
        self.assertEqual(
            data["H11_character_dimensions"],
            {"trivial": 10, "zeta3": 11, "zeta3_squared": 11},
        )
        self.assertEqual(data["mordell_weil"]["rank_ceiling"], 30)
        self.assertEqual(
            data["mordell_weil"]["rank30_equivalent_picard_number"], 32
        )


if __name__ == "__main__":
    unittest.main()
