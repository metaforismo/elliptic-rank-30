import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "eisenstein_packet_determinant_certificate.py"
CERTIFICATE = ROOT / "certificates" / "eisenstein_packet_determinant.json"


class EisensteinPacketDeterminantTest(unittest.TestCase):
    def test_exact_certificate(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertTrue(all(data["exact_arithmetic_checks"].values()))
        self.assertEqual(data["one_orbit"]["difference_determinant"], 108)
        self.assertEqual(data["one_orbit"]["projected_determinant"], 12)
        self.assertEqual(data["general_formulas"]["rank30_n11"], "Norm(det H)/27")


if __name__ == "__main__":
    unittest.main()
