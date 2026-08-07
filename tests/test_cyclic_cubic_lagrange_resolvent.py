import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "cyclic_cubic_lagrange_resolvent_certificate.py"
CERTIFICATE = ROOT / "certificates" / "cyclic_cubic_lagrange_resolvent.json"


class CyclicCubicLagrangeResolventTest(unittest.TestCase):
    def test_symbolic_certificate(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertTrue(all(data["symbolic_checks"].values()))
        self.assertEqual(data["lagrange_cubes"]["sum"], "-27q")
        self.assertEqual(data["lagrange_cubes"]["product"], "-27p^3")


if __name__ == "__main__":
    unittest.main()
