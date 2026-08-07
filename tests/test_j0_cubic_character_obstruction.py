import json
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "research" / "j0_cubic_character_obstruction_certificate.py"
CERTIFICATE = ROOT / "certificates" / "j0_cubic_character_obstruction.json"


class J0CubicCharacterObstructionTest(unittest.TestCase):
    def test_symbolic_certificate(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertTrue(all(data["symbolic_checks"].values()))
        self.assertTrue(data["character_one"]["no_Q_rational_eigen_section"])
        self.assertIn("empty", data["character_two"]["conclusion"])
        self.assertEqual(data["cm_compatible_locus"]["rational_point"], ["2", "1"])


if __name__ == "__main__":
    unittest.main()
