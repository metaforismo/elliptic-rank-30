import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / "certificates" / "j0_six_channel_capacity.json"
SCRIPT = ROOT / "research" / "j0_six_channel_capacity_certificate.py"


class J0SixChannelCapacityTest(unittest.TestCase):
    def test_capacity_vector(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERT.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertEqual(data["capacity_vector"], [4, 6, 4, 4, 6, 6])
        self.assertEqual(data["total_geometric_rank_capacity"], 30)
        self.assertEqual(data["automatic_rational_surface_rank_sum"], 14)
        self.assertEqual(data["k3_capacity_vector"], [4, 6, 6])


if __name__ == "__main__":
    unittest.main()
