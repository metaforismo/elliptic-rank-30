import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "certificates" / "cyclic_trisection_trace_code.json"
SCRIPT = ROOT / "research" / "cyclic_trisection_trace_code_certificate.py"


class CyclicTrisectionTraceCodeTest(unittest.TestCase):
    def test_compact_certificate(self):
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["status"], "pass")
        self.assertEqual(data["e8_norm6_vector_count"], 6720)
        self.assertEqual(data["nonzero_isotropic_vectors_mod3"], 2240)
        self.assertEqual(data["norm6_lifts_per_nonzero_isotropic_vector"], 3)
        self.assertEqual(data["maximal_totally_isotropic_4_spaces"], 2240)
        self.assertEqual(data["norm6_trisection_classes_per_maximal_space"], 240)
        self.assertEqual(data["orbit_height_formula"]["projected_gram"], [[4, -2], [-2, 4]])
        self.assertEqual(data["orbit_height_formula"]["difference_gram"], [[12, -6], [-6, 12]])

    def test_dependency_free_regeneration(self):
        subprocess.run([sys.executable, str(SCRIPT)], cwd=ROOT, check=True)
        data = json.loads(CERTIFICATE.read_text())
        self.assertEqual(data["shell_counts_for_intersection_m_0_through_7"]["0"], 6720)
        self.assertEqual(data["shell_counts_for_intersection_m_0_through_7"]["1"], 60480)


if __name__ == "__main__":
    unittest.main()
