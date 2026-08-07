import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "search" / "enumerate_e8_mod3_trace_codes.py"


class E8Mod3TraceCodesTest(unittest.TestCase):
    def test_complete_enumeration(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "codes.json"
            subprocess.run(
                [sys.executable, str(SCRIPT), "--output", str(output)],
                cwd=ROOT,
                check=True,
            )
            data = json.loads(output.read_text())
            self.assertEqual(data["status"], "pass")
            self.assertEqual(data["projective_isotropic_points"], 1120)
            self.assertEqual(
                data["subspace_counts_by_vector_dimension"],
                {"1": 1120, "2": 36400, "3": 44800, "4": 2240},
            )
            self.assertEqual(data["maximal_space_count"], 2240)
            self.assertEqual(len(data["canonical_rref_bases_base3_encoded"]), 2240)
            self.assertEqual(len(data["basis_payload_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
