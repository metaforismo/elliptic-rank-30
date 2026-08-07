import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "search" / "e8_polynomial_twist_denominator_sieve.py"


class E8PolynomialTwistDenominatorSieveTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("e8_denominator_sieve", SCRIPT)
        cls.module = importlib.util.module_from_spec(spec)
        assert spec.loader is not None
        spec.loader.exec_module(cls.module)

    def test_prime_41_obstructs_all_45_branch_pairs(self):
        result = self.module.scan_prime(41, 4)
        self.assertEqual(result["pair_count"], 45)
        self.assertTrue(result["all_pairs_empty"])
        self.assertEqual(result["nonempty_pair_count"], 0)
        self.assertEqual(result["total_x_classes"], 0)


if __name__ == "__main__":
    unittest.main()
