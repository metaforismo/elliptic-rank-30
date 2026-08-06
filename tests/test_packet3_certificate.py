from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "research" / "packet3_certificate.py"
SPEC = importlib.util.spec_from_file_location("packet3_certificate", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Packet3CertificateTests(unittest.TestCase):
    def test_exact_certificate(self) -> None:
        result = MODULE.verify_packet()
        self.assertEqual(result["status"], "pass")
        self.assertEqual(result["characters"]["cover_genus"], 0)
        self.assertEqual(result["characters"]["geometric_dimension"], 2)
        self.assertEqual(
            result["rank_conclusion"],
            "rank E(L) >= rank E(Q(t)) + 3",
        )
        self.assertEqual(
            len(result["nontorsion_specialization"]["witnesses"]),
            3,
        )


if __name__ == "__main__":
    unittest.main()
