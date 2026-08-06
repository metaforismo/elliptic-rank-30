import json
import subprocess
import sys
import unittest
from pathlib import Path


class RiemannRochVisibilityTest(unittest.TestCase):
    def test_certificate(self):
        subprocess.run([sys.executable,'research/riemann_roch_visibility_certificate.py'],check=True)
        d=json.loads(Path('certificates/riemann_roch_visibility_barrier.json').read_text())
        self.assertEqual(d['status'],'pass')
        self.assertEqual(d['near_square_degree_15']['rank_upper_bound_from_these_points'],29)
        self.assertEqual(d['minimal_single_function_rank30_target']['pole_order'],31)
        self.assertEqual(d['minimal_single_function_rank30_target']['maximum_possible_span_after_forced_relation'],30)


if __name__=='__main__': unittest.main()
