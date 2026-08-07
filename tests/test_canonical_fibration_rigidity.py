import json
import subprocess
import sys
import unittest
from pathlib import Path


class CanonicalFibrationRigidityTest(unittest.TestCase):
    def test_certificate(self):
        subprocess.run([sys.executable,'research/canonical_fibration_rigidity_certificate.py'],check=True)
        d=json.loads(Path('certificates/canonical_fibration_rigidity.json').read_text())
        self.assertEqual(d['status'],'pass')
        self.assertFalse(d['four_II_star_example']['alternative_fibration_allowed'])
        self.assertEqual(d['four_II_star_example']['root_rank'],32)
        self.assertEqual(d['exception'],'chi=2 (K3), where K_S=0 and neighbour fibrations may exist')


if __name__=='__main__': unittest.main()
