import json
import subprocess
import sys
import unittest
from pathlib import Path


class TorsionPacketCapacityTest(unittest.TestCase):
    def test_certificate(self):
        subprocess.run([sys.executable,'research/torsion_packet_capacity_certificate.py'],check=True)
        data=json.loads(Path('certificates/torsion_packet_capacity.json').read_text())
        self.assertEqual(data['status'],'pass')
        self.assertEqual(data['two_torsion']['smooth_packet_rank_upper_bound'],22)
        self.assertEqual(data['three_torsion']['smooth_packet_rank_upper_bound'],14)


if __name__=='__main__': unittest.main()
