#!/usr/bin/env python3
"""Arithmetic certificate for the torsion packet-capacity inequalities."""
import json
from pathlib import Path


def minimum_root_cost(order: int, divisor_degree: int) -> dict:
    # For multiplicities m_i>0 with sum d, sum(order*m_i-1)
    # is minimized when the support has maximal size d.
    support_max=divisor_degree
    cost=order*divisor_degree-support_max
    return {
        'torsion_order':order,
        'divisor_degree':divisor_degree,
        'maximal_support_size':support_max,
        'minimum_root_cost':cost,
        'base_rank_upper_bound':8-cost,
        'smooth_packet_rank_upper_bound':4*(8-cost)+6,
    }


def main():
    two=minimum_root_cost(2,4)
    three=minimum_root_cost(3,3)
    assert two['minimum_root_cost']==4
    assert two['base_rank_upper_bound']==4
    assert two['smooth_packet_rank_upper_bound']==22
    assert three['minimum_root_cost']==6
    assert three['base_rank_upper_bound']==2
    assert three['smooth_packet_rank_upper_bound']==14
    out={
      'status':'pass',
      'theorem':'small rational torsion prevents rank 30 in a smooth three-branch V4 packet on a rational elliptic surface',
      'two_torsion':two,
      'three_torsion':three,
      'scope':'restricted packet obstruction, not a universal rank bound',
    }
    Path('certificates').mkdir(exist_ok=True)
    Path('certificates/torsion_packet_capacity.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__': main()
