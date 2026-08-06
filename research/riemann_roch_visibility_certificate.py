#!/usr/bin/env python3
"""Machine-readable arithmetic part of the one-function visibility barrier."""
import json
from pathlib import Path


def pole_order(deg_a: int, deg_b: int) -> int:
    return max(2*deg_a, 2*deg_b+3)


def main():
    near_square_order=pole_order(15,-10**9)  # B=0, A degree 15
    norm31_order=pole_order(15,14)
    assert near_square_order==30
    assert norm31_order==31
    out={
      'status':'pass',
      'principal_divisor_relation':'sum of N rational zeros equals N*O, hence sum(P_i)=O in E(k)',
      'near_square_degree_15':{
        'pole_order':near_square_order,
        'visible_point_count':30,
        'forced_relation_count_at_least':1,
        'rank_upper_bound_from_these_points':29,
      },
      'minimal_single_function_rank30_target':{
        'deg_A_max':15,
        'deg_B':14,
        'pole_order':norm31_order,
        'norm_degree':31,
        'visible_point_count':31,
        'forced_relation_count_at_least':1,
        'maximum_possible_span_after_forced_relation':30,
        'norm_form':'A(x)^2-R(x)B(x)^2',
      },
      'truth_note':'the 31-zero norm identity is a construction target, not an independence proof',
    }
    Path('certificates').mkdir(exist_ok=True)
    Path('certificates/riemann_roch_visibility_barrier.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__': main()
