#!/usr/bin/env python3
"""Machine-readable consequences of canonical fibre-class rigidity.

The geometric theorem is proved in canonical_fibration_rigidity.md.  This
script checks the numerical specialisations used by the failed rank-32 lattice
route and records their exact scope.
"""
import json
from pathlib import Path


def main():
    chi=4
    canonical_multiple=chi-2
    ii_star_root_rank=8
    four_fibre_root_rank=4*ii_star_root_rank
    trivial_rank=2+four_fibre_root_rank
    h11=10*chi
    apparent_mw_room=h11-trivial_rank
    assert canonical_multiple==2
    assert four_fibre_root_rank==32
    assert trivial_rank==34
    assert h11==40
    assert apparent_mw_room==6
    out={
      'status':'pass',
      'theorem':'for chi(O_S) != 2 the canonical class fixes the primitive genus-one fibre ray',
      'chi':chi,
      'canonical_class':'K_S = 2 F',
      'four_II_star_example':{
        'root_lattice':'E8^4',
        'root_rank':four_fibre_root_rank,
        'trivial_lattice_rank':trivial_rank,
        'h11':h11,
        'same_fibration_MW_upper_bound':apparent_mw_room,
        'alternative_fibration_allowed':False,
      },
      'exception':'chi=2 (K3), where K_S=0 and neighbour fibrations may exist',
      'scope':'restricted obstruction to recycling vertical root lattices; not a universal rank bound',
    }
    Path('certificates').mkdir(exist_ok=True)
    Path('certificates/canonical_fibration_rigidity.json').write_text(
        json.dumps(out,indent=2,sort_keys=True)+'\n'
    )
    print(json.dumps(out,sort_keys=True))

if __name__=='__main__': main()
