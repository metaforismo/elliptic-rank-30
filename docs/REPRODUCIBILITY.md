# Reproducibility

The final verifier must depend only on the curve equation and point list.
Search databases and cached rankings are never proof inputs.

Core commands:

```bash
python3 verify_exact.py
python3 rank_packet_obstruction_tests.py
python3 -m unittest discover -s tests -v
sage -python verify_sage.py
magma verify_magma.m
```

Every retained result records exact inputs, commands, software versions, proof
flags, SHA-256 hashes, and provenance.
