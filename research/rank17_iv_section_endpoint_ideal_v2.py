#!/usr/bin/env python3
"""Corrected launcher for rank17_iv_section_endpoint_ideal.py.

The underlying implementation is reused verbatim, but the one-variable time
ring is forced to be a genuine univariate polynomial ring.  This guarantees
that H[k], quo_rem, degree, and resultant have their intended univariate
meaning.  Multivariate coefficient rings continue to use the original Sage
PolynomialRing constructor.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

from sage.all import PolynomialRing as SagePolynomialRing

HERE = Path(__file__).resolve().parent
IMPLEMENTATION = HERE / "rank17_iv_section_endpoint_ideal.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "rank17_iv_section_endpoint_ideal", IMPLEMENTATION
    )
    if spec is None or spec.loader is None:
        raise ImportError(IMPLEMENTATION)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def corrected_polynomial_ring(base_ring, *args, **kwargs):
    names = kwargs.get("names")
    if names == ("t",) and "order" not in kwargs and not args:
        return SagePolynomialRing(base_ring, "t")
    return SagePolynomialRing(base_ring, *args, **kwargs)


def main() -> None:
    module = load_module()
    module.PolynomialRing = corrected_polynomial_ring
    module.main()


if __name__ == "__main__":
    main()
