from __future__ import annotations

import importlib.util
import unittest

SYMPY_AVAILABLE = importlib.util.find_spec("sympy") is not None

if SYMPY_AVAILABLE:
    import sympy as sp

    from research.degree47_target_elimination import (
        c,
        d,
        b,
        q_coefficients,
        rational_cube_root,
        reduced_equations,
    )


@unittest.skipUnless(SYMPY_AVAILABLE, "SymPy is an optional exact research dependency")
class Degree47TargetEliminationTests(unittest.TestCase):
    def test_wronskian_recursion_and_reduced_equations(self) -> None:
        v = sp.symbols("v")
        qs = q_coefficients(b, d)
        L = 1 + b * v**2 + c * v**3 + d * v**4
        Q = sum(qs[index] * v**index for index in range(8))
        wronskian = sp.Poly(
            sp.expand(2 * L * Q + 3 * v * sp.diff(L, v) * Q - 2 * v * L * sp.diff(Q, v)),
            v,
            domain=sp.QQ.frac_field(b, d),
        )
        self.assertEqual(wronskian.nth(0), 2)
        for degree in range(1, 8):
            self.assertEqual(sp.factor(wronskian.nth(degree)), 0)
        self.assertEqual(sp.factor(wronskian.nth(11)), 0)

        equations = reduced_equations()
        for degree, equation in zip((8, 9, 10), equations):
            coefficient = sp.factor(wronskian.nth(degree))
            self.assertEqual(sp.factor(coefficient / equation).free_symbols, set())

    def test_low_target_coefficient(self) -> None:
        qs = q_coefficients(b, d)
        self.assertEqual(sp.factor(2 * qs[3] + 3 * qs[2]), 1)

    def test_exact_rational_cube_root(self) -> None:
        self.assertEqual(rational_cube_root(sp.Rational(8, 27)), sp.Rational(2, 3))
        self.assertEqual(rational_cube_root(sp.Rational(-125, 64)), sp.Rational(-5, 4))
        self.assertIsNone(rational_cube_root(sp.Rational(2, 3)))


if __name__ == "__main__":
    unittest.main()
