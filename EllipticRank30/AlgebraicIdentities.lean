import Mathlib

namespace EllipticRank30

theorem cubicBaseBirationalIdentity (c t u : ℚ)
    (h : u ^ 3 = c * t * (t - 1)) :
    (4 * c ^ 2 * (2 * t - 1)) ^ 2 =
      (4 * c * u) ^ 3 + 16 * c ^ 4 := by
  calc
    (4 * c ^ 2 * (2 * t - 1)) ^ 2 =
        16 * c ^ 4 * (4 * t * (t - 1) + 1) := by ring
    _ = 64 * c ^ 3 * u ^ 3 + 16 * c ^ 4 := by
      rw [h]
      ring
    _ = (4 * c * u) ^ 3 + 16 * c ^ 4 := by ring

theorem markedQuadraticPositive (mu : ℝ) :
    0 < mu ^ 2 - mu + 1 := by
  nlinarith [sq_nonneg (mu - (1 / 2 : ℝ))]

theorem minimalCharacterOneRealContradiction (mu e : ℝ)
    (h : e ^ 2 = -(mu ^ 2 - mu + 1) ^ 3) : False := by
  have hq : 0 < mu ^ 2 - mu + 1 := markedQuadraticPositive mu
  have hc : 0 < (mu ^ 2 - mu + 1) ^ 3 := pow_pos hq 3
  have hs : 0 <= e ^ 2 := sq_nonneg e
  nlinarith

theorem muTwoCoefficientFactorization (v : ℚ) :
    (729 / 4 : ℚ) * v ^ 2 - 27 * (v + 1) ^ 3 =
      -(27 / 4 : ℚ) * (v - 2) ^ 2 * (4 * v + 1) := by
  ring

theorem muTwoTwistFactorization (v : ℚ) :
    (1 + 4 * v) ^ 3 *
        ((729 / 4 : ℚ) * v ^ 2 - 27 * (v + 1) ^ 3) =
      -(27 / 4 : ℚ) * (v - 2) ^ 2 * (4 * v + 1) ^ 4 := by
  ring

end EllipticRank30
