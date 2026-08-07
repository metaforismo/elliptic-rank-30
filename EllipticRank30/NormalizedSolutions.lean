import Mathlib

namespace EllipticRank30

def normalizedSystem (rho s g d : ℚ) : Prop :=
  -3 * g ^ 2 = rho ^ 3 * s ∧
  -2 * g * d = rho ^ 2 * s - 9 ∧
  -3 * d ^ 2 - 18 * g = 3 * rho * s + 405 / 4 ∧
  s = 81 - 18 * d

theorem normalizedSolutionOne :
    normalizedSystem (-1 / 2) 54 (-3 / 2) (3 / 2) := by
  norm_num [normalizedSystem]

theorem normalizedSolutionTwo :
    normalizedSystem (-1 / 2) 216 3 (-15 / 2) := by
  norm_num [normalizedSystem]

theorem bZeroCoefficientMismatch :
    (-315 / 4 : ℚ) ≠ 405 / 4 := by
  norm_num

end EllipticRank30
