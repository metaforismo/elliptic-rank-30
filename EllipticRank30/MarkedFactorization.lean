import Mathlib

namespace EllipticRank30

theorem markedCoefficientFactorization (mu v : ℚ)
    (hmu0 : mu ≠ 0) (hmu1 : mu ≠ 1) :
    let a := mu ^ 2 - mu + 1
    let b := a ^ 3 / (mu * (mu - 1))
    b ^ 2 * v ^ 2 - a ^ 3 * (v + 1) ^ 3 =
      -(a ^ 3) / (mu ^ 2 * (mu - 1) ^ 2) *
        (v - mu * (mu - 1)) *
        (mu ^ 2 * v + mu - 1) *
        ((mu - 1) ^ 2 * v - mu) := by
  dsimp
  field_simp [hmu0, sub_ne_zero.mpr hmu1]
  ring

theorem markedRootsProduct (mu : ℚ)
    (hmu0 : mu ≠ 0) (hmu1 : mu ≠ 1) :
    (mu * (mu - 1)) * ((1 - mu) / mu ^ 2) *
        (mu / (mu - 1) ^ 2) = -1 := by
  field_simp [hmu0, sub_ne_zero.mpr hmu1]
  ring

theorem markedRootsPairwiseSum (mu : ℚ)
    (hmu0 : mu ≠ 0) (hmu1 : mu ≠ 1) :
    let r1 := mu * (mu - 1)
    let r2 := (1 - mu) / mu ^ 2
    let r3 := mu / (mu - 1) ^ 2
    r1 * r2 + r1 * r3 + r2 * r3 = 3 := by
  dsimp
  field_simp [hmu0, sub_ne_zero.mpr hmu1]
  ring

end EllipticRank30
