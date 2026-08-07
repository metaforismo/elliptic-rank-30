import Mathlib

namespace EllipticRank30

def rankCapacity (chi rootRank : Nat) : Nat :=
  10 * chi - 2 - rootRank

def sixChannelCapacityVector : List Nat :=
  [rankCapacity 1 4, rankCapacity 1 2, rankCapacity 2 14,
   rankCapacity 1 4, rankCapacity 2 12, rankCapacity 2 12]

theorem sixChannelCapacityVector_eq :
    sixChannelCapacityVector = [4, 6, 4, 4, 6, 6] := by
  native_decide

theorem sixChannelTotalCapacity :
    sixChannelCapacityVector.sum = 30 := by
  native_decide

theorem rationalSurfaceRankSum : 4 + 6 + 4 = 14 := by
  norm_num

theorem k3RankCapacitySum : 4 + 6 + 6 = 16 := by
  norm_num

theorem rankThirtySplit : 14 + 16 = 30 := by
  norm_num

theorem muTwoInvariantEuler : 3 * 4 + 3 * 2 + 6 = 24 := by
  norm_num

theorem muTwoInvariantRootRank : 3 * 2 + 4 = 10 := by
  norm_num

theorem muTwoInvariantCapacity : rankCapacity 2 10 = 8 := by
  norm_num [rankCapacity]

theorem muTwoTwistEuler : 3 * 4 + 3 * 8 = 36 := by
  norm_num

theorem muTwoTwistRootRank : 3 * 2 + 3 * 6 = 24 := by
  norm_num

theorem muTwoTwistCapacity : rankCapacity 3 24 = 4 := by
  norm_num [rankCapacity]

theorem muTwoCombinedCeiling : 8 + 4 = 12 := by
  norm_num

theorem muTwoRequiredJump : 30 - 12 = 18 := by
  norm_num

theorem k20Euler : 3 * 2 + 8 + 10 = 24 := by
  norm_num

theorem k20RootRank : 6 + 8 = 14 := by
  norm_num

theorem k20Capacity : rankCapacity 2 14 = 4 := by
  norm_num [rankCapacity]

theorem k11Euler : 3 * 2 + 4 + 6 + 8 = 24 := by
  norm_num

theorem k11RootRank : 2 + 4 + 6 = 12 := by
  norm_num

theorem k11Capacity : rankCapacity 2 12 = 6 := by
  norm_num [rankCapacity]

theorem k21Euler : 3 * 2 + 8 + 6 + 4 = 24 := by
  norm_num

theorem k21RootRank : 6 + 4 + 2 = 12 := by
  norm_num

theorem k21Capacity : rankCapacity 2 12 = 6 := by
  norm_num [rankCapacity]

theorem nonzeroIsotropicVectorCount :
    (3 ^ 3 + 1) * (3 ^ 4 - 1) = 2240 := by
  norm_num

theorem maximalIsotropicFourSpaceCount :
    2 * (3 + 1) * (3 ^ 2 + 1) * (3 ^ 3 + 1) = 2240 := by
  norm_num

theorem normSixShellCount :
    240 * (1 + 3 ^ 3) = 6720 := by
  norm_num

theorem traceCodeClassCount : 80 * 3 = 240 := by
  norm_num

end EllipticRank30
