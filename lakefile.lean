import Lake

open Lake DSL

package EllipticRank30 where
  leanOptions := #[
    ⟨`autoImplicit, false⟩,
    ⟨`warningAsError, true⟩
  ]
  moreServerOptions := #[⟨`linter.unusedVariables, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.32.2"

@[default_target]
lean_lib EllipticRank30 where
  roots := #[`EllipticRank30]
