# Canonical rigidity outside the K3 case

## Theorem

Let

\[
\pi:S\longrightarrow \mathbf P^1
\]

be a relatively minimal elliptic surface over a characteristic-zero field,
with a section and no multiple fibres.  Put

\[
d=\chi(\mathcal O_S).
\]

If `d != 2`, then every genus-one fibration on `S` has the same primitive
fibre class as `pi`.  In particular, up to an automorphism of the base, `pi`
is the unique genus-one fibration on `S`.

### Proof

The canonical bundle formula gives

\[
K_S\sim (d-2)F,
\]

where `F` is the fibre class of `pi`.  Let `F'` be the primitive nef class of
the general fibre of another genus-one fibration.  Adjunction gives

\[
(F')^2=0,
\qquad K_S\cdot F'=0.
\]

When `d != 2`, the canonical bundle formula implies

\[
F\cdot F'=0.
\]

The two classes `F` and `F'` are nonzero nef isotropic classes.  By the Hodge
index theorem, two such classes with zero intersection are proportional in
`NS(S)_R`.  Since both are primitive integral fibre classes, they are equal.
The complete linear systems therefore define the same fibration, up to a
change of coordinate on the base.

The exceptional case is `d=2`, when `K_S=0`: elliptic K3 surfaces may admit
many inequivalent elliptic fibrations and neighbour moves are genuinely
available.

## A long and useful failed route

For `d=4`, consider the explicit family

\[
y^2=x^3+c\,h(t)^4x+h(t)^5q(t),
\qquad \deg h=4,\quad \deg q\le4.
\]

When the four roots of `h` are distinct and `q` is nonzero there, the surface
has four fibres of type `II*`; generically the remaining discriminant factor
has degree eight and is square-free.  The trivial lattice contains

\[
U\oplus E_8(-1)^{\oplus4}.
\]

Because

\[
U\oplus E_8(-1)^{\oplus4}
\]

is the unique even unimodular lattice of signature `(1,33)`, it is abstractly
isometric to `U+L(-1)` for any even unimodular positive-definite rank-32
lattice `L`, including rootless examples.  This makes it look as though a new
elliptic fibration could convert the four vertical `E8` lattices into a
rank-32 Mordell--Weil lattice.

The theorem shows exactly why that argument fails.  Since `d=4`,

\[
K_S=2F,
\]

so the canonical class fixes the fibre ray.  There is no alternative elliptic
fibration in which the four `E8` root lattices become sections.

## Additional lattice check

On an elliptic surface with `d=4` and no reducible fibres, every nonzero
section satisfies

\[
\langle P,P\rangle=2d+2(P\cdot O)\ge8.
\]

Hence a hypothetical unimodular rank-32 Mordell--Weil lattice on such a
surface would be an even unimodular rank-32 lattice of minimum at least eight.
The modular-form bound for even unimodular lattices in dimension 32 gives
minimum at most four.  Thus the particular unimodular lattice conversion is
also impossible from the height side.

## Search consequence

Root lattices on a high-`d` elliptic surface cannot be recycled into
Mordell--Weil directions by changing fibrations.  Inverse lattice engineering
must construct the desired section lattice in the canonical elliptic
fibration itself.  The freedom to change elliptic fibrations is a special
advantage of K3 surfaces, not a general feature of higher-arithmetic-genus
elliptic surfaces.

This is a restricted obstruction to a proposed rank-32 construction.  It is
not an upper bound on ranks over `Q`.
