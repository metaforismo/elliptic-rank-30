# Minimal C3-character sections on the marked j=0 cubic base

**Author:** Francesco Giannicola  
**Truth status:** new restricted obstruction and field-of-definition theorem;
no rank-30 curve is claimed.

## The cubic-base surface

Put

\[
q(t)=t^2-t+1,
\qquad
r(t)=t(t-1)=q(t)-1.
\]

For a rational marking parameter `mu`, define

\[
a=q(\mu),
\qquad
b=\frac{q(\mu)^3}{r(\mu)}.
\]

The marked rational `j=0` E8 surface is

\[
y^2=x^3+b^2r(t)^2-a^3q(t)^3.
\]

After the cyclic cubic base change

\[
u^3=r(t),
\]

the curve descends to the `u`-line:

\[
\boxed{
E_u:\ y^2=x^3+b^2u^6-a^3(u^3+1)^3.
}
\]

The automorphism `u -> zeta3 u` decomposes polynomial sections into three
characters.

## Character two is geometrically impossible in the minimal polynomial box

A section disjoint from the zero section in character two must have

\[
x=A u^2,
\qquad
y=C u^6+D u^3+E.
\]

The coefficient of `u^12` in the curve equation gives

\[
C^2=0,
\]

so `C=0`. The coefficient of `u^9` then gives

\[
0=-a^3.
\]

Since `a=q(mu)` is nonzero for every rational `mu`, this is impossible even
over an algebraic closure.

Therefore

\[
\boxed{
\text{the character-two minimal polynomial section scheme is empty.}
}
\]

Any character-two direction must have poles, meet the zero section, or arise
only after changing the geometric representation.

## Character one has a sharp field-of-definition obstruction

A character-one minimal polynomial section has

\[
x=u(Au^3+B),
\qquad
y=Cu^6+Du^3+E.
\]

Coefficient comparison gives

\[
\begin{aligned}
C^2&=A^3,\\
2CD&=3A^2B-a^3,\\
D^2+2CE&=3AB^2+b^2-3a^3,\\
2DE&=B^3-3a^3,\\
E^2&=-a^3.
\end{aligned}
\]

The final equation is unavoidable. Hence every such eigen-section is defined
over a field containing

\[
\boxed{\sqrt{-a}=\sqrt{-q(\mu)}.}
\]

For rational `mu`,

\[
q(\mu)=\left(\mu-\frac12\right)^2+\frac34>0,
\]

so no character-one eigen-section in the minimal polynomial box is defined
over `Q`.

This is not an obstruction to rational rank on the original cubic cover. The
rational trace-zero representation is obtained by descending a conjugate pair
of nontrivial character eigenspaces. The theorem identifies the exact constant
field in which the eigenvectors must first be constructed.

## CM-compatible parameters

The most favorable parameters satisfy

\[
q(\mu)=3s^2,
\qquad s\in\mathbf Q,
\]

because then

\[
\sqrt{-q(\mu)}=s\sqrt{-3}
\]

lies in the natural CM field `Q(sqrt(-3))`. The conic

\[
\mu^2-\mu+1=3s^2
\]

has the rational point `(mu,s)=(2,1)` and is rationally parametrizable.
Thus there is an infinite rational parameter family in which every minimal
character-one eigensection, if it exists, is defined over the same Eisenstein
field used by the cyclic packet.

## Decisive change of setting

The failed rational-coefficient search should not be enlarged blindly.
Instead:

1. restrict to the CM-compatible conic `q(mu)=3s^2`;
2. solve the five exact character-one equations over `Q(sqrt(-3))`;
3. compute the conjugate eigensection and its Galois descent to the rational
   trace-zero Mordell--Weil representation;
4. test whether the descended orbit contributes a new Hermitian norm-six
   direction;
5. move to controlled-pole sections only after the complete minimal scheme is
   exhausted.

The coefficient obstruction is verified by
`research/j0_cubic_character_obstruction_certificate.py`.
