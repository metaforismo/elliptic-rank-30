# The two degenerate degree-(4,7) remainders are impossible

**Status:** `proved` intermediate obstruction.

For

\[
Q(v)^2-v^2L(v)^3=v^3-Sv^2+3v+1, \tag{1}
\]

with \(\deg L=4\), \(\deg Q=7\), the cubic on the right is repeated-root
only for

\[
S=-3\quad\text{or}\quad S=\frac{15}{4}.
\]

Neither value admits a solution.

## Triple root: \(S=-3\)

Here

\[
R=(v+1)^3.
\]

Let

\[
q=\operatorname{ord}_{v=-1}Q,\qquad
\ell=\operatorname{ord}_{v=-1}L.
\]

If \(q=\ell=0\), the order three on the right is produced by cancellation.
The three Mason--Stothers radical supports then have degrees at most

\[
7,\quad5,\quad1.
\]

But the two large terms have degree \(14\), so Mason would require

\[
14\le7+5+1-1=12,
\]

a contradiction.

For a positive local order, the valuation of the difference can equal three
only when

\[
\ell=1,\qquad q\ge2.
\]

The full common divisor of the two large terms is then \((v+1)^3=R\).
After division, both terms have degree \(11\), while the radical has degree at
most

\[
(1+7-q)+4=12-q.
\]

Mason would give

\[
11\le(12-q)-1=11-q,
\]

impossible for \(q\ge2\).

## Double root: \(S=15/4\)

Here

\[
R=(v-2)^2\left(v+\frac14\right).
\]

If \(Q,L\) are units at \(v=2\), the radical degrees are at most

\[
7,\quad5,\quad2,
\]

and Mason would require

\[
14\le7+5+2-1=13,
\]

again impossible.

For a positive local order, an order-two remainder can occur only when

\[
\operatorname{ord}_{2}Q=1,\qquad
\operatorname{ord}_{2}L=\ell\ge1.
\]

Removing the full common divisor \((v-2)^2\) leaves two degree-12 terms.  Their
radical together with the remaining linear remainder has degree at most

\[
6+(2+4-\ell)+1=13-\ell.
\]

Mason would give

\[
12\le12-\ell,
\]

which is impossible for \(\ell\ge1\).

The simple root \(v=-1/4\) cannot contribute a common factor: a common zero of
\(Q\) and \(L\) would make the left side vanish to order at least two, whereas
\(R\) has order one there.

## Consequence

Every degree-(4,7) solution for the target coefficient pattern must therefore
belong to the squarefree nondegenerate Belyi passport

\[
(2^7),\quad(3^4\,2),\quad(11\,1^3).
\]

Once that finite passport is arithmetically resolved, Mason--Stothers and the
complete degree-(2,4) classification close **all polynomial identities of this
shape**.  Controlled-pole rational functions remain a different research
space.

## Reproduction

```bash
python3 research/degree47_degenerate_obstruction.py \
  --check certificates/degree47_degenerate_obstruction.json
python3 -m unittest tests.test_degree47_degenerate_obstruction -v
```
