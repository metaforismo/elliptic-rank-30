# Local nilpotence at the rank-17 section-incidence point over F7

## Claim status

**Proved in characteristic 7.** This note does not assert a p-adic lift, a
characteristic-zero section, or a rank-30 curve.

## Setup

The normalized split semistable surface-section incidence system has 19 exact
polynomial equations in 17 variables. At the certified F7 point, its Jacobian
has rank 14. The rows and columns recorded in
`certificates/rank17_section_local_branch_f7.json` contain an invertible
14-by-14 minor.

The formal implicit-function theorem therefore eliminates those fourteen
regular coordinates uniquely in terms of

```text
a = y1 - 1
u = y2 - 4
v = y3 - 2.
```

The five remaining equations are computed in the completed local ring
`F7[[a,u,v]]` by
`research/rank17_local_three_variable_elimination.py`.

## Adapted coordinates

Set

```text
b = a - 2(u + v)
w = u + v
t = u.
```

Equivalently,

```text
a = b + 2w
u = t
v = w - t.
```

The tangent-cone calculation had already exposed the initial terms `b^3` and
`4 w^4`. The exact finite membership calculation in
`research/rank17_local_nilpotence_f7.py` additionally constructs polynomial
multipliers for the five residual series whose combinations have initial
terms

```text
b^3,
w^4,
t^12.
```

The committed JSON certificate stores every multiplier and independently
recomputes the eliminated residual series before checking the three
identities.

## Consequence

The m-adic initial ideal contains

```text
(b^3, w^4, t^12).
```

It is therefore `(b,w,t)`-primary. Hence the local F7 incidence algebra is
zero-dimensional and its length is at most

```text
3 * 4 * 12 = 144.
```

In particular, the visible tangent direction is nilpotent: it does not define
a reduced positive-dimensional branch of the special fibre.

## What this does not prove

A zero-dimensional special fibre can still be the reduction of a vertical or
non-flat mixed-characteristic component. Thus this result does **not** exclude
an isolated Z7-point or a point over a ramified extension of Q7. Such a lift
requires a separate mixed-characteristic Hensel or deflation certificate.

The result also says nothing by itself about the existence of the desired
rank-one seed over Q, the rank-17 terminal fibration, or rank at least 30.

## Reproduction

```bash
python3 research/rank17_local_nilpotence_f7.py \
  --compare certificates/rank17_local_nilpotence_f7.json

python3 -m unittest tests.test_rank17_local_nilpotence_f7 -v
```
