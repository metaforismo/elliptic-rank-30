#!/usr/bin/env python3
"""Idempotently integrate the 2026-08-07 exact frontier into the living paper."""
from pathlib import Path
import re

PAPER=Path('paper/paper.tex')
MARKER='% BEGIN FRONTIER 2026-08-07'

SECTION=r'''
% BEGIN FRONTIER 2026-08-07
\section{The current structural frontier}

The author of this continuing report is \textbf{Francesco Giannicola}.  The
truth status remains a new intermediate theorem and an improved exact search
method: no curve with thirty certified independent rational points is claimed.

\subsection{A torsion-capacity obstruction}

Let a rational elliptic surface have geometric Mordell--Weil rank $r_0$, and
let $r_1,r_2,r_3$ be the three twist ranks in a connected split
$(\mathbf Z/2\mathbf Z)^2$ base change ramified at three smooth fibres.  The
packet bound is
\[
 r_0+r_1+r_2+r_3\le 4r_0+6.
\]
If the surface has a rational two-torsion section, it admits a minimal model
\[
 y^2=x^3+A(t)x^2+B(t)x,
 \qquad \deg A\le2,\quad \deg B\le4,
\]
with discriminant $16B^2(A^2-4B)$.  A zero of $B$ of multiplicity $m$
contributes root-lattice rank at least $2m-1$.  Since the zero divisor of $B$
has degree four, the total root cost is at least four.  Shioda--Tate gives
$r_0\le4$, and hence
\[
 \boxed{r_0+r_1+r_2+r_3\le22.}
\]
Likewise, Tate normal form for a rational three-torsion section has
$a_3\in H^0(\mathcal O(3))$ and discriminant
$a_3^3(a_1^3-27a_3)$.  Its root cost is at least six, so $r_0\le2$ and
\[
 \boxed{r_0+r_1+r_2+r_3\le14.}
\]
Thus the apparently convenient split $2$- and $3$-descent models cannot carry
a rank-thirty smooth packet.  The high-capacity base should instead have
trivial small rational torsion and, ideally, twelve irreducible $I_1$ fibres.

\subsection{The one-function visibility barrier}

Let $f\in\mathbf Q(E)^\times$ have its only pole at $O$, of exact order $N$.
If its zero divisor consists of rational points $P_1,\ldots,P_N$, counted with
multiplicity, then
\[
 (f)=P_1+\cdots+P_N-NO.
\]
The divisor is principal, so under $\operatorname{Pic}^0(E)\simeq E$,
\[
 \boxed{P_1+\cdots+P_N=O.}
\]
In particular, on $y^2=R(x)$, a near-square construction
\[
 Q(x)^2-R(x)=c\prod_{i=1}^{30}(x-r_i),\qquad \deg Q=15,
\]
produces thirty visible rational points with one forced relation.  It can span
rank at most twenty-nine, regardless of how impressive the point count looks.

The smallest one-function norm target still compatible with rank thirty is
\[
 \boxed{
 A(x)^2-R(x)B(x)^2
 =c\prod_{i=1}^{31}(x-r_i),
 \qquad \deg A\le15,\quad \deg B=14.
 }
\]
The function $A(x)+yB(x)$ has pole order $31$.  Its thirty-one rational zeros
have one unavoidable relation, leaving room---but not a proof---for a
rank-thirty span.  This is now an independent constructive workstream.

\subsection{Exact computational experiments}

The current branch runs several fail-closed experiments: complete
minimal-section schemes for the marked cyclic $j=0$ family; a proved rank and
torsion computation for the auxiliary CM curve $Y^2=X^3-432$; finite-field
function-field rank bounds; deterministic construction of a split $E_8$
rational elliptic surface from a general-position cubic pencil; exact scans of
all three quadratic packet channels; and finite-field seeds for the
thirty-one-zero norm equation.  Characteristic-$p$ output is treated only as
discovery evidence until explicit characteristic-zero equations, sections,
and certificates are produced.
% END FRONTIER 2026-08-07
'''


def main():
    text=PAPER.read_text(encoding='utf-8')
    text=re.sub(r'\\author\{[^{}]*\}',r'\\author{Francesco Giannicola}',text,count=1,flags=re.S)
    if '\\author{' not in text:
        text=text.replace('\\begin{document}','\\author{Francesco Giannicola}\n\\begin{document}',1)
    if MARKER not in text:
        text=text.replace('\\end{document}',SECTION+'\n\\end{document}',1)
    PAPER.write_text(text,encoding='utf-8')

if __name__=='__main__': main()
