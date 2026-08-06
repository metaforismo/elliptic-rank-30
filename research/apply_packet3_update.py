#!/usr/bin/env python3
"""Idempotently integrate the three-channel packet theorem into the project.

The script updates the living README, STATUS record, evidence/provenance JSON,
and the single maintained LaTeX paper.  It deliberately does not create a
versioned copy of the paper.  The machine certificate is regenerated from the
independent exact verifier.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT_MODULE_PATH = ROOT / "research" / "packet3_certificate.py"


def load_certificate_module():
    spec = importlib.util.spec_from_file_location("packet3_certificate", CERT_MODULE_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def update_status() -> None:
    path = ROOT / "STATUS.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc["truth_status"] = "new intermediate theorem"
    doc["rank_30_curve_found"] = False
    doc["latest_intermediate_result"] = {
        "name": "genus-zero three-channel quadratic rank packet",
        "statement": (
            "For E_t: y^2=x^3-(t^2+1)/(t+1)*x+t and "
            "L=Q(t)(sqrt(t),sqrt(2t/(t+1))), rank E_t(L) is at least "
            "rank E_t(Q(t))+3."
        ),
        "auxiliary_base_genus": 0,
        "forced_nontrivial_character_directions": 3,
        "exact_certificate": "certificates/packet3_certificate.json",
        "human_proof": "research/genus_zero_three_channel_packet.md",
        "record_implication": (
            "Every pair of quadratic multisections must be audited for a hidden "
            "product-twist section; two visible channels can force three independent "
            "directions without increasing the composite-base genus."
        ),
        "record_status": "does not produce a rank-30 curve by itself",
    }
    path.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def update_readme() -> None:
    path = ROOT / "README.md"
    text = path.read_text(encoding="utf-8")
    marker = "## Latest intermediate theorem: a genus-zero three-channel packet"
    if marker in text:
        return
    text += f"""

{marker}

The project now contains an exact family in which **two independent quadratic
characters force all three nontrivial character directions over a genus-zero
base**:

```text
E_t: y^2 = x^3 - (t^2+1)/(t+1) x + t
d1 = t
d2 = 2t/(t+1)
```

The `d1`, `d2`, and `d1*d2` twists contain non-torsion sections at
`x = 0, 1, -1`, respectively. Their total branch support is
`{{0, -1, infinity}}`, so the composite cover is rational and

```text
rank E(L) >= rank E(Q(t)) + 3.
```

This changes the K3 audit target: every pair of quadratic multisections must be
scored by the rank of the **entire character packet**, including the hidden
product twist. See
[`research/genus_zero_three_channel_packet.md`](research/genus_zero_three_channel_packet.md)
and run:

```bash
python3 research/packet3_certificate.py \\
  --output certificates/packet3_certificate.json
```

This is an intermediate theorem, not a rank-30 discovery.
"""
    path.write_text(text, encoding="utf-8")


def update_evidence_index() -> None:
    path = ROOT / "evidence_index.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    entry = {
        "id": "packet3-genus-zero-v1",
        "kind": "exact intermediate theorem certificate",
        "claim": (
            "two independent quadratic characters force three non-torsion "
            "character directions over a genus-zero composite base"
        ),
        "human_proof": "research/genus_zero_three_channel_packet.md",
        "machine_certificate": "certificates/packet3_certificate.json",
        "verifier": "research/packet3_certificate.py",
        "test": "tests/test_packet3_certificate.py",
        "truth_status": "new intermediate theorem",
        "rank_30_claim": False,
    }
    if isinstance(doc, dict):
        bucket = doc.setdefault("intermediate_results", [])
    elif isinstance(doc, list):
        bucket = doc
    else:
        raise TypeError("unsupported evidence_index schema")
    if not any(isinstance(x, dict) and x.get("id") == entry["id"] for x in bucket):
        bucket.append(entry)
    path.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def update_provenance() -> None:
    path = ROOT / "search_provenance.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    entry = {
        "date": "2026-08-06",
        "workstream": "Galois-character height packets",
        "result": "genus-zero three-channel packet",
        "mechanism": (
            "the hidden product twist contributes a third independent character "
            "direction for the same biquadratic base"
        ),
        "next_experiment": (
            "for every pair of quadratic multisections, compute the full product-"
            "twist Mordell-Weil group before scoring or rejecting the pair"
        ),
        "certificate": "certificates/packet3_certificate.json",
    }
    if isinstance(doc, dict):
        bucket = doc.setdefault("research_updates", [])
    elif isinstance(doc, list):
        bucket = doc
    else:
        raise TypeError("unsupported search_provenance schema")
    if not any(isinstance(x, dict) and x.get("result") == entry["result"] for x in bucket):
        bucket.append(entry)
    path.write_text(json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def update_paper() -> None:
    path = ROOT / "paper" / "paper.tex"
    text = path.read_text(encoding="utf-8")
    marker = r"\section{A genus-zero three-channel packet}"
    if marker in text:
        return
    section = r'''
\section{A genus-zero three-channel packet}
\label{sec:packet-three}

The branch-code obstruction limits the number of independent quadratic
\emph{characters} that can be carried by a low-genus composite base, but it
does not limit each character to one Mordell--Weil direction.  The following
explicit construction shows more: two independent characters can force all
three nontrivial eigenspaces while the composite base remains rational.

\begin{theorem}[Three-channel genus-zero packet]
Let $K=\mathbb Q(t)$ and
\[
 E_t:\quad y^2=x^3-\frac{t^2+1}{t+1}x+t.
\]
Put $d_1=t$, $d_2=2t/(t+1)$, and
$L=K(\sqrt{d_1},\sqrt{d_2})$.  Then the connected base curve of $L/K$
has genus zero and
\[
 \operatorname{rank}E_t(L)\geq
 \operatorname{rank}E_t(K)+3.
\]
\end{theorem}

\begin{proof}
Writing $f_t(x)=x^3-(t^2+1)(t+1)^{-1}x+t$, direct substitution gives
\[
 f_t(0)=t,\qquad
 f_t(1)=\frac{2t}{t+1},\qquad
 f_t(-1)=\frac{2t^2}{t+1}.
\]
Hence the $d_1$, $d_2$, and $d_1d_2$ twists have sections with
$x$-coordinates $0$, $1$, and $-1$, respectively.  At the places
$0,-1,\infty$, the branch-parity vectors are
\[
 (1,0,1),\qquad (1,1,0),
\]
so the squareclasses are geometrically independent and the total branch
support has size three.  Riemann--Hurwitz yields
\[
 g=1+2^{2-2}(3-4)=0.
\]

It remains to rule out torsion.  At $t=-3$, the three standard integral twist
models and points are
\[
\begin{array}{c|c|c}
 d & E^{(d)} & P_d\\ \hline
 -3 & y^2=x^3+45x+81 & (0,9)\\
  3 & y^2=x^3+45x-81 & (3,9)\\
 -9 & y^2=x^3+405x+2187 & (9,81).
\end{array}
\]
Exact addition gives
\[
 x(2P_{-3})=\frac{25}{4},\qquad
 x(3P_3)=\frac{1479}{49},\qquad
 x(3P_{-9})=\frac{13077}{121}.
\]
By Lutz--Nagell, each point is non-torsion.  The three generic sections lie in
distinct nontrivial characters of $\operatorname{Gal}(L/K)$; Galois
invariance of the canonical height pairing makes those character spaces
orthogonal.  They are therefore independent.
\end{proof}

The rationality is explicit.  If $u^2=t$, $v^2=2t/(1+t)$ and $w=v/u$, then
$v^2+w^2=2$.  One parametrisation is
\[
 v=\frac{r^2-2r-1}{r^2+1},\qquad
 w=\frac{1-2r-r^2}{r^2+1},\qquad
 u=\frac vw,\qquad t=u^2.
\]
The points $(0,u)$, $(1,v)$, and $(-1,uv)$ are then rational over
$\mathbb Q(r)$.

For the rank-$30$ search, the operational consequence is exact.  A pair of
quadratic multisections must not be credited with only the two sections used
to construct it.  The product twist $E^{d_1d_2}$ is a third channel and can
raise a rank-$17$ rational family directly to rank at least $20$ without a
genus penalty.  The executable certificate is
\texttt{research/packet3\_certificate.py}.

'''
    if r"\end{document}" not in text:
        raise RuntimeError("paper/paper.tex has no end-document marker")
    path.write_text(text.replace(r"\end{document}", section + r"\end{document}", 1), encoding="utf-8")


def write_certificate() -> None:
    module = load_certificate_module()
    result = module.verify_packet()
    path = ROOT / "certificates" / "packet3_certificate.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    write_certificate()
    update_status()
    update_readme()
    update_evidence_index()
    update_provenance()
    update_paper()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
