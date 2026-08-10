from __future__ import annotations

import argparse
import hashlib
import json
import math
import signal
import traceback
from pathlib import Path

from sage.all import (
    EllipticCurve,
    EllipticCurve_from_cubic,
    GF,
    PolynomialRing,
    QQ,
    matrix,
    prime_range,
    version,
)

# The rational quartic parameterizing complete splitting of the common
# degree-three trisection cover.
QA = QQ(172682352793305664)
QB = QQ(-537299076125816880)
QC = QQ(501757910513324641)
QD = QQ(-138642489541559748)
QE = QQ(-1600508800995452)
S0 = QQ(34) / 19

# A second exact rational point on the quartic, obtained by the tangent
# construction already certified in the repository.
S1 = QQ(-74260562264603258716762922548661741) / QQ(
    7211741447506945055617213747917800
)
Y1_HOM = QQ(
    -2643117350435980352735359399184756839959943134725449336765884451449564992737248
)
S1_DEN = QQ(7211741447506945055617213747917800)
Y1 = Y1_HOM / (S1_DEN**2)

T_NUM = [QQ(1932817518), QQ(-3751492583), QQ(7878955272), QQ(-4386960832)]
T_DEN = [QQ(-1342798922), QQ(-12609842723), QQ(6173506772), QQ(658473248)]

BASE_POINTS = [
    [QQ(-2), QQ(-1), QQ(1)],
    [QQ(-4), QQ(2), QQ(1)],
    [QQ(1), QQ(-3), QQ(1)],
    [QQ(4), QQ(-3), QQ(1)],
    [QQ(2), QQ(-2), QQ(1)],
    [QQ(0), QQ(3), QQ(1)],
    [QQ(3), QQ(4), QQ(1)],
    [QQ(-3), QQ(0), QQ(1)],
    [QQ(5264) / 13547, QQ(-35295) / 13547, QQ(1)],
]

FORMULAS = {
    "C1": {
        "den": [-118794, 154067, 179873, 437604],
        "x": [72912, 289296, -2367708],
        "y": [-356382, 535113, 828915, -1054896],
    },
    "C2": {
        "den": [-17044, 6464, -23081, 13816],
        "x": [46142, 104745, -91352],
        "y": [-51132, 65534, 35502, -49904],
    },
    "Q": {
        "den": [5590868988, -13577828351, 8426176988],
        "x": [-12673653216, 25187117592, -7840260876],
        "y": [-4236741228, 15734911491, -15419156388],
    },
}


class SearchTimeout(Exception):
    pass


def alarm_handler(signum, frame):
    raise SearchTimeout("bounded Sage generator search timed out")


def evaluate(coefficients, x):
    return sum((QQ(c) * x**i for i, c in enumerate(coefficients)), QQ(0))


def quartic(s):
    return QA * s**4 + QB * s**3 + QC * s**2 + QD * s + QE


def normalize_projective(point):
    values = [QQ(v) for v in point]
    if all(v == 0 for v in values):
        raise ValueError("zero projective vector")
    for v in values:
        if v != 0:
            values = [w / v for w in values]
            break
    return tuple(values)


def gradient(F, P):
    X, Y, Z = F.parent().gens()
    return (
        QQ(F.derivative(X)(*P)),
        QQ(F.derivative(Y)(*P)),
        QQ(F.derivative(Z)(*P)),
    )


def collinear(P, Q):
    return (
        P[1] * Q[2] - P[2] * Q[1] == 0
        and P[2] * Q[0] - P[0] * Q[2] == 0
        and P[0] * Q[1] - P[1] * Q[0] == 0
    )


def tangent_direction(F, P):
    gx, gy, gz = gradient(F, P)
    candidates = [(0, gz, -gy), (-gz, 0, gx), (gy, -gx, 0)]
    for D in candidates:
        D = tuple(QQ(v) for v in D)
        if D != (0, 0, 0) and not collinear(P, D):
            return D
    raise ArithmeticError("could not choose tangent direction")


def third_intersection(F, P, Q):
    P = normalize_projective(P)
    Q = normalize_projective(Q)
    lam = PolynomialRing(QQ, "lam").gen()
    if P == Q:
        D = tangent_direction(F, P)
        polynomial = F(*(P[i] + lam * D[i] for i in range(3)))
        coefficients = [QQ(polynomial[i]) for i in range(4)]
        if coefficients[3] == 0:
            raise ArithmeticError("tangent line has unexpected cubic coefficient")
        root = -coefficients[2] / coefficients[3]
        return normalize_projective([P[i] + root * D[i] for i in range(3)])

    polynomial = F(*(P[i] + lam * Q[i] for i in range(3)))
    coefficients = [QQ(polynomial[i]) for i in range(4)]
    # F(P)=F(Q)=0, hence c0=c3=0 and the third point is the other
    # root of lam*(c1+c2*lam).
    if coefficients[2] != 0:
        root = -coefficients[1] / coefficients[2]
        return normalize_projective([P[i] + root * Q[i] for i in range(3)])

    # Swap the affine chart on the line if Q was the repeated point.
    polynomial = F(*(Q[i] + lam * P[i] for i in range(3)))
    coefficients = [QQ(polynomial[i]) for i in range(4)]
    if coefficients[2] == 0:
        raise ArithmeticError("degenerate secant calculation")
    root = -coefficients[1] / coefficients[2]
    return normalize_projective([Q[i] + root * P[i] for i in range(3)])


def cubic_add(F, origin, P, Q):
    if P == origin:
        return normalize_projective(Q)
    if Q == origin:
        return normalize_projective(P)
    R = third_intersection(F, P, Q)
    return third_intersection(F, R, origin)


def cubic_mul(F, origin, n, P):
    if n < 0:
        return cubic_mul(F, origin, -n, third_intersection(F, P, origin))
    result = origin
    addend = P
    while n:
        if n & 1:
            result = cubic_add(F, origin, result, addend)
        addend = cubic_add(F, origin, addend, addend)
        n >>= 1
    return result


def quartic_group_model():
    q0 = quartic(S0)
    if not q0.is_square():
        raise AssertionError("the base quartic value is not a square")
    Y0 = q0.sqrt()
    if Y0 < 0:
        Y0 = -Y0

    # q(S0+u)=q0+q1*u+...+q4*u^4.
    U = PolynomialRing(QQ, "u").gen()
    expansion = (QA * (S0 + U) ** 4 + QB * (S0 + U) ** 3 + QC * (S0 + U) ** 2 + QD * (S0 + U) + QE)
    q = [QQ(expansion[i]) for i in range(5)]

    R = PolynomialRing(QQ, names=("X", "V", "Z"))
    X, V, Z = R.gens()
    C = X * V**2 + 2 * Y0 * V * Z**2 - q[1] * Z**3 - q[2] * X * Z**2 - q[3] * X**2 * Z - q[4] * X**3
    O = normalize_projective((0, q[1] / (2 * Y0), 1))
    u1 = S1 - S0
    if u1 == 0:
        raise AssertionError("second quartic point equals the origin")
    v1 = (Y1 - Y0) / u1
    P1 = normalize_projective((u1, v1, 1))
    if C(*O) != 0 or C(*P1) != 0:
        # The other sign of Y0 may be the one used by the tangent-derived point.
        Y0 = -Y0
        C = X * V**2 + 2 * Y0 * V * Z**2 - q[1] * Z**3 - q[2] * X * Z**2 - q[3] * X**2 * Z - q[4] * X**3
        O = normalize_projective((0, q[1] / (2 * Y0), 1))
        v1 = (Y1 - Y0) / u1
        P1 = normalize_projective((u1, v1, 1))
    if C(*O) != 0 or C(*P1) != 0:
        raise AssertionError("failed to place the certified quartic points on the cubic")
    return C, O, P1, Y0


def point_formula(s, kind):
    formula = FORMULAS[kind]
    denominator = evaluate(formula["den"], s)
    if denominator == 0:
        raise ZeroDivisionError(f"{kind} has a pole")
    return [
        evaluate(formula["x"], s) / denominator,
        evaluate(formula["y"], s) / denominator,
        QQ(1),
    ]


def split_fibre_from_quartic_point(P, Y0):
    X, V, Z = P
    if Z == 0:
        raise ValueError("quartic point maps to infinity")
    u = X / Z
    v = V / Z
    s = S0 + u
    y = Y0 + u * v
    if y**2 != quartic(s):
        raise AssertionError("quartic point failed exact substitution")
    denominator = evaluate(T_DEN, s)
    if denominator == 0:
        raise ZeroDivisionError("cover parameter has a pole")
    t = evaluate(T_NUM, s) / denominator
    R = PolynomialRing(QQ, "r")
    r = R.gen()
    polynomial = sum((T_NUM[i] - t * T_DEN[i]) * r**i for i in range(4))
    factors = polynomial.factor()
    roots = []
    for factor, multiplicity in factors:
        if factor.degree() != 1:
            raise ValueError("the cubic cover is not completely split")
        root = -factor[0] / factor[1]
        roots.extend([QQ(root)] * multiplicity)
    if len(roots) != 3:
        raise ValueError("unexpected number of source roots")
    roots.sort()
    return s, y, t, roots


def build_specialized_curve(t, roots):
    R = PolynomialRing(QQ, names=("X", "Y", "Z"))
    X, Y, Z = R.gens()
    F = (
        3972 * X**3 + 8080 * X**2 * Y - 65622 * X**2 * Z
        + 31679 * X * Y**2 - 104467 * X * Y * Z - 232614 * X * Z**2
        + 24484 * Y**3 - 15556 * Y**2 * Z - 173688 * Y * Z**2
    )
    G = (
        33084 * X**3 + 44912 * X**2 * Y - 62778 * X**2 * Z
        + 24409 * X * Y**2 - 70613 * X * Y * Z - 138714 * X * Z**2
        - 36220 * Y**3 - 122924 * Y**2 * Z + 347376 * Y * Z**2
        + 1042128 * Z**3
    )
    cubic = F + t * G
    origin = BASE_POINTS[8]
    phi = EllipticCurve_from_cubic(cubic, origin, morphism=True)
    E = phi.codomain()
    C = phi.domain()
    base = [phi(C(P)) for P in BASE_POINTS[:8]]
    D1 = []
    D2 = []
    packet = []
    for source in roots:
        P1 = point_formula(source, "C1")
        P2 = point_formula(source, "C2")
        PQ = point_formula(source, "Q")
        if any(cubic(*P) != 0 for P in (P1, P2, PQ)):
            raise AssertionError("trisection point is off the specialized cubic")
        M1, M2, MQ = (phi(C(P)) for P in (P1, P2, PQ))
        D1.append(M1 - MQ)
        D2.append(M2 - MQ)
        packet.append((M1, M2, MQ))
    selected = base + [D1[0], D2[0], D1[1], D2[1]]
    return E, selected, {"base": base, "D1": D1, "D2": D2, "packet": packet}


def rank_mod(rows, ell):
    return 0 if not rows else int(matrix(GF(ell), rows).rank())


def local_quotient(E, points, p, ell):
    Fp = GF(p)
    Ep = EllipticCurve(Fp, [Fp(a.numerator()) / Fp(a.denominator()) for a in E.a_invariants()])
    elements = list(Ep)
    H = {ell * P for P in elements}
    unassigned = set(elements)
    cosets = []
    which = {}
    while unassigned:
        representative = next(iter(unassigned))
        coset = {representative + h for h in H}
        index = len(cosets)
        for P in coset:
            which[P] = index
        unassigned -= coset
        cosets.append(coset)
    zero = which[Ep(0)]
    representatives = [next(iter(C)) for C in cosets]
    coordinates = {zero: ()}
    while len(coordinates) < len(cosets):
        basis = next(i for i in range(len(cosets)) if i not in coordinates)
        new = {}
        for i, vector in coordinates.items():
            for c in range(ell):
                new[which[representatives[i] + c * representatives[basis]]] = vector + (c,)
        coordinates = new
    vectors = []
    for P in points:
        if P.is_zero():
            Q = Ep(0)
        else:
            x, y = P[0], P[1]
            if x.denominator() % p == 0 or y.denominator() % p == 0:
                raise ZeroDivisionError
            Q = Ep(
                Fp(x.numerator()) / Fp(x.denominator()),
                Fp(y.numerator()) / Fp(y.denominator()),
            )
        vectors.append(coordinates[which[Q]])
    return vectors, int(Ep.cardinality()), len(next(iter(coordinates.values())))


def finite_reduction_certificate(E, points, prime_bound=500):
    discriminant = E.discriminant()
    for ell in (2, 3, 5, 7, 11, 13):
        rows = [[] for _ in points]
        records = []
        torsion_witness = None
        for p0 in prime_range(5, prime_bound):
            p = int(p0)
            try:
                if discriminant.numerator() % p == 0 or discriminant.denominator() % p == 0:
                    continue
                vectors, order, dimension = local_quotient(E, points, p, ell)
            except Exception:
                continue
            if torsion_witness is None and order % ell:
                torsion_witness = {"prime": p, "group_order": order}
            if dimension:
                for i, vector in enumerate(vectors):
                    rows[i].extend(vector)
                records.append({"prime": p, "group_order": order, "quotient_dimension": dimension})
                rank = rank_mod(rows, ell)
                if rank == len(points) and torsion_witness is not None:
                    return {
                        "ell": ell,
                        "rank": rank,
                        "rows": rows,
                        "local_records": records,
                        "torsion_witness": torsion_witness,
                    }
    return None


def coords(P):
    return ["0"] if P.is_zero() else [str(P[0]), str(P[1]), str(P[2])]


def curve_height(E):
    values = []
    for a in E.a_invariants():
        values.extend([abs(int(a.numerator())), abs(int(a.denominator()))])
    return max((v.bit_length() for v in values if v), default=0)


def generate_candidates(count):
    cubic, origin, generator, Y0 = quartic_group_model()
    candidates = []
    seen_t = set()
    for n in range(1, count + 1):
        try:
            P = cubic_mul(cubic, origin, n, generator)
            s, y, t, roots = split_fibre_from_quartic_point(P, Y0)
            if t in seen_t:
                continue
            seen_t.add(t)
            E, known, packet = build_specialized_curve(t, roots)
            cert12 = finite_reduction_certificate(E, known, prime_bound=350)
            if cert12 is None or cert12["rank"] != 12:
                continue
            candidates.append({
                "multiple": n,
                "quartic_s": s,
                "quartic_y": y,
                "fibre_parameter": t,
                "cover_roots": roots,
                "curve": E,
                "known": known,
                "packet": packet,
                "rank12_certificate": cert12,
                "coefficient_height_bits": curve_height(E),
            })
        except Exception:
            continue
    candidates.sort(key=lambda item: (item["coefficient_height_bits"], abs(item["multiple"])))
    return candidates


def bounded_sage_generators(E, timeout_seconds):
    old = signal.signal(signal.SIGALRM, alarm_handler)
    signal.alarm(timeout_seconds)
    try:
        rank_estimate = None
        try:
            rank_estimate = E.rank(proof=False)
        except Exception:
            pass
        generators = E.gens(proof=False)
        return rank_estimate, generators, None
    except SearchTimeout as exc:
        return None, [], repr(exc)
    except Exception as exc:
        return None, [], repr(exc)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def search(args):
    candidates = generate_candidates(args.multiples)
    records = []
    winners = []
    for candidate in candidates[: args.curves]:
        E = candidate["curve"]
        known = candidate["known"]
        record = {
            "multiple": candidate["multiple"],
            "quartic_s": str(candidate["quartic_s"]),
            "quartic_y": str(candidate["quartic_y"]),
            "fibre_parameter": str(candidate["fibre_parameter"]),
            "cover_roots": [str(r) for r in candidate["cover_roots"]],
            "a_invariants": [str(a) for a in E.a_invariants()],
            "discriminant": str(E.discriminant()),
            "coefficient_height_bits": candidate["coefficient_height_bits"],
            "known_rank12_labels": [
                "B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7",
                "D1_0", "D2_0", "D1_1", "D2_1",
            ],
            "known_rank12_points": [coords(P) for P in known],
            "rank12_certificate": candidate["rank12_certificate"],
        }
        rank_estimate, generators, error = bounded_sage_generators(E, args.timeout)
        record["sage_rank_estimate_unproved"] = None if rank_estimate is None else int(rank_estimate)
        record["sage_generator_search_error"] = error
        record["sage_generators"] = [coords(P) for P in generators]

        accepted = []
        current = list(known)
        for generator in generators:
            if generator.is_zero() or any(generator == P or generator == -P for P in current):
                continue
            certificate = finite_reduction_certificate(E, current + [generator], prime_bound=args.prime_bound)
            if certificate is not None and certificate["rank"] == len(current) + 1:
                accepted.append({"point": coords(generator), "certificate": certificate})
                current.append(generator)
        record["certified_new_generators"] = accepted
        record["certified_displayed_rank"] = len(current)
        records.append(record)
        if accepted:
            winners.append(record)
            if len(current) >= args.stop_rank:
                break

    payload = {
        "schema_version": 1,
        "status": "completed",
        "sage_version": str(version()),
        "quartic_multiples_requested": args.multiples,
        "rank12_fibres_constructed": len(candidates),
        "curves_deep_searched": len(records),
        "generator_timeout_seconds_per_curve": args.timeout,
        "records": records,
        "winner_count": len(winners),
        "maximum_certified_displayed_rank": max(
            (record["certified_displayed_rank"] for record in records), default=0
        ),
        "solved_rank30": any(record["certified_displayed_rank"] >= 30 for record in records),
        "truth_note": (
            "Sage rank estimates are heuristic search metadata only. A rank increment is reported "
            "only when the exact displayed points have a full finite-reduction certificate."
        ),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["record_sha256"] = hashlib.sha256(raw).hexdigest()
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--multiples", type=int, default=40)
    parser.add_argument("--curves", type=int, default=12)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--prime-bound", type=int, default=700)
    parser.add_argument("--stop-rank", type=int, default=14)
    args = parser.parse_args()
    try:
        payload = search(args)
    except Exception as exc:
        payload = {
            "schema_version": 1,
            "status": "error",
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "solved_rank30": False,
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload.get("status"),
        "rank12_fibres_constructed": payload.get("rank12_fibres_constructed"),
        "curves_deep_searched": payload.get("curves_deep_searched"),
        "winner_count": payload.get("winner_count"),
        "maximum_certified_displayed_rank": payload.get("maximum_certified_displayed_rank"),
        "solved_rank30": payload.get("solved_rank30"),
        "error": payload.get("error"),
        "record_sha256": payload.get("record_sha256"),
    }, indent=2, sort_keys=True))
    if payload.get("status") == "error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
