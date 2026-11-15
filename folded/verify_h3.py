#!/usr/bin/env python3
"""Exact verification of the explicit H3 semisimplicity calculation.

Requires Python 3 and SymPy. Run:
    python3 verify_h3.py > h3_verification.json

The abstract Coxeter group is W(H3)=A5 x <c>, with c central of order 2.
Its reflections are (a,c), for a a double transposition in A5. The parameter
``t`` is defined by E_i**2=t*E_i; in Liu's folded normalization t=delta**2.
We check the intended Coxeter presentation: the last letters of (2.12) and
(2.13) in arXiv:1305.6528v1 are corrected to r2, consistently with its later
conjugation identity and Coxeter group. No floating-point arithmetic is used.
"""
import itertools
import json
import sympy as sp

POINTS = tuple(range(5))
IDENTITY = POINTS
T = sp.Symbol("t")


def multiply(a, b):
    """Permutation composition a after b."""
    return tuple(a[b[i]] for i in POINTS)


def inverse(a):
    return tuple(a.index(i) for i in POINTS)


def conjugate(g, a):
    return multiply(multiply(g, a), inverse(g))


def parity(a):
    return sum(a[i] > a[j] for i in POINTS for j in range(i + 1, 5)) % 2


def order(a):
    current = IDENTITY
    for n in range(1, 61):
        current = multiply(current, a)
        if current == IDENTITY:
            return n
    raise AssertionError("Permutation order exceeds |A5|.")


def check_relations(r, e):
    """Check every defining family, with both ordered simple-bond cases."""
    count = 0
    identity = sp.eye(r[0].rows)

    def require(lhs, rhs, label):
        nonlocal count
        assert lhs == rhs, f"Failed defining relation: {label}"
        count += 1

    for i in range(3):
        require(r[i] * r[i], identity, f"r{i + 1}^2=1")
        require(r[i] * e[i], e[i], f"r{i + 1}e{i + 1}=e{i + 1}")
        require(e[i] * r[i], e[i], f"e{i + 1}r{i + 1}=e{i + 1}")
        require(e[i] * e[i], T * e[i], f"e{i + 1}^2=t e{i + 1}")
    for i, j in [(0, 2), (2, 0)]:
        require(r[i] * r[j], r[j] * r[i], "nonadjacent r-r")
        require(e[i] * r[j], r[j] * e[i], "nonadjacent e-r")
        require(e[i] * e[j], e[j] * e[i], "nonadjacent e-e")
    for i, j in [(1, 2), (2, 1)]:
        require(r[i] * r[j] * r[i], r[j] * r[i] * r[j], "three-bond braid")
        require(r[j] * r[i] * e[j], e[i] * e[j], "three-bond mixed product")
        require(r[i] * e[j] * r[i], r[j] * e[i] * r[j], "three-bond conjugation")
    i, j = 0, 1
    require(r[i] * r[j] * r[i] * r[j] * r[i],
            r[j] * r[i] * r[j] * r[i] * r[j], "five-bond braid")
    require(r[i] * r[j] * e[i] * r[j] * r[i],
            r[j] * r[i] * e[j] * r[i] * r[j], "five-bond conjugation")
    require(e[i] * e[j] * e[i], e[i], "five-bond e-e-e")
    require(e[i] * r[j] * e[i], e[i], "five-bond length-one sandwich")
    require(e[i] * r[j] * r[i] * r[j] * e[i], e[i],
            "five-bond length-three sandwich")
    return count


def verify():
    group = [p for p in itertools.permutations(POINTS) if parity(p) == 0]
    reflections = [p for p in group if p != IDENTITY and multiply(p, p) == IDENTITY]
    assert len(group) == 60 and len(reflections) == 15
    a = reflections[0]
    klein = [g for g in group if multiply(g, a) == multiply(a, g)]
    assert len(klein) == 4
    b = next(g for g in klein if g not in [IDENTITY, a])
    ab = multiply(a, b)
    representatives = {r: next(g for g in group if conjugate(g, a) == r)
                       for r in reflections}
    middle = next(s for s in reflections
                  if order(multiply(a, s)) == 5 and order(multiply(b, s)) == 3)
    simple = [a, middle, b]
    expected = {
        (1, 1): T**10 * (T - 3)**4 * (T + 12),
        (-1, 1): (T - 4)**6 * (T + 1)**4 * (T + 4)**5,
        (1, -1): (T - 4)**5 * (T - 1)**4 * (T + 4)**6,
        (-1, -1): (T - 4)**5 * (T - 1)**4 * (T + 4)**6,
    }
    results = []
    for x, y in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
        # sigma(s0)=1, sigma(s2)=x, sigma(c)=y; restrict sigma to A5's V4.
        character = {IDENTITY: 1, a: y, b: x * y, ab: x}
        matrix = sp.zeros(15)
        for i, ri in enumerate(reflections):
            for j, rj in enumerate(reflections):
                if i == j:
                    matrix[i, j] = T
                elif multiply(ri, rj) == multiply(rj, ri):
                    matrix[i, j] = 0
                else:
                    choices = [s for s in reflections if conjugate(s, rj) == ri]
                    assert len(choices) == 1
                    h = multiply(multiply(inverse(representatives[ri]), choices[0]),
                                 representatives[rj])
                    assert h in klein
                    matrix[i, j] = y * character[h]
        assert matrix == matrix.T
        determinant = sp.factor(matrix.det(method="domain-ge"))
        assert sp.expand(determinant - expected[x, y]) == 0
        r_matrices, e_matrices = [], []
        for s in simple:
            r_matrix = sp.zeros(15)
            for j, rj in enumerate(reflections):
                ri = conjugate(s, rj)
                i = reflections.index(ri)
                h = multiply(multiply(inverse(representatives[ri]), s), representatives[rj])
                r_matrix[i, j] = y * character[h]
            e_matrix = sp.zeros(15)
            i = reflections.index(s)
            e_matrix[i, :] = matrix[i, :]
            r_matrices.append(r_matrix)
            e_matrices.append(e_matrix)
        count = check_relations(r_matrices, e_matrices)
        for A in r_matrices + e_matrices:
            assert (A.T * matrix - matrix * A).applyfunc(sp.expand) == sp.zeros(15)
        results.append({"dimension": 15, "sigma_s2": x, "sigma_c": y,
                        "determinant": str(determinant),
                        "defining_relation_checks": count})
    r_matrices, e_matrices = [], []
    for s in simple:
        r_matrix = sp.zeros(5)
        for j in POINTS:
            r_matrix[s[j], j] = 1
        fixed = next(i for i in POINTS if s[i] == i)
        e_matrix = sp.zeros(5)
        for j in POINTS:
            e_matrix[fixed, j] = T if j == fixed else 1
        r_matrices.append(r_matrix)
        e_matrices.append(e_matrix)
    count = check_relations(r_matrices, e_matrices)
    small = (T - 1) * sp.eye(5) + sp.ones(5)
    for A in r_matrices + e_matrices:
        assert (A.T * small - small * A).applyfunc(sp.expand) == sp.zeros(5)
    determinant = sp.factor(small.det())
    assert sp.expand(determinant - (T - 1)**4 * (T + 4)) == 0
    results.append({"dimension": 5, "determinant": str(determinant),
                    "defining_relation_checks": count})
    square_sum = 120 + 4 * 15**2 + 5**2
    assert square_sum == 1045
    payload = {
        "status": "all exact assertions passed",
        "parameter": "e_i^2=t e_i; t=delta^2 in Liu's folded normalization",
        "coxeter_generator_A5_permutations_one_based":
            [[i + 1 for i in s] for s in simple],
        "modules": results,
        "total_defining_relation_checks": sum(r["defining_relation_checks"] for r in results),
        "dimension_square_sum": square_sum,
        "contravariance_identities": 30,
        "bad_t": [-12, -4, -1, 0, 1, 3, 4],
        "nonzero_delta_criterion": "delta^2 not in {-12,-4,-1,1,3,4}",
        "scope": "Determinants and relation identities verified here; the proof also uses the published 1045-dimensional basis theorem.",
    }
    return payload


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
