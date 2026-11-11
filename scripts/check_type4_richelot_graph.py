#!/usr/bin/env python3
"""Exact polynomial and Galois-action checks for PROOF.md Section 5.38.

Uses SymPy for rational polynomial arithmetic. No class-group, GRH,
point-search, rank, or rational-point-completeness assumptions enter.
"""
from itertools import combinations, product
from math import isqrt
import sympy as sp

x = sp.symbols('x')
A = [x, 240*x*x+40*x+167, 256*x*x+97*x+466]
B = [128*x*x-233, 240*x*x-167, 13040*x*x+138176*x+2441]
f = 3*sp.prod(A)
g = -6477*sp.prod(B)


def determinant(qs):
    return sp.det(sp.Matrix([[sp.Poly(q, x).nth(i) for i in (2, 1, 0)]
                            for q in qs]))


def bracket(u, v):
    return sp.expand(sp.diff(u, x)*v-u*sp.diff(v, x))


def dual(c, qs):
    return sp.cancel(-c*sp.prod(bracket(qs[i], qs[j])
                     for i, j in [(1, 2), (2, 0), (0, 1)])/determinant(qs))


def squarefree_part(n):
    return sp.prod(p for p, e in sp.factorint(n).items() if e % 2)


def independent(ds):
    for size in range(1, len(ds)+1):
        for subset in combinations(ds, size):
            n = int(sp.prod(subset))
            assert n < 0 or isqrt(n)**2 != n


def matchings(points):
    if not points:
        yield ()
        return
    a = points[0]
    for b in points[1:]:
        for rest in matchings([p for p in points if p not in (a, b)]):
            yield tuple(sorted(((a, b),)+rest))


def action(matching, permutation):
    return tuple(sorted(tuple(sorted((permutation[a], permutation[b])))
                        for a, b in matching))


def group(swaps):
    result = []
    for bits in product((0, 1), repeat=len(swaps)):
        p = list(range(6))
        for bit, (a, b) in zip(bits, swaps):
            if bit:
                p[a], p[b] = p[b], p[a]
        result.append(p)
    return result


assert determinant(A) == -69088
assert determinant(B) == -4773151744
assert sp.cancel(g/dual(3, A)) == 8636**2
assert sp.cancel(dual(-6477, B)/f) == 17272**2
assert all(sp.Poly(h, x, domain=sp.QQ).gcd(
    sp.Poly(sp.diff(h, x), x, domain=sp.QQ)).degree() == 0 for h in (f, g))
ad = [squarefree_part(sp.discriminant(a, x)) for a in A[1:]]
bd = [squarefree_part(sp.discriminant(b, x)) for b in B]
assert ad == [-155, -231]
assert bd == [466, 2505, 32925841]
independent(ad)
independent(bd)
all_matchings = list(matchings(list(range(6))))
assert len(all_matchings) == 15
expected = ((0, 1), (2, 3), (4, 5))
for name, swaps in [('R4', [(2, 3), (4, 5)]),
                    ('dual', [(0, 1), (2, 3), (4, 5)])]:
    G = group(swaps)
    fixed = [m for m in all_matchings if all(action(m, p) == m for p in G)]
    assert fixed == [expected]
    orbit_sizes = sorted({len({action(m, p) for p in G}) for m in all_matchings})
    assert all(n & (n-1) == 0 for n in orbit_sizes)
    print(name, 'Galois order', len(G), 'unique rational matching', fixed[0],
          'matching orbit sizes', orbit_sizes)
print('PASS: exact Richelot models, dual composition, quadratic independence, and unique kernels.')
print('The two-group argument excludes indecomposable double Richelot edges.')
print('No new rank bound or rational-point classification is asserted.')
