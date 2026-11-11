#!/usr/bin/env python3
"""Exact arithmetic certificate for the all-height exclusion in PROOF.md 5.40.

The proof uses two-isogeny descent and good reduction, not a height search.
This checker verifies the birational map, every local obstruction, the
torsion bound, and the two remaining center configurations.
"""
from math import isqrt
import sympy as s

Q = s.Rational
x, v, t = s.symbols('x v t')
f = Q(2, 3)*(v*v-6)*(4*v*v-9)
g = x**3+44*x*x+100*x
assert s.discriminant(f, v) != 0
assert s.expand(f.subs(v, t+Q(3, 2))) == (
    -30*t+14*t*t+16*t**3+Q(8, 3)*t**4)
v_of_x = Q(3, 2)-30/(x+10)
assert s.cancel(f.subs(v, v_of_x)-900*g/(x+10)**4) == 0
assert s.cancel(-10-30/(v_of_x-Q(3, 2))-x) == 0
assert g.subs(x, -10) == 2400


def square(q):
    q = Q(q)
    return q >= 0 and isqrt(q.p)**2 == q.p and isqrt(q.q)**2 == q.q


assert not square(2400) and not square(Q(8, 3))


def chart_residues(A, B, d, p, e):
    """All primitive pairs, up to fourth-power scaling modulo p**e."""
    m = p**e
    assert B % d == 0
    pairs = [(u, 1) for u in range(m)]
    pairs += [(1, w) for w in range(0, m, p)]
    return {(d*u**4+A*u*u*w*w+(B//d)*w**4) % m for u, w in pairs}


certificates = [
    (44, 100, -10, 3, 2, {6, 8}),
    (44, 100, -1, 3, 2, {6, 8}),
    (44, 100, 2, 3, 2, {2, 5, 6}),
    (44, 100, 5, 3, 2, {2, 5, 6}),
    (44, 100, -5, 2, 4, {3, 11, 12}),
    (44, 100, -2, 2, 4, {8, 14}),
    (-88, 1536, 2, 5, 1, {2, 3}),
    (-88, 1536, 3, 5, 1, {2, 3}),
]
for A, B, d, p, e, expected in certificates:
    residues = chart_residues(A, B, d, p, e)
    assert residues == expected
    assert residues.isdisjoint({w*w % (p**e) for w in range(p**e)})
    print(f'A={A}, B={B}, d={d}: excluded modulo {p**e}')


def signed_squarefree_divisors(B):
    return {sign*d for d in s.divisors(B)
            if all(e == 1 for e in s.factorint(d).values())
            for sign in (-1, 1)}


first = signed_squarefree_divisors(100)
first -= {d for A, B, d, *_ in certificates if A == 44}
assert first == {1, 10}
second = signed_squarefree_divisors(1536)
real_excluded = {d for d in second if d < 0}
assert all(d < 0 and -88 < 0 and 1536//d < 0 for d in real_excluded)
second -= real_excluded
second -= {d for A, B, d, *_ in certificates if A == -88}
assert second == {1, 6}
assert 44**2-4*100 == 1536 and square(Q(1536, 6))
assert g.subs(x, 10) == 80**2
# Both descent images have dimension one; both isogeny kernels have order two.
assert len(first)*len(second)//4 == 1

discriminant = 16*100**2*(44**2-4*100)
assert discriminant == 2**17*3*5**4
for p, expected in ((7, 4), (11, 12)):
    assert discriminant % p != 0
    count = 1+sum((yy*yy-xx**3-44*xx*xx-100*xx) % p == 0
                  for xx in range(p) for yy in range(p))
    assert count == expected
    print(f'#E(F_{p})={count}')
assert s.gcd(4, 12) == 4
slope = Q(3*10**2+88*10+100, 2*80)
x_double = slope*slope-44-20
assert x_double == 0 and slope*(10-x_double)-80 == 0
quartic_points = {(Q(3, 2), Q(0))}
for xx, yy in ((0, 0), (10, 80), (10, -80)):
    vv = v_of_x.subs(x, xx)
    y = -Q(30*yy, (xx+10)**2)
    assert f.subs(v, vv) == y*y
    quartic_points.add((vv, y))
assert quartic_points == {(Q(3, 2), 0), (-Q(3, 2), 0), (0, 6), (0, -6)}
assert {vv*vv for vv, _ in quartic_points} == {0, Q(9, 4)}

S = lambda c: (195-20*c)/12
T = lambda c: (160*c*c-1335*c+2304)/36
U = lambda c: 4*(9-c)/3
V = lambda c: 2*(c-6)*(4*c-9)/3
disc_zero = S(Q(0))**2-4*T(Q(0))
assert disc_zero == Q(129, 16) and not square(disc_zero)
c = Q(9, 4)
assert s.expand((x-Q(49, 4))*(x-Q(1, 4))) == x*x-S(c)*x+T(c)
assert U(c) == 9 and V(c) == 0
central_chain = {Q(k, 2) for k in (-5, -3, -1, 1, 3, 5)}
singleton_pair = {-Q(1, 2), Q(1, 2)}
assert singleton_pair <= central_chain
print('TYPE 16 EXCLUDED AT ALL RATIONAL HEIGHTS: two-isogeny descent and collisions')
