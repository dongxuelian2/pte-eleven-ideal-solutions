"""Exact rank-at-least-two certificate for the necessary type-16 curve.

Uses the rational factors of the odd-degree Kummer map. No class-group,
Selmer upper bound, saturation or rational-point completeness is claimed.
"""
from fractions import Fraction
from math import isqrt

import sympy as s

x = s.symbols("x")
f = 6*x*(x-6)*(4*x-9)*(160*x*x-1335*x+2304)
h = s.LC(s.Poly(f, x))
assert h == 3840 and s.discriminant(f, x) != 0
disc = 1335**2-4*160*2304
assert disc == 307665 and isqrt(disc)**2 != disc
roots = [s.Rational(0), s.Rational(6), s.Rational(9,4)]
points = [(s.Rational(3,2),s.Rational(567,2)),
          (s.Rational(48,5),s.Rational(127008,25))]
for a,b in points:
    assert f.subs(x,a) == b*b

F = s.Poly(s.expand(h**4*f.subs(x,x/h)), x)
assert F.LC() == 1 and all(c.q == 1 for c in F.all_coeffs())
expected = [[15,-10,-15], [10,6,1], [15,-1,-6],
            [10,-30,-5], [1,6,1]]


def rational_square(a):
    a = Fraction(a)
    return (a >= 0 and isqrt(a.numerator)**2 == a.numerator
            and isqrt(a.denominator)**2 == a.denominator)


values = []
for a in roots+[p[0] for p in points]:
    values.append([h**3*s.diff(f,x).subs(x,r) if a==r else h*(a-r)
                   for r in roots])
for actual, reduced in zip(values,expected):
    assert all(rational_square(v/c) for v,c in zip(actual,reduced))
print("PASS: exact points, monic model and all fifteen rational Kummer coordinates.")


def valuation(a,p):
    a = Fraction(a)
    n,d = abs(a.numerator),a.denominator
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    while d % p == 0:
        d //= p
        out -= 1
    return out


# A square-class character is a homomorphism; this minor suffices even
# without evaluating the quadratic factor of the étale algebra.
characters = [(2,2),(0,3),(1,5),(0,2),(1,3)]
minor = s.Matrix([[valuation(row[j],p) % 2 for j,p in characters] for row in values])
assert minor == s.Matrix([[0,1,1,0,0], [0,0,0,1,1], [1,1,0,0,0],
                         [0,0,1,1,1], [0,0,0,0,1]])
assert minor.det() % 2 == 1
assert all(not rational_square(a) for a,_ in points)
print("PASS: five independent Kummer classes; rational two-torsion dimension is three.")
print("Hence rank Jac(R16)(Q) >= 2. Neither displayed abscissa is a center square.")
print("No rank upper bound, rational-point classification or PTE solution asserted.")
