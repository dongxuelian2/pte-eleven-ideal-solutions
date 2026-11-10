"""Exact Kummer rank lower bounds for the four necessary shift curves.

Only rational components of the odd-degree Kummer map are needed.
No Selmer upper bound, saturation or point completeness is claimed.
"""
from fractions import Fraction
from math import isqrt

import sympy as s

x = s.symbols("x")
Q = s.Rational


def rational_square(a):
    a = Fraction(a)
    return (a >= 0 and isqrt(a.numerator)**2 == a.numerator
            and isqrt(a.denominator)**2 == a.denominator)


def character(a,p):
    a = Fraction(a)
    if p == -1:
        return int(a < 0)
    n,d = abs(a.numerator),a.denominator
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    while d % p == 0:
        d //= p
        out -= 1
    return out % 2


families = {
    5: (x*(80*x*x-840*x+589)*(4*x-49)*(4*x-9),
        [Q(0),Q(49,4),Q(9,4)],
        [(Q(1,4),192),(Q(25,4),1920),(Q(81,4),27648)],
        [[2945,-5,-5],[5,2,2],[5,-2,7],[5,-15,-10],[5,-30,5],[5,10,10]],
        [(0,5),(0,19),(1,-1),(1,2),(1,3),(1,5)]),
    15: (2*x*(2*x-23)*(4*x-31)*(160*x*x-20*x+49),
         [Q(0),Q(23,2),Q(31,4)],
         [(Q(1,4),135),(Q(9,4),1221),(Q(25,4),2205),(Q(2,5),Q(4662,25))],
         [[3565,-115,-310],[115,17871,6],[310,-6,-341],
          [10,-2,-3],[10,-370,-55],[10,-210,-15],[1,-111,-6]],
         [(0,2),(0,5),(0,31),(1,-1),(1,3),(1,5),(2,2)]),
    16: (6*x*(x-6)*(4*x-9)*(160*x*x-1335*x+2304),
         [Q(0),Q(6),Q(9,4)],
         [(Q(3,2),Q(567,2)),(Q(48,5),Q(127008,25))],
         [[15,-10,-15],[10,6,1],[15,-1,-6],[10,-30,-5],[1,6,1]],
         [(2,2),(0,3),(1,5),(0,2),(1,3)]),
}

for typ,(f,roots,points,expected,characters) in families.items():
    h = s.LC(s.Poly(f,x))
    assert s.discriminant(f,x) != 0
    F = s.Poly(s.expand(h**4*f.subs(x,x/h)),x)
    assert F.LC() == 1 and all(c.q == 1 for c in F.all_coeffs())
    # Three rational branch points and one irreducible quadratic give
    # rational two-torsion dimension three.
    factors = s.factor_list(f,x)[1]
    assert sorted(s.degree(g,x) for g,e in factors) == [1,1,1,2]
    assert all(e == 1 for g,e in factors)
    for a,b in points:
        assert f.subs(x,a) == b*b
    values = []
    for a in roots+[p[0] for p in points]:
        values.append([h**3*s.diff(f,x).subs(x,r) if a==r else h*(a-r)
                       for r in roots])
    for actual,reduced in zip(values,expected):
        assert all(rational_square(v/c) for v,c in zip(actual,reduced))
    minor = s.Matrix([[character(row[j],p) for j,p in characters] for row in values])
    assert minor.det() % 2 == 1
    print(f"PASS type {typ}: exact points, Kummer classes and odd determinant; rank >= {len(values)-3}.")

# R4 has two irreducible quadratic factors and only one rational factor.
# That one component already distinguishes two torsion classes and a point.
f = 3*x*(240*x*x+40*x+167)*(256*x*x+97*x+466)
h = s.LC(s.Poly(f,x))
assert h == 184320 and s.discriminant(f,x) != 0
assert sorted(s.degree(g,x) for g,e in s.factor_list(f,x)[1]) == [1,2,2]
assert f.subs(x,Q(1,4)) == 270**2
values = [h**3*s.diff(f,x).subs(x,0), h*h*Q(167,240), h/4]
expected = [1167330,2505,5]
assert all(rational_square(v/c) for v,c in zip(values,expected))
minor = s.Matrix([[character(v,p) for p in [233,167,5]] for v in values])
assert minor.det() % 2 == 1
print("PASS type 4: two branch-divisor classes and one point are independent; rank >= 1.")
print("Classical rank < genus fails directly for types 5, 15 and 16.")
print("No rank upper bound, complete rational-point set or valid PTE output asserted.")
