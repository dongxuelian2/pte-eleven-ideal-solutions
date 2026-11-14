"""Exact identities and all final lifting cases for Theorem 17, PROOF 5.43.

Complete elliptic rational-image certification is in the separate Magma
script. This checker verifies the elementary unit argument's exact data
and reconstructs the original rational centers and chains.
"""
from fractions import Fraction
from math import isqrt
import sympy as s

c, u, z = s.symbols('c u z')
Q = s.Rational
S = (20*c+15)/8
T = (160*c*c-20*c+49)/96
U = (69-4*c)/8
V = (2*c-23)*(4*c-31)/48
qT = 6*(160*c*c-20*c+49)
qA = 3*(-80*c*c+1960*c+283)
qB = 3*(-80*c*c+808*c+2875)
qV = 3*(2*c-23)*(4*c-31)
assert s.expand(qT-576*T) == 0
assert s.expand(qA-576*(S*S-4*T)) == 0
assert s.expand(qB-576*(U*U-4*V)) == 0
assert s.expand(qV-144*V) == 0
rt = s.sqrt(15)
r = (245+64*rt)/20
rc = (245-64*rt)/20
eps = 4+rt
assert s.simplify(qA.subs(c, r)) == 0
assert s.simplify((c-r)*(c-rc)+qA/240) == 0
assert s.simplify((3+rt)**2-6*eps) == 0
assert s.simplify(eps*(4-rt)) == 1
assert s.simplify((3+rt)*(3-rt)) == -6
assert s.simplify((Q(49, 4)-r)/(6*(Q(1, 4)-r))-(2/(3+rt))**2) == 0
assert 15 < 4**2  # Minkowski bound sqrt(15)<4.
square5 = {x*x % 5 for x in range(5)}
assert 2 not in square5 and 3 not in square5 and 4 in square5
assert 2 not in {x*x % 3 for x in range(3)}  # No norm -1 unit.
assert qA.subs(c, Q(1, 4)) == 48**2
assert qB.subs(c, Q(1, 4)) == 96**2
assert qV.subs(c, Q(49, 4)) == 9**2


def valuation(n, p):
    n = Fraction(n)
    result = 0
    for sign, v in ((1, abs(n.numerator)), (-1, n.denominator)):
        while v % p == 0:
            result += sign
            v //= p
    return result


assert [valuation(t, 5) for t in (283, 1960, -80)] == [0, 1, 1]
assert s.factorint(abs(s.discriminant(qA, c))) == {2: 18, 3: 3, 5: 1}
assert [valuation(t, 5) for t in (294, -120, 960)] == [0, 1, 1]

# Rebuild both cubic elliptic changes u=L*c, v=L*z.
rr, gam = s.symbols('rr gam')
for q, c0 in ((qB, Q(1, 4)), (qV, Q(49, 4))):
    a, b, d = s.Poly(q, c).all_coeffs()
    LL = a*gam
    xx, yy = LL*c, LL*z
    rhs = (xx**3+gam*(b-a*rr)*xx**2
           +LL*gam*(d-b*rr)*xx-LL**2*gam*d*rr)
    assert s.cancel((rhs-yy**2)/LL**2-(gam*(c-rr)*q-z*z)) == 0


def square(q):
    q = Q(q)
    return q >= 0 and isqrt(q.p)**2 == q.p and isqrt(q.q)**2 == q.q


assert not square(Q(31, 4)) and not square(Q(23, 2))


def centers(total, prod):
    x = s.symbols('x')
    roots = s.solve(x*x-total*x+prod, x)
    assert len(roots) == 2 and all(square(t) for t in roots)
    return [s.sqrt(t) for t in roots]


def paired(centers, length):
    return [sign*center+Q(2*j-length+1, 2)
            for center in centers for sign in (-1, 1) for j in range(length)]


central = [Q(k, 2) for k in range(-5, 6, 2)]
for c0 in (Q(1, 4), Q(49, 4)):
    ac = centers(S.subs(c, c0), T.subs(c, c0))
    bc = centers(U.subs(c, c0), V.subs(c, c0))
    aa = central+paired(ac, 1)
    bb = paired(bc, 1)+paired([s.sqrt(c0)], 3)
    assert len(aa) == len(bb) == 10
    assert len(set(aa)) < 10 or len(set(bb)) < 10 or set(aa) & set(bb)
    if c0 == Q(49, 4):
        assert sorted(aa) == sorted(bb) == [Q(k, 2) for k in range(-9, 10, 2)]
print('PASS: type15 unit/norm identities, both elliptic models and all surviving collisions.')
