"""Exact local normalization data and final chain collision for PROOF 5.44.

The separate Magma script certifies fake covers, rank, saturation and all
rational images, including the transported inverse minimal-model map.
"""
from fractions import Fraction
from itertools import product
import sympy as s

c, rr, gamma, z = s.symbols('c rr gamma z')
Q = s.Rational
S = (20*c+25)/6
T = (240*c*c+40*c+167)/48
U = (44*c+91)/12
V = (256*c*c+97*c+466)/36
qT = 3*(240*c*c+40*c+167)
qA = -80*c*c+220*c+31
qB = 3*(-720*c*c+2152*c+275)
assert s.expand(qT-144*T) == 0
assert s.expand(qA-9*(S*S-4*T)) == 0
assert s.expand(qB-144*(U*U-4*V)) == 0
assert [s.discriminant(q, c) for q in (qT, qA, qB)] == [
    -2**10*3**2*5*31, 2**4*3**6*5, 2**14*3**2*331]
assert s.discriminant(qT*qA*qB, c) != 0
assert [q.subs(c, Q(1, 4)) for q in (qT, qA, qB)] == [24**2, 9**2, 48**2]


def valuation(n, p):
    n = Fraction(n)
    result = 0
    for sign, v in ((1, abs(n.numerator)), (-1, n.denominator)):
        while v % p == 0:
            result += sign
            v //= p
    return result


assert [valuation(t, 3) for t in (275, 2152, -720)] == [0, 0, 2]
assert [valuation(t, 5) for t in (275, 2152, -720)] == [2, 0, 1]
assert [valuation(t, 3) for t in (501, 120, 720)] == [1, 1, 2]
assert [valuation(t, 5) for t in (501, 120, 720)] == [0, 1, 1]
assert 80 % 3 == 2


def local_square(n, p):
    vp = valuation(n, p)
    if vp % 2:
        return False
    n //= p**vp
    return pow(n % p, (p-1)//2, p) == 1


possible = [2**a*3**b*5**d for a, b, d in product((0, 1), repeat=3)]
assert [n for n in possible if local_square(n, 3) and local_square(n, 5)] == [1]
ss = (1+s.sqrt(5))/2
r = (14+27*ss)/20
assert s.simplify(ss*ss-ss-1) == 0
assert s.simplify(qA.subs(c, r)) == 0
assert qA.subs(c, Q(1, 4)) > 0

# Reconstruct the original cubic-to-elliptic map before minimalization.
a, b, d = 720, 120, 501
LL = a*gamma
xx, yy = LL*c, LL*z
rhs = (xx**3+gamma*(b-a*rr)*xx**2
       +LL*gamma*(d-b*rr)*xx-LL**2*gamma*d*rr)
assert s.cancel((rhs-yy**2)/LL**2-(gamma*(c-rr)*qT-z*z)) == 0

# The only finite c-image repeats entries inside the A seed.
c0 = Q(1, 4)
assert S.subs(c, c0) == 5 and T.subs(c, c0) == 4
assert s.expand((c-1)*(c-4)) == c*c-5*c+4
aa = [-Q(1, 2), Q(1, 2)]
aa += [sign*center+offset for center in (1, 2) for sign in (-1, 1)
       for offset in (-Q(1, 2), Q(1, 2))]
assert len(aa) == 10 and len(set(aa)) == 6
assert aa.count(Q(1, 2)) > 1 and aa.count(Q(3, 2)) > 1
print('PASS: type4 normalization, component-square identities, raw elliptic map and final collision.')
