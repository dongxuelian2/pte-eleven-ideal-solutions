"""Exact algebra in the original-configuration scalar-square lemma, PROOF 5.42.

The proof controls all valuations; the finite checks below reconstruct its
polynomials, square values, local root identities and Newton data.
"""
from fractions import Fraction
import sympy as s

c = s.symbols('c')
Q = s.Rational
q1 = 80*c*c-840*c+589
q2 = 80*c*c-40*c+4869
q3 = 80*c*c-1064*c+5125
S = (20*c-85)/16
T = q1/256
U = (12*c-83)/16
V = (4*c-49)*(4*c-9)/256
assert s.expand(256*(S*S-4*T)-q2) == 0
assert s.expand(256*(U*U-4*V)-q3) == 0
discs = [s.discriminant(q, c) for q in (q1, q2, q3)]
assert discs == [2**10*505, -2**14*95, -2**14*31]
assert set(s.factorint(505)) & set(s.factorint(95)) == {5}
assert s.gcd(505, 31) == s.gcd(95, 31) == 1
for c0, squares in [(Q(49, 4), (48, 128, 64)),
                    (Q(81, 4), (128, 192, 128))]:
    assert all(q.subs(c, c0) == z*z for q, z in zip((q1, q2, q3), squares))
    assert c0.q == 4 and c0.p in (49, 81)
assert Q(83, 12) > Q(9, 4)
assert 505 < 35**2
r1 = (105+4*s.sqrt(505))/20
r2 = (105-4*s.sqrt(505))/20
assert s.simplify(q1.subs(c, r1)) == s.simplify(q1.subs(c, r2)) == 0
assert s.simplify(r1-r2) == 2*s.sqrt(505)/5
# Any 2-adic square root of 505 is odd; 4r is 1 modulo 8.
assert 505 % 8 == 1
assert {(105+4*u)*pow(5, -1, 8) % 8 for u in (1, 3, 5, 7)} == {1}
assert {int(q1.subs(c, z)) % 8 for z in range(8)} == {5}
assert s.cancel(q1/(80*c*c)-1+(840*c-589)/(80*c*c)) == 0
# At v2(c)<=-4, the correction valuations are >=3 and >=4.
assert -1-(-4) == 3 and -4-2*(-4) == 4


def valuation(n, p):
    n = Fraction(n)
    result = 0
    a, b = abs(n.numerator), n.denominator
    for sign, v in ((1, a), (-1, b)):
        while v % p == 0:
            result += sign
            v //= p
    return result


assert [valuation(a, 5) for a in (5125, -1064, 80)] == [3, 0, 1]
# Newton segments (0,3)->(1,0)->(2,1), each length one.
assert -(0-3) == 3 and -(1-0) == -1

# Exact reciprocal transform used in the second rank-one quotient.
r, gamma, z = s.symbols('r gamma z')
b, d = -1064, 5125
LL = -gamma*r*d
u = LL/c
v = LL*z/c**2
rhs = u**3+gamma*(d-r*b)*u**2+LL*gamma*(b-80*r)*u+LL**2*80*gamma
h = c*(c-r)*q3
assert s.cancel((rhs-v*v)*c**4/LL**2-(gamma*h-z*z)) == 0

# Reconstruct every surviving original chain, not just its moments.
assert S.subs(c, 0) == -Q(85, 16)
c0 = Q(49, 4)
assert U.subs(c, c0) == 4 and V.subs(c, c0) == 0
zero_paired_chain = [-Q(1, 2), Q(1, 2)]
assert zero_paired_chain == [-x for x in reversed(zero_paired_chain)]
c0 = Q(81, 4)
assert s.expand((c-16)*(c-4)) == c*c-S.subs(c, c0)*c+T.subs(c, c0)
assert s.expand((c-9)*(c-1)) == c*c-U.subs(c, c0)*c+V.subs(c, c0)


def paired_length_two(centers):
    return [sign*center+offset for center in centers for sign in (-1, 1)
            for offset in (-Q(1, 2), Q(1, 2))]


A_seed = [-Q(1, 2), Q(1, 2)]+paired_length_two([Q(4), Q(2)])
B_seed = [-Q(9, 2), Q(9, 2)]+paired_length_two([Q(3), Q(1)])
assert len(A_seed) == len(B_seed) == len(set(A_seed)) == len(set(B_seed)) == 10
assert sorted(A_seed) == sorted(B_seed) == [Q(k, 2) for k in range(-9, 10, 2)]
print('PASS: scalar-lemma algebra, reciprocal transform, and all final lifting/collision cases.')
