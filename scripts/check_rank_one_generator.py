"""Exact finite part of the elementary height-saturation proof.

Read PROOF.md Section 5.8 for the uniform real and p-adic
height inequalities and the standard canonical-height theorem inputs.
This does not certify simultaneous lifting or solve the PTE problem.
"""
from fractions import Fraction as Q
from math import gcd, isqrt

A, B = -560, 72900
P, T = (Q(180), Q(900)), (Q(0), Q(0))

def rhs(x):
    return x*(x*x+A*x+B)

def add(U, V):
    if U is None: return V
    if V is None: return U
    x, y = U
    z, w = V
    if x == z and y == -w: return None
    m = (w-y)/(z-x) if x != z else (3*x*x+2*A*x+B)/(2*y)
    xx = m*m-A-x-z
    return xx, -y+m*(x-xx)

def rational_sqrt(q):
    if q < 0: return None
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    return Q(n,d) if n*n == q.numerator and d*d == q.denominator else None

def count_mod(p):
    return 1 + sum(sum(y*y % p == rhs(x) % p for y in range(p)) for x in range(p))

Delta = 16*B*B*(A*A-4*B)
assert Delta == 2**12*3**12*5**7*11
assert Delta % 7 and Delta % 13
assert (count_mod(7), count_mod(13)) == (8,12)
assert gcd(count_mod(7), count_mod(13)) == 4
assert isqrt(A*A-4*B)**2 != A*A-4*B
assert rhs(Q(270)) == Q(270)**2*(-20) < 0
assert rhs(Q(-270)) == Q(270)**2*(-1100) < 0
assert rhs(P[0]) == P[1]**2 and P[1] != 0
assert B == 270**2 and P[0] == 5*6**2

# Archimedean and gcd-bound constants from the written proof.
assert Q(27)*(1-Q(5,6)**2)**2 == Q(121,48)
assert -27*Q(5,6)**2+56*Q(5,6)-27 == Q(11,12)
assert 4*Q(5,6)*Q(11,12) == Q(55,18) > Q(121,48)
assert 4*Q(5,6)*27 == 90 > Q(121,48)
assert 27+56+27 == 110 and 4*110 == 440
assert 2**3*3**3*5*11 == 11880
assert Q(11880)*Q(48,121) == Q(51840,11)
# This strict integer inequality gives H(Q)<24 for every proper generator.
assert 27*440*51840**9 < 24**27*11**9

points = {}
for b in range(1,25):
    for a in range(-24,25):
        if gcd(a,b) != 1: continue
        u, x = Q(a,b), 270*Q(a,b)
        y = rational_sqrt(rhs(x))
        if y is not None: points[u] = (x,y)
expected = {
    Q(0): (Q(0),Q(0)),
    Q(2,3): (Q(180),Q(900)),
    Q(3,2): (Q(405),Q(2025)),
    Q(15,8): (Q(2025,4),Q(38475,8)),
    Q(8,15): (Q(144),Q(1368)),
}
assert points == expected
assert add(P,T) == (Q(405),Q(-2025))
assert add(P,P) == (Q(2025,4),Q(38475,8))
assert add(add(P,P),T) == (Q(144),Q(-1368))
for (x,y) in points.values():
    if y == 0: continue
    u = x/270
    F, G = 27*(u*u-1)**2, 4*u*(27*u*u-56*u+27)
    assert add((x,y),(x,y))[0]/270 == F/G

print('Verified exact torsion data, strict saturation height cutoff, and all rational points below it.')
print('With the proved rank-one descent and the standard height theorem: E(Q) = Z*(180,900) + Z/2*(0,0).')
print('This does not classify the remaining simultaneous lifts or solve the PTE problem.')
