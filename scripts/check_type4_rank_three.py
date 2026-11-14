#!/usr/bin/env python3
"""Exact rational certificate for the rank-three result in PROOF.md 5.39.

The divisor search used Fisher--Yan's Magma package. This certificate
needs only SymPy; it does not rerun that search or its pairing algorithm.
The rank upper bound is the separate two-descent in Section 5.37.
"""
from math import isqrt
import sympy as s

x = s.symbols('x')
Q = s.Rational
f = 3*x*(240*x*x+40*x+167)*(256*x*x+97*x+466)
h = 184320
q = x*x+Q(1, 6)*x+Q(167, 240)
divisors = [
    (x*x+Q(4243480,243765769)*x+Q(520159249,975063076),
     Q(896286184264773,3805914951397)*x+Q(269660901168378,3805914951397)),
    (x*x+x/4+Q(25,16), Q(855,4)-72*x),
    (x*x+Q(12779,40500)*x+Q(7456,10125),
     Q(1314152,151875)*x+Q(22476112,151875)),
    (x*x+Q(2465261,3168750)*x+Q(1042247,1014000),
     -Q(4410362474,171640625)*x+Q(4413599413,13731250)),
]


def square(a):
    a = Q(a)
    return a >= 0 and isqrt(a.p)**2 == a.p and isqrt(a.q)**2 == a.q


def parity(a, p):
    a = Q(a)
    assert a != 0
    n, d, result = abs(int(a.p)), int(a.q), 0
    for value in (n, d):
        while value % p == 0:
            value //= p
            result ^= 1
    return result


def residue(a, p):
    a = Q(a)
    return int(a.p) % p * pow(int(a.q), -1, p) % p


assert s.LC(s.Poly(f, x)) == h
assert sorted(s.degree(g, x) for g, e in s.factor_list(f, x)[1]) == [1,2,2]
assert s.discriminant(f, x) != 0
assert f.subs(x, Q(1,4)) == 270**2
values = []
for u, v in divisors:
    assert s.degree(u, x) == 2 and s.LC(s.Poly(u, x)) == 1
    assert s.discriminant(u, x) != 0 and s.gcd(u, f) == 1
    assert s.rem(f-v*v, u, x) == 0
    values.append((u.subs(x,0), s.resultant(q,u,x)))
values.append((Q(h,4), h*h*q.subs(x,Q(1,4))))
expected = [(1,32925841), (1,1), (2330,164629205),
            (2505,82479231705), (5,5)]
assert all(square(a/c) and square(b/d)
           for (a,b),(c,d) in zip(values,expected))

# The root alpha=2 mod 17 specifies a split prime of Q(sqrt(-155)).
assert residue(q.subs(x,2),17) == 0
assert residue(s.diff(q,x).subs(x,2),17) != 0
residues = [residue(u.subs(x,2),17) for u,v in divisors]
residues.append(residue(h*(Q(1,4)-2),17))
assert residues == [9,5,9,10,15]
assert all(a != 0 for a in residues)
matrix = s.Matrix([
    [parity(a,2), parity(a,3), parity(a,5), parity(b,13),
     int(pow(r,8,17) == 16)]
    for (a,b),r in zip(values,residues)
])
assert matrix == s.Matrix([[0,0,0,1,0],[0,0,0,0,1],
                          [1,0,1,1,0],[0,1,1,1,1],[0,0,1,0,0]])
assert matrix.det() == 1
print('PASS: four exact rational degree-two divisors and the rational point (1/4,270).')
print('PASS: five independent Kummer classes; rank >= 5 - 2 = 3.')
print('Together with Section 5.37: rank = 3 and Sha[2] = 0.')
print('No saturated Mordell-Weil basis or rational-point classification is asserted.')
