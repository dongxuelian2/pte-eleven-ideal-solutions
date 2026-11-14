#!/usr/bin/env python3
"""Exact 13-adic logarithm-rank certificate; dependency-free Python 3.
See PROOF.md Section 5.16 for tiny-integral precision and saturation arguments.
No heights, quadratic-Chabauty roots, or point-completeness sieve are computed.
"""
from fractions import Fraction
P = 13
MOD = P * P
class E:
    def __init__(self, a=0, b=0):
        if isinstance(a, E):
            a, b = a.a, a.b
        self.a, self.b = a % MOD, b % MOD
    def __add__(self, other):
        other = E(other)
        return E(self.a + other.a, self.b + other.b)
    __radd__ = __add__
    def __neg__(self):
        return E(-self.a, -self.b)
    def __sub__(self, other):
        return self + (-E(other))
    def __rsub__(self, other):
        return E(other) + (-self)
    def __mul__(self, other):
        other = E(other)
        return E(self.a * other.a + 2 * self.b * other.b,
                 self.a * other.b + self.b * other.a)
    __rmul__ = __mul__
    def inverse(self):
        norm = (self.a * self.a - 2 * self.b * self.b) % MOD
        inverse_norm = pow(norm, -1, MOD)
        return E(self.a * inverse_norm, -self.b * inverse_norm)
    def __truediv__(self, other):
        return self * E(other).inverse()
    def __bool__(self):
        return bool(self.a or self.b)
    def __eq__(self, other):
        other = E(other)
        return self.a == other.a and self.b == other.b
    def __repr__(self):
        return str((self.a, self.b))
    def residue(self):
        return self.a % P, self.b % P
    def divide_by_p(self):
        assert self.a % P == self.b % P == 0
        return E(self.a // P, self.b // P)

def trim(a):
    a = list(a)
    while a and not a[-1]:
        a.pop()
    return a

def add(a, b):
    return trim([(a[i] if i < len(a) else E()) +
                 (b[i] if i < len(b) else E())
                 for i in range(max(len(a), len(b)))])

def negative(a):
    return [-x for x in a]

def subtract(a, b):
    return add(a, negative(b))

def multiply(a, b):
    c = [E() for _ in range(len(a) + len(b) - 1)] if a and b else []
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)

def scale(a, c):
    return trim([x * c for x in a])

def divide(a, b):
    a, b = trim(a), trim(b)
    q = [E() for _ in range(max(0, len(a) - len(b) + 1))]
    inverse_leading = b[-1].inverse()
    while a and len(a) >= len(b):
        i = len(a) - len(b)
        c = a[-1] * inverse_leading
        q[i] = c
        a = subtract(a, [E()] * i + scale(b, c))
    return trim(q), a

def remainder(a, b):
    return divide(a, b)[1]

def extended_gcd(a, b):
    x, y, xx, yy = [E(1)], [], [], [E(1)]
    while b:
        q, r = divide(a, b)
        a, b = b, r
        x, xx = xx, subtract(x, multiply(q, xx))
        y, yy = yy, subtract(y, multiply(q, yy))
    inverse_leading = a[-1].inverse()
    return scale(a, inverse_leading), scale(x, inverse_leading), scale(y, inverse_leading)

def evaluate(a, x):
    answer = E()
    for c in reversed(a):
        answer = answer * x + c
    return answer

def derivative(a):
    return [a[i] * i for i in range(1, len(a))]

def add_divisors(A, B, f):
    u1, v1 = A
    u2, v2 = B
    assert len(u1) == 3
    d, a, b = extended_gcd(u1, u2)
    assert d == [E(1)]
    assert add(multiply(a, u1), multiply(b, u2)) == [E(1)]
    u = multiply(u1, u2)
    v = remainder(add(v1, multiply(u1, remainder(multiply(a, subtract(v2, v1)), u2))), u)
    assert not remainder(subtract(f, multiply(v, v)), u)
    while len(u) > 3:
        new_u, r = divide(subtract(f, multiply(v, v)), u)
        assert not r
        # A degree drop may discard a point near infinity: reject it.
        assert len(new_u) == 3
        u = scale(new_u, new_u[-1].inverse())
        v = remainder(negative(v), u)
    assert len(u) == 3 and u[-1] == E(1)
    assert not remainder(subtract(f, multiply(v, v)), u)
    return u, v

def translation(points, f):
    (x1, y1), (x2, y2) = points
    assert y1 * y1 == evaluate(f, x1)
    assert y2 * y2 == evaluate(f, x2)
    assert y1.residue() != (0, 0) and y2.residue() != (0, 0)
    slope = (y2 - y1) / (x2 - x1)
    return [x1 * x2, -x1 - x2, E(1)], [y1 - slope * x1, slope]

assert pow(2, 6, 13) == 12
for k in (-15, -3):
    f = [E(1)]
    for factor in ([0, 1], [-k, 1], [5 * k, 1], [-275 * k * k, -68 * k, 1]):
        f = multiply(f, list(map(E, factor)))
    g,ga,gb=extended_gcd(f,derivative(f))
    assert g == [E(1)]
    assert add(multiply(ga,f),multiply(gb,derivative(f))) == [E(1)]
    # Independent counts in F_13 and F_169, not stored Frobenius data.
    squares13 = {a*a % P for a in range(P)}
    squares169 = {(E(a,b)*E(a,b)).residue() for a in range(P) for b in range(P)}
    n1 = 1
    for a in range(P):
        value = evaluate(f,E(a)).residue()
        assert value[1] == 0
        n1 += 1 if value[0] == 0 else (2 if value[0] in squares13 else 0)
    n2 = 1
    for a in range(P):
        for b in range(P):
            value = evaluate(f,E(a,b)).residue()
            n2 += 1 if value == (0,0) else (2 if value in squares169 else 0)
    assert (n1,n2) == (14,182)
    assert (n1*n1+n2)//2-P == 176 and 176 % P != 0
    if k == -15:
        x, y = E(225) / E(4), E(320625) / E(32)
        Q = ([E(81225), E(1020), E(1)], [E(-1564650), E(-13590)])
        points = [(E(6), E(20)), (E(4, 8), E(35, 83))]
        expected = [(8, 7), (7, 12)]
    else:
        x, y = E(33), E(10692)
        Q = ([E(3249), E(60), E(1)], [E(-40014), E(-1026)])
        points = [(E(12, 12), E(127, 49)), (E(1, 7), E(11, 166))]
        expected = [(11, 7), (5, 11)]
    point_class = ([-x, E(1)], [y])
    R = translation(points, f)
    rows = []
    for D in (point_class, Q):
        assert not remainder(subtract(f, multiply(D[1], D[1])), D[0])
        T = R
        for _ in range(176):
            T = add_divisors(T, D, f)
        assert list(map(len, T)) == list(map(len, R))
        assert all(a.residue() == b.residue()
                   for aa, bb in zip(T, R) for a, b in zip(aa, bb))
        first, second = E(), E()
        for xx, yy in points:
            xp = xx - evaluate(T[0], xx) / evaluate(derivative(T[0]), xx)
            yp = evaluate(T[1], xp)
            assert not evaluate(T[0], xp)
            assert yp * yp == evaluate(f, xp) and yp.residue() == yy.residue()
            delta = (xp - xx).divide_by_p()
            first += delta / (2 * yy)
            second += xx * delta / (2 * yy)
        assert first.b % P == second.b % P == 0
        rows.append((first.a % P, second.a % P))
    assert rows == expected
    determinant = (rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]) % P
    assert determinant == 8
    normalized = [tuple(a*pow(176,-1,P)%P for a in row) for row in rows]
    expected_normalized = [(3,1),(1,11)] if k==-15 else [(9,1),(10,9)]
    assert normalized == expected_normalized
    assert (normalized[0][0]*normalized[1][1]-normalized[0][1]*normalized[1][0])%P == 6
    print(k, 'log([176]P)/13, log([176]Q)/13 mod13:', rows)
    print('log(P)/13, log(Q)/13 mod13:', normalized, '; determinant 6')

# Exact rational elliptic arithmetic independently calibrates the third log.
def elliptic_add(U,V):
    if U is None:return V
    if V is None:return U
    x,y=U;u,v=V
    if x==u and y==-v:return None
    slope=(3*x*x-1120*x+72900)/(2*y) if U==V else (v-y)/(u-x)
    w=slope*slope+560-x-u
    return w,-y+slope*(x-w)

def elliptic_multiple(n,U):
    V=None
    while n:
        if n&1:V=elliptic_add(U,V)
        U=elliptic_add(U,U);n//=2
    return V

def valuation13(q):
    q=Fraction(q);a,b=q.numerator,q.denominator;v=0
    assert a
    while a%13==0:a//=13;v+=1
    while b%13==0:b//=13;v-=1
    return v

elliptic_n1 = 1 + sum(1 for a in range(13) for b in range(13)
                     if (b*b-a*(a*a-560*a+72900))%13 == 0)
assert elliptic_n1 == 12
elliptic_point=(Fraction(180),Fraction(900))
xx,yy=elliptic_multiple(12,elliptic_point)
formal_parameter=-xx/yy
assert formal_parameter == Fraction(
    4061199830818098935573944975439863124829948,
    31468337648588340182480022250211888181370801)
assert (valuation13(xx),valuation13(yy),valuation13(formal_parameter))==(-2,-3,1)
q=formal_parameter/13
assert q.numerator*pow(q.denominator,-1,13)%13 == 12
print('Elliptic log(P)/13 = 1 mod13; full product logarithm rank 3.')
print('By the proved prime-to-13 isogenies, both genus-three logarithm ranks are 3.')
print('Known generator subgroups are 13-saturated, not claimed fully saturated.')
print('NOT COMPUTED: heights, QC roots, or a complete rational-point sieve.')
