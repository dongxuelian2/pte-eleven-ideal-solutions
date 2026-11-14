#!/usr/bin/env python3
"""Exact arithmetic certificate for the genus-two Richelot target.

No third-party Python packages. This proves the identities and finite counts
used in PROOF.md Section 5.12. It does not compute a Selmer group,
a Mordell-Weil rank upper bound, or a complete rational-point set.
"""
from fractions import Fraction
from functools import reduce
from math import gcd, isqrt

# Polynomial coefficients are in ascending order.
def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a

def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])

def scale(a, c):
    return trim([c * t for t in a])

def mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)

def derivative(a):
    return trim([i * a[i] for i in range(1, len(a))] or [0])

def bracket(a, b):
    return add(mul(derivative(a), b), scale(mul(a, derivative(b)), -1))

def evaluate(a, x):
    y = 0
    for c in reversed(a):
        y = y * x + c
    return y

def determinant3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))

def remainder_mod(a, b, p):
    a = trim([x % p for x in a])
    b = trim([x % p for x in b])
    while a != [0] and len(a) >= len(b):
        c = a[-1] * pow(b[-1], -1, p) % p
        shift = len(a) - len(b)
        for i, t in enumerate(b):
            a[i + shift] = (a[i + shift] - c * t) % p
        a = trim(a)
    return a

def squarefree_mod(a, p):
    a = trim([x % p for x in a])
    b = trim([x % p for x in derivative(a)])
    while b != [0]:
        a, b = b, remainder_mod(a, b, p)
    return len(a) == 1

def legendre(a, p):
    a %= p
    return 0 if a == 0 else 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def counts(a, p):
    """Return #C(F_p), #C(F_{p^2}), a1, a2, #J(F_p).

    For F_{p^2}=F_p[s]/(s^2-d), character(w) = Legendre(Norm(w),p).
    Infinity contributes 1 for odd degree, 1+chi(leading) for even degree.
    """
    assert p > 2 and a[-1] % p and squarefree_mod(a, p)
    n1 = p + sum(legendre(evaluate(a, x), p) for x in range(p))
    n1 += 1 if len(a) % 2 == 0 else 1 + legendre(a[-1], p)
    d = next(t for t in range(2, p) if legendre(t, p) == -1)
    character_sum = 0
    for x in range(p):
        for y in range(p):
            u = v = 0
            for c in reversed(a):
                u, v = (u * x + v * y * d + c) % p, (u * y + v * x) % p
            character_sum += legendre(u * u - d * v * v, p)
    n2 = p * p + character_sum + (1 if len(a) % 2 == 0 else 2)
    a1 = p + 1 - n1
    assert (n2 - p * p - 1 + a1 * a1) % 2 == 0
    a2 = (n2 - p * p - 1 + a1 * a1) // 2
    order = p * p + 1 - a1 * (p + 1) + a2
    return n1, n2, a1, a2, order

def prime_factors(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out

def valuation(n, p):
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k

def torsion_bound(orders):
    """Use only prime-to-p torsion injection, not full p-torsion injection."""
    candidates = set(q for n in orders.values() for q in prime_factors(n))
    bound = 1
    for ell in candidates:
        other_orders = [n for p, n in orders.items() if p != ell]
        assert other_orders
        bound *= ell ** min(valuation(n, ell) for n in other_orders)
    return bound

def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]

def irreducible_monic_quartic(a):
    """Exact Gauss-lemma test, by all possible integer linear/quadratic factors."""
    assert len(a) == 5 and a[-1] == 1 and a[0] != 0
    candidates = [s * d for d in divisors(abs(a[0])) for s in [-1, 1]]
    if any(evaluate(a, r) == 0 for r in candidates):
        return False
    # (x^2+A*x+B)(x^2+C*x+D), with B*D=constant and A+C=cubic.
    for B in candidates:
        D = a[0] // B
        if D != B:
            numerator = a[1] - B * a[3]
            if numerator % (D - B):
                continue
            A = numerator // (D - B)
            C = a[3] - A
            if A * C + B + D == a[2]:
                return False
        elif a[1] == B * a[3]:
            disc = a[3] ** 2 - 4 * a[2] + 8 * B
            if disc >= 0 and isqrt(disc) ** 2 == disc:
                if (a[3] + isqrt(disc)) % 2 == 0:
                    return False
    return True

q = [[5, 0, 1], [80, 30, 4], [185, 36, 1], [275, 0, 1]]
f = [reduce(mul, [q[i] for i in range(4) if i != j], [1]) for j in range(4)]
F = reduce(mul, [[0, 1], [-1, 1], [5, 1], [-275, -68, 1]], [1])
triple = [q[0], q[1], q[3]]
Delta = determinant3(triple)
assert Delta == -8100
L = [bracket(triple[(i + 1) % 3], triple[(i + 2) % 3]) for i in range(3)]
assert L == [[8250, 2040, -30], [0, -540], [-150, 120, 30]]
assert reduce(mul, L, [1]) == scale(F, 486000)
assert scale(F, 486000) == scale(scale(F, -60), Delta)
# Twisting f by delta multiplies its Richelot codomain by delta.
for delta in [1, 5]:
    gt = [scale(q[0], delta), q[1], q[3]]
    Lt = [bracket(gt[(i + 1) % 3], gt[(i + 2) % 3]) for i in range(3)]
    assert determinant3(gt) == delta * Delta
    assert reduce(mul, Lt, [1]) == scale(F, delta * delta * 486000)
# -60 F is 4(-15 F); -300 F is 100(-3 F).
assert scale(F, -60) == scale(scale(F, -15), 4)
assert scale(F, -300) == scale(scale(F, -3), 100)

expected = {
    1: {7: (8, 62, 0, 6, 56), 13: (14, 182, 0, 6, 176),
        17: (22, 334, -4, 30, 392), 23: (26, 586, -2, 30, 608)},
    5: {7: (8, 62, 0, 6, 56), 13: (14, 182, 0, 6, 176),
        17: (14, 334, 4, 30, 248), 23: (22, 586, 2, 30, 512)},
}
for delta, k in [(1, -15), (5, -3)]:
    for p, row in expected[delta].items():
        assert counts(scale(F, k), p) == row
        # delta*y^2=f and Y^2=delta*f are Q-isomorphic.
        assert counts(scale(f[2], delta), p) == row
    assert torsion_bound({p: row[-1] for p, row in expected[delta].items()}) == 8
    print('j=2, delta=', delta, 'counts=', expected[delta])

P1 = (Fraction(-15, 4), Fraction(1425, 32))
P5 = (Fraction(-11), Fraction(1188))
for k, (x, y) in [(-15, P1), (-3, P5)]:
    assert y != 0 and y * y == k * evaluate(F, x)
assert (-68) ** 2 - 4 * (-275) == 36 * 159
assert isqrt(159) ** 2 != 159

# Non-Weierstrass rational controls on every original twist.
for j, ordinate in enumerate([Fraction(3375, 2), Fraction(3375, 8),
                             Fraction(675, 2), Fraction(675, 2)]):
    assert ordinate != 0 and 5 * ordinate * ordinate == evaluate(f[j], Fraction(5, 2))
    assert isqrt(f[j][-1]) ** 2 == f[j][-1]  # two rational infinities for delta=1

# Original delta=5 finite controls and delta=1 rational infinities.
for j, y in enumerate([Fraction(3375,2),Fraction(3375,8),Fraction(675,2),Fraction(675,2)]):
    assert y != 0 and 5*y*y == evaluate(f[j],Fraction(5,2))
    assert isqrt(f[j][-1])**2 == f[j][-1]
# Additional certificates: rational torsion bounds and Q-simplicity for all j.
primes = {0: [7, 13, 17], 1: [7, 13, 53], 2: [7, 13, 17], 3: [7, 11, 47]}
for j in range(4):
    for delta in [1, 5]:
        orders = {p: counts(scale(f[j], delta), p)[-1] for p in primes[j]}
        bound = torsion_bound(orders)
        assert bound == [4, 4, 8, 4][j]
        print('j=', j, 'delta=', delta, 'orders=', orders, 'torsion bound=', bound)
    p = 11 if j == 3 else 7
    n1, n2, a1, a2, order = counts(f[j], p)
    charpoly = [p * p, -p * a1, a2, -a1, 1]
    assert irreducible_monic_quartic(charpoly)
    _,_,twist_a1,twist_a2,_=counts(scale(f[j],5),p)
    assert [p*p,-p*twist_a1,twist_a2,-twist_a1,1]==charpoly
    print('j=', j, 'Q-simplicity prime=', p, 'Frobenius coefficients=', charpoly)
print('PASS: Richelot identities, exact finite-field counts, torsion bounds,')
print('Q-simplicity tests, and infinite-order point coordinates verified.')
print('NOT COMPUTED: rank upper bounds, Selmer groups, saturation, complete rational points.')
