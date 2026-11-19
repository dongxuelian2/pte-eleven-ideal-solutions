"""Exact translated Cantor and tiny-integral certificate for Type 11 at 37.

The rank bound is verified by check_type11_genus3_rank_bounds.m.
PROOF.md §5.49 gives the analytic argument and the back-substitution.
"""
from functools import partial

import sympy as s
from sympy.polys import galoistools as gf
from sympy.polys.domains import ZZ

prime = 37
modulus = prime**2
scalar = 58624
x = s.symbols("x")
f = s.Poly(x*(x*x - 11*x - 1)*(7*x*x - 2*x + 18)*(9*x*x - 69*x - 49),
           x, domain=s.QQ)
fp = s.Poly(f.as_expr(), x, modulus=prime)
assert fp.degree() == 7 and fp.gcd(fp.diff()).degree() == 0

add = partial(gf.gf_add, p=modulus, K=ZZ)
subtract = partial(gf.gf_sub, p=modulus, K=ZZ)
multiply = partial(gf.gf_mul, p=modulus, K=ZZ)
scale = partial(gf.gf_mul_ground, p=modulus, K=ZZ)
remainder = partial(gf.gf_rem, p=modulus, K=ZZ)
evaluate = partial(gf.gf_eval, p=modulus, K=ZZ)


def coefficients(polynomial):
    result = []
    for q in polynomial.all_coeffs():
        num, den = map(int, q.as_numer_denom())
        result.append(num * pow(den, -1, modulus) % modulus)
    return gf.gf_trunc(result, modulus)


curve = coefficients(f)


def exact_quotient(a, u):
    q, rem = gf.gf_div(a, u, modulus, ZZ)
    assert not rem
    return q


def inverse(a, u):
    b, _, common = gf.gf_gcdex(a, u, modulus, ZZ)
    assert common == [1] and remainder(multiply(a, b), u) == [1]
    return b


def reduce_divisor(u, v):
    while len(u) > 4:
        if len(u) > 5:
            # A unit top coefficient certifies the degree before reduction.
            assert len(v) == len(u) - 1 and v[0] % prime
        q = exact_quotient(subtract(curve, multiply(v, v)), u)
        u = scale(q, pow(int(q[0]), -1, modulus))
        v = remainder(scale(v, -1), u)
    assert not remainder(subtract(curve, multiply(v, v)), u)
    return u, v


def sum_divisors(first, second):
    u1, v1 = first
    u2, v2 = second
    correction = remainder(multiply(subtract(v2, v1), inverse(u1, u2)), u2)
    return reduce_divisor(multiply(u1, u2), add(v1, multiply(u1, correction)))


def double_divisor(divisor):
    u, v = divisor
    q = exact_quotient(subtract(curve, multiply(v, v)), u)
    correction = remainder(multiply(q, inverse(scale(v, 2), u)), u)
    return reduce_divisor(multiply(u, u), add(v, multiply(u, correction)))


def chart_point(disks):
    points = []
    for a, b in disks:
        difference = int(evaluate(curve, a)) - b*b
        assert difference % prime == 0 and b % prime
        y = b + difference * pow(2*b, -1, modulus)
        points.append(([1, -a % modulus], [y % modulus]))
    return sum_divisors(sum_divisors(points[0], points[1]), points[2])


def log_mod_prime(divisor, disks, expected):
    # Keep B + [n]D in a regular degree-three chart throughout.
    base = chart_point(disks)
    negative_base = (base[0], scale(base[1], -1))
    shifted = sum_divisors(base, divisor)
    for bit in bin(scalar)[3:]:
        shifted = sum_divisors(double_divisor(shifted), negative_base)
        if bit == "1":
            shifted = sum_divisors(shifted, divisor)
    assert shifted == expected
    for a, b in zip(shifted, base):
        assert gf.gf_trunc(a, prime) == gf.gf_trunc(b, prime)

    u, _ = shifted
    derivative = gf.gf_diff(u, modulus, ZZ)
    log = [0, 0, 0]
    for a, b in disks:
        assert evaluate(u, a) % prime == 0
        dx_over_p = -(int(evaluate(u, a)) // prime) * pow(int(evaluate(derivative, a)), -1, prime)
        for j in range(3):
            log[j] += dx_over_p * pow(2*b, -1, prime) * a**j
    return [entry % prime for entry in log]


P = (s.Poly(x - 8, x, domain=s.QQ), s.Poly(1500, x, domain=s.QQ))
D = (s.Poly(x**3 - s.Rational(51, 4)*x*x + s.Rational(23, 4)*x + 8, x, domain=s.QQ),
     s.Poly(s.Rational(415, 8)*x*x + s.Rational(785, 8)*x + 20, x, domain=s.QQ))
for u, v in (P, D):
    assert u.is_monic and u.degree() <= 3 and v.degree() < u.degree()
    assert (f - v*v).rem(u).is_zero
    assert u.gcd(u.diff()).degree() == 0 and u.gcd(f).degree() == 0

log_p = log_mod_prime(tuple(map(coefficients, P)), [(1, -7), (3, 8), (6, 3)],
                      ([1, 693, 101, 1129], [1324, 576, 276]))
log_d = log_mod_prime(tuple(map(coefficients, D)), [(1, 7), (6, -3), (8, -17)],
                      ([1, 540, 876, 1173], [808, 1113, 417]))
assert log_p == [25, 28, 33] and log_d == [3, 30, 12]
annihilator = [int(t) % prime for t in s.Matrix(log_p).cross(s.Matrix(log_d))]
assert annihilator == [12, 21, 0]
assert [sum(t*a**i for i, t in enumerate(annihilator)) % prime for a in (0, 8)] == [12, 32]
assert fp.diff().eval(0) % prime == 31
assert {int(b) for b in range(prime) if (f.eval(0) - b*b) % prime == 0} == {0}
assert {int(b) for b in range(prime) if (f.eval(8) - b*b) % prime == 0} == {17, 20}
print("PASS: independent 37-adic logs and unit differential on the three Type 11 disks.")
