"""Exact input to the one-prime Chabauty argument for Type 14.

The script computes a rational Mumford representative of 25[P-∞].
The rank upper bound and the analytic argument are in PROOF.md §5.46.
"""
import sympy as s

x = s.symbols("x")
prime = 7
scale = -2700


def poly(expression):
    return s.Poly(expression, x, domain=s.QQ)


original = -3 * (900*x**5 - 7100*x**4 - 21675*x**3 + 71658*x**2 + 3321*x - 93312)
f = poly(scale**4 * original.subs(x, x / s.Integer(scale)))
assert f.LC() == 1
reduction = s.Poly(f.as_expr(), x, modulus=prime)
assert reduction.gcd(reduction.diff()).degree() == 0
affine = {(a, b) for a in range(prime) for b in range(prime)
          if (b*b - f.eval(a)) % prime == 0}
assert affine == {(4, 3), (4, 4)}


def add(D, E):
    u1, v1 = D
    u2, v2 = E
    a, b, d0 = s.gcdex(u1, u2)
    c, e, d = s.gcdex(d0, v1 + v2)
    u = (u1 * u2).exquo(d*d).monic()
    v = (a*c*u1*v2 + b*c*u2*v1 + e*(v1*v2 + f)).exquo(d).rem(u)
    assert (f - v*v).rem(u).is_zero
    while u.degree() > 2:
        u = (f - v*v).exquo(u).monic()
        v = (-v).rem(u)
    assert u.LC() == 1 and v.degree() < u.degree()
    assert (f - v*v).rem(u).is_zero
    return u, v


def multiply(n, divisor):
    result = poly(1), poly(0)
    while n:
        if n % 2:
            result = add(result, divisor)
        n //= 2
        if n:
            divisor = add(divisor, divisor)
    return result


def residue_after_valuation(value):
    numerator, denominator = map(int, s.Rational(value).as_numer_denom())
    assert numerator
    valuation = 0
    while numerator % prime == 0:
        numerator //= prime
        valuation += 1
    while denominator % prime == 0:
        denominator //= prime
        valuation -= 1
    return valuation, numerator * pow(denominator, -1, prime) % prime


assert f.eval(9 * scale) == (3240 * scale**2)**2
u, v = multiply(25, (poly(x - 9 * scale), poly(3240 * scale**2)))
assert [residue_after_valuation(a) for a in u.all_coeffs()] == [(0, 1), (0, 4), (0, 4)]
assert [residue_after_valuation(a) for a in v.all_coeffs()] == [(-1, 5), (-1, 3)]
assert residue_after_valuation(u.discriminant()) == (2, 5)
assert f.eval(5) % prime == 5

print("PASS: smooth mod 7, two affine residue points and the point at infinity.")
print("PASS: exact Cantor addition gives U=(X−5)^2 mod 7, v7(disc U)=2, v7(lc V)=−1.")
