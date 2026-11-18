"""Exact Kummer and map checks for the Type 13 product covers in PROOF.md."""
from itertools import product
import sympy as s


def squareclass(q):
    numerator, denominator = s.Rational(q).as_numer_denom()
    result = -1 if numerator < 0 else 1
    for prime, exponent in s.factorint(abs(numerator * denominator)).items():
        if exponent % 2:
            result *= prime
    return result


c, theta, delta = s.symbols("c theta delta")
N = -240*c**3 + 1600*c**2 - 260*c + 475
alpha, beta = s.Rational(25, 4), s.Rational(1, 4)
f = 6*c*(25 - 4*c)*N
a = s.Poly(f, c).LC()
assert a == 5760
assert f.subs(c, beta) == 135**2
assert f.subs(c, s.Rational(19, 4)) == s.Rational(2565, 2)**2
assert N.subs(c, 0) and N.subs(c, alpha) and -240 % 31
assert all(N.subs(c, j) % 31 for j in range(31))

# At a rational branch point the vanishing component is replaced by f'_monic.
fp = s.diff(f, c)
components = [
    (a*beta, a*(beta - alpha)),
    (a*s.Rational(19, 4), a*(s.Rational(19, 4) - alpha)),
    (a**3*fp.subs(c, 0), -a*alpha),
    (a*alpha, a**3*fp.subs(c, alpha)),
]
classes = [tuple(map(squareclass, pair)) for pair in components]
assert classes == [(10, -15), (190, -15), (285, -10), (10, -15)]
primes = (2, 3, 19)
characters = [tuple(int(first % p == 0) for p in primes) for first, _ in classes]
by_character = {}
for bits in product((0, 1), repeat=4):
    character = tuple(
        sum(bit*row[i] for bit, row in zip(bits, characters)) % 2
        for i in range(len(primes))
    )
    by_character.setdefault(character, set()).add(bits)
assert by_character[(0, 0, 0)] == {(0, 0, 0, 0), (1, 0, 0, 1)}
assert by_character[characters[0]] == {(1, 0, 0, 0), (0, 0, 0, 1)}
assert N.subs(c, 6) % 11 == 0 and s.diff(N, c).subs(c, 6) % 11 == 5
ratio = ((1 - 4*6)*pow(25 - 4*6, -1, 11)) % 11
assert ratio == 10 and ratio not in {j*j % 11 for j in range(11)}
assert squareclass(-8640) == -15

# The degree-four quotient of each genus-one curve.
d = -6*delta
rhs = d*c*(d*c - d*alpha)*(d*c - d*theta)
lhs = 324*delta**4*c*((25 - 4*c)/6)*((c - theta)/delta)
assert s.cancel(lhs - rhs) == 0
for value, v in [(alpha, 0), (beta, 2)]:
    assert s.sqrt(value).is_Rational
    assert 4*value + 6*v*v == 25

# The alpha-cover coordinate map after the 2-isogeny with kernel (d*theta,0).
X = s.symbols("X")
d = -6*(alpha - theta)
A = d*(2*theta - alpha)
B = d*d*theta*(theta - alpha)
new_a, new_b = -2*A, A*A - 4*B
assert s.cancel(new_b - d*d*alpha*alpha) == 0
dual_x = d*theta + (X - 2*A + new_b/X)/4
c_map = alpha + d*alpha*(alpha - theta)/(dual_x - d*alpha)
assert s.cancel(c_map - alpha*(X*X + new_a*X + new_b)/(X - d*alpha)**2) == 0
print("PASS: Type 13 Kummer independence, two possible classes, and covering maps.")
