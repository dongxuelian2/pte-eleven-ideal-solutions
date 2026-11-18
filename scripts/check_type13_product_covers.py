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

# The original singleton cubic, including the multiplicity of reduced roots.
G = (144*X**3 - (240*c + 1200)*X**2 + (-80*c*c + 2080*c + 1816)*X - N)/144
b = (25 - 4*c)/6


def reduce_rational(value, p):
    numerator, denominator = map(int, s.Rational(value).as_numer_denom())
    return (numerator * pow(denominator, -1, p)) % p


for p, expected, exclusions in [
    (11, [3, 9], [(9, 6, 10)]),
    (23, [6, 8, 12], [(8, 11, 19), (12, 13, 10)]),
]:
    squares = {j*j % p for j in range(p)}
    assert -2764800 % p not in squares
    admissible = []
    for r in sorted(squares):
        if reduce_rational(b.subs(c, r), p) not in squares:
            continue
        coefficients = s.Poly(G.subs(c, r), X).all_coeffs()
        residual = s.Poly.from_list([reduce_rational(z, p) for z in coefficients], X, modulus=p)
        for root in squares:
            factor = s.Poly(X - root, X, modulus=p)
            while residual.degree() > 0 and residual.eval(root) == 0:
                residual = residual.exquo(factor)
        if residual.degree() == 0:
            admissible.append(r)
    assert admissible == expected
    for r, root, nonsquare in exclusions:
        assert N.subs(c, root) % p == 0 and s.diff(N, c).subs(c, root) % p != 0
        assert reduce_rational((r - root)/(beta - root), p) == nonsquare
        assert nonsquare not in squares
        admissible.remove(r)
    assert admissible == [reduce_rational(beta, p)]

# Dyadic local solvability near the collision, using integral shifted cubics.
def v2(q):
    numerator, denominator = map(int, s.Rational(q).as_numer_denom())
    if not numerator:
        return s.oo
    numerator = abs(numerator)
    return (numerator & -numerator).bit_length() - (denominator & -denominator).bit_length()


h, x = s.symbols('h x')
roots = [beta, s.Rational(9, 4), alpha]
assert s.expand(G.subs(c, beta) - s.prod(X - r for r in roots)) == 0
for r, expected_valuation in zip(roots, [2, 3, 3]):
    shifted = s.Poly(G.subs({c: beta + h, X: r + x}), x)
    for i in range(4):
        coefficient = s.Poly(shifted.nth(i), h)
        bounds = [v2(coefficient.nth(j)) + 11*j for j in range(4)]
        assert min(bounds) >= 0
        if i == 0:
            assert min(bounds) > 2*expected_valuation
            assert min(bounds) - expected_valuation >= 4
        if i == 1:
            assert bounds[0] == expected_valuation < min(bounds[1:])
assert s.expand((beta+h)/beta) == 1+4*h
assert s.expand(b.subs(c, beta+h)/4) == 1-h/6

# The signed sum of the three singleton centers gives a quartic.
t, Z = s.symbols('t Z')
e1, e2, e3 = [(-1)**j*s.Poly(G, X).nth(3-j) for j in [1, 2, 3]]
Q2 = t*t - (10*c + 50*Z)*t/3 + 5*c*c - 30*c*Z + 19*Z*Z
assert s.expand(Q2.subs(Z, 1) - ((t-e1)**2 - 4*e2)) == 0
quartic = s.expand(Q2.subs(Z, 1)**2 - 64*t*e3)
assert list(s.groebner([quartic, s.diff(quartic, c), s.diff(quartic, t)], c, t)) == [1]
infinity = sum(
    coefficient*c**i*t**j for (i, j), coefficient in s.Poly(quartic, c, t).terms() if i+j == 4
)
assert infinity.subs({c: 0, t: 1}) == 1 and s.discriminant(infinity.subs(c, 1), t) != 0

# The Prym determinant over Q(theta) for the beta double cover.
F = -N/240
Ftheta = F.subs(c, theta)
q = s.div(F, c-theta, c)[0]
qhom = s.expand(Z*Z*q.subs(c, c/Z))
delta_beta = beta-theta
form1 = t*(c-theta*Z)
form3 = -s.Rational(320, 3)*qhom
assert s.rem(s.expand(form1*form3 + s.Rational(320, 3)*t*Z**3*F.subs(c, c/Z)), Ftheta, theta) == 0
forms = [form1, Q2, form3]
matrices = [s.hessian(form, (c, t, Z))/2 for form in forms]
M1, M2, M3 = matrices
prym = s.Poly(-(M1 + 2*x*M2 + x*x*M3).det(), x)
reduced = s.Poly.from_list([s.rem(co, Ftheta, theta) for co in prym.all_coeffs()], x)
assert reduced.degree() == 5 and reduced.nth(0) == 0
assert s.rem(s.discriminant(reduced, x), Ftheta, theta) != 0
twist = s.Poly(s.cancel(-delta_beta**4*(M1/delta_beta + 2*x*M2 + x*x*delta_beta*M3).det()
                       - delta_beta*reduced.as_expr().subs(x, delta_beta*x)), x)
assert all(s.rem(co, Ftheta, theta) == 0 for co in twist.all_coeffs())

print("PASS: Type 13 Kummer classes, beta local conditions, signed sum quartic, and Prym determinant.")
