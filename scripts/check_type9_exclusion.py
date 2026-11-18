"""Exact plane model, global height bounds and split fibers for Type 9 (§5.56)."""

from collections import defaultdict
from fractions import Fraction
from itertools import accumulate, permutations

import sympy as s


R, b, w = s.symbols("R b w")
cubic = (64*R**3 - (352*b + 688*w)*R**2
         + (432*b**2 + 1840*b*w + 2124*w**2)*R
         + 216*b**3 - 492*b**2*w - 1974*b*w**2 - 2025*w**3)

# Generate the singleton polynomial directly from the moments of the chains.
u, v = s.symbols("u v")
powers = []
for k in range(1, 4):
    central = sum(s.Rational(j, 2)**(2*k) for j in (1, 3))
    long_pair = sum((u+j)**(2*k) for j in (-1, 0, 1))
    short_pair = sum((v+s.Rational(j, 2))**(2*k) for j in (-1, 1))
    moment = s.Poly(central + long_pair - short_pair, u, v)
    assert all(i % 2 == j % 2 == 0 for (i, j), _ in moment.terms())
    powers.append(sum(c*((10*b+9)/4)**(i//2)*b**(j//2)
                      for (i, j), c in moment.terms()))
p1, p2, p3 = powers
e2 = (p1**2-p2)/2
e3 = (p3-p1*p2+e2*p1)/3
assert s.expand(64*(R**3-p1*R**2+e2*R-e3) - cubic.subs(w, 1)) == 0

forward = s.Matrix([[-44, -84, 347], [-1088, 3102, 4064], [16, 6, -28]])
back = s.Matrix([[2884, 7, 36757], [-896, 112, 5152], [1456, 28, 5908]])
product = forward * back
assert product[0, 0] != 0 and product == product[0, 0]*s.eye(3)
X, Y, Z = forward * s.Matrix([R, b, w])
weierstrass = Y**2*Z + Y*Z**2 - X**3 - X**2*Z - 295*X*Z**2 - 1109*Z**3
assert s.cancel(weierstrass/cubic).is_Integer and weierstrass.expand() != 0

# Row-sum norms and the coefficient bounds in the Cauchy root estimate.
def row_norm(matrix):
    return max(sum(abs(entry) for entry in matrix.row(i)) for i in range(3))


assert row_norm(forward) < 9000 and row_norm(back) < 40000
coeff_norms = [sum(abs(c) for c in s.Poly(coef, b, w).coeffs())
               for coef in s.Poly(cubic, R).all_coeffs()]
assert coeff_norms == [64, 1040, 4396, 4707]
assert 64 + max(coeff_norms[1:]) < 80*64
assert 1000*40*9000*40000*64*80 < 2**60

# Degree-seven Bezout identities for the duplication polynomials. In each
# block the four coefficients multiply B^3, AB^2, A^2 B, A^3, in that order.
A, B = s.symbols("A B")
F = A**4 - 590*A**2*B**2 - 8874*A*B**3 + 82588*B**4
D = 4*A**3*B + 4*A**2*B**2 + 1180*A*B**3 + 4437*B**4
den = 2079016875
certificates = [
    (B, [18864, 32, 48, 0, 117439, 5900, 4, -12]),
    (A, [1821314552958, 471289250586, -371149600, 2079016875,
         -33900997588392, 3886113560932, 188739889016, 92787400]),
]
basis = [f*A**i*B**(3-i) for f in (F, D) for i in range(4)]
for coordinate, coeffs in certificates:
    assert sum(abs(c) for c in coeffs) < 2**46
    assert s.expand(sum(c*f for c, f in zip(coeffs, basis)) - den*coordinate**7) == 0
assert max(sum(abs(c) for c in s.Poly(f, A, B).coeffs()) for f in (F, D)) < 2**17

# Check x(2P) against the tangent formula from the Weierstrass equation.
x = s.symbols("x")
disc_y = 4*(x**3 + x**2 + 295*x + 1109) + 1
double_x = (3*x**2+2*x+295)**2/disc_y - 1 - 2*x
assert s.cancel(double_x - (F/D).subs({A: x, B: 1})) == 0


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    z, t = right
    if x == z:
        if y+t == -1:
            return None
        slope = (3*x*x+2*x+295)/(2*y+1)
    else:
        slope = (t-y)/(z-x)
    new_x = slope*slope-1-x-z
    return new_x, slope*(x-new_x)-y-1


G = (Fraction(31), Fraction(202))
assert G[1]**2+G[1] == G[0]**3+G[0]**2+295*G[0]+1109
multiples = [None, *accumulate([G]*32, add)]
x32 = multiples[32][0]
assert max(abs(x32.numerator), x32.denominator).bit_length() == 149
assert Fraction(148*3-46, 3*32**2) > Fraction(1, 8)
assert Fraction(3, 4)*26**2/8 > 40 + Fraction(17+46, 3)

# A finite b has W != 0. Group the actual rational projections of all
# 51 multiples in the proven range, including the projective origin.
fibers = defaultdict(list)
for n in range(-25, 26):
    if n == 0:
        point = (0, 1, 0)
    else:
        px, py = multiples[abs(n)]
        point = (px, py if n > 0 else -py-1, 1)
    rr, bb, ww = (sum(int(c)*t for c, t in zip(row, point)) for row in back.tolist())
    if ww:
        fibers[Fraction(bb, ww)].append((n, Fraction(rr, ww)))

assert {value: roots for value, roots in fibers.items() if len(roots) > 1} == {
    0: [(-2, Fraction(25, 4)), (1, Fraction(9, 4))],
    4: [(-3, Fraction(49, 4)), (0, Fraction(1, 4)), (3, Fraction(81, 4))],
}
assert s.expand(cubic.subs({b: 4, w: 1}) - (4*R-1)*(4*R-49)*(4*R-81)) == 0
assert s.expand(cubic.subs({b: 0, w: 1}) - (4*R-9)**2*(4*R-25)) == 0
h2 = 675*b**5-7950*b**4-26125*b**3-35704*b**2-23248*b-6144
assert s.expand(s.discriminant(cubic.subs(w, 1)/64, R) + b*h2/4) == 0
assert 675 % 7 and [h2.subs(b, k) % 7 for k in range(7)] == [2, 1, 6, 1, 3, 2, 6]
branch = b*h2
assert s.degree(s.gcd(branch, branch.diff(b)), b) == 0
assert s.discriminant(cubic.subs({b: 1, w: 0})/64, R) == -s.Rational(675, 4)
other = s.symbols("other")
quadric = s.cancel((cubic - cubic.subs(R, other))/(R-other))
ordered_points = set()
for value, numerators in [(0, (9, 9, 25)), (4, (1, 49, 81))]:
    for first, second, _ in set(permutations(numerators)):
        coords = {R: s.Rational(first, 4), other: s.Rational(second, 4), b: value, w: 1}
        assert cubic.subs(coords) == quadric.subs(coords) == 0
        ordered_points.add((first, second, value))
assert len(ordered_points) == 9
print("The genus-four ordered-root curve has nine rational points.")
print("Type 9: the only finite fiber with three distinct rational roots is b=4.")
print("Its length-two B chain and centered A chain share 3/2.")
