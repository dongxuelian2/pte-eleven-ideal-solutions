"""Type 3 moment elimination, squareclasses and the two local tests (§5.58)."""

from itertools import chain, combinations, product
from math import prod

import sympy as s


a, b, c, S, T, w, z = s.symbols("a b c S T w z")
Q = s.Rational


def chain_moment(length, k, center_square):
    offsets = [Q(2*j-length+1, 2) for j in range(length)]
    return sum(s.binomial(2*k, 2*i)*center_square**(k-i)*sum(h**(2*i) for h in offsets)
               for i in range(k+1))


p = [2, S]
for i in range(2, 5):
    p.append(s.expand(S*p[-1]-T*p[-2]))
moments = []
for k in range(1, 5):
    paired = sum(s.binomial(2*k, 2*i)*p[k-i]*2*Q(1, 2)**(2*i) for i in range(k+1))
    moments.append(Q(1, 4)**k + a**k + chain_moment(3, k, b) - c**k - paired)

a0 = 2*S-3*b+c-Q(5, 4)
t0 = -(16*S**2-96*S*b+32*S*c-64*S+96*b**2-48*b*c+156*b-20*c+27)/32
residuals = [s.together(f.subs({a: a0, T: t0}, simultaneous=True)).as_numer_denom()[0].expand()
             for f in moments]
assert residuals[:2] == [0, 0]
H = 1024*b**2-1520*b*c+636*b+640*c**2-1860*c+2025
s0 = (20*b-8*c+27)/6
quotients = []
for f in residuals[2:]:
    assert s.rem(f.subs(S, s0), H, b) == 0
    quotients.append(s.cancel((f-f.subs(S, s0))/(S-s0)))
exception = s.resultant(H, s.resultant(*quotients, S), b)
quadratics = [80*c**2-280*c+281, 144*c**2-616*c+3025, 144*c**2-328*c+4489]
factor = (4*c-25)*(4*c-9)*(4*c-1)**4*prod(quadratics)
assert s.cancel(exception/factor).is_Rational and exception != 0
assert all(s.discriminant(f, c) < 0 for f in quadratics)
for value, roots in [(Q(9, 4), [Q(9, 4), Q(15, 32)]), (Q(25, 4), [Q(25, 4), Q(77, 32)])]:
    assert set(s.solve(H.subs(c, value), b)) == set(roots)
    assert s.sqrt(roots[1]).is_rational is False

# Recover the parameter by the slope of the line through (25/4,25/4).
L = z**2-57*z+1656
A = 49*z**2-393*z+39744
B = 9*z**2+487*z+14904
C = 9*z**2+247*z+39744
V = 2*z**2-14*z+1587
D = 81*z**4-1194*z**3+140521*z**2-1016784*z+43877376
aa, bb, cc = (f/(4*L) for f in (A, B, C))
assert s.cancel(H.subs({b: bb, c: cc})) == 0
assert s.cancel((bb-Q(25, 4))/(cc-Q(25, 4)) - (z-16)/(z-1)) == 0
# The equation for c-25/4 on a nonvertical line is linear after removing the base point.
h = s.symbols('h')
assert s.expand(H.subs({b: Q(25, 4)+w*h, c: Q(25, 4)+h})
                - 16*h*((64*w**2-95*w+40)*h+246*w-210)) == 0
assert s.discriminant(64*w**2-95*w+40, w) < 0
sub = {b: bb, c: cc}
tt = 300*z*V/L**2
assert s.cancel(a0.subs(S, s0).subs(sub)-aa) == 0
assert s.cancel(t0.subs(S, s0).subs(sub)-tt) == 0
assert s.cancel(s0.subs(sub)**2-4*tt-D/L**2) == 0
assert s.limit(tt, z, s.oo) == 0
assert s.discriminant(L, z) < 0 and s.discriminant(V, z) < 0

assert s.factorint(s.resultant(L, A, z)) == {2: 6, 3: 5, 5: 6, 23: 1}
assert s.factorint(V.subs(z, 0)) == {3: 1, 23: 2} and s.LC(s.Poly(V, z)) == 2


def squarefree_products(primes):
    return {prod(subset) for n in range(len(primes)+1) for subset in combinations(primes, n)}


def binary_values(u, v):
    u2, uv, v2 = u*u, u*v, v*v
    return (u2-57*uv+1656*v2, 49*u2-393*uv+39744*v2,
            9*u2+487*uv+14904*v2, 9*u2+247*uv+39744*v2,
            2*u2-14*uv+1587*v2,
            81*u2*u2-1194*u2*uv+140521*u2*v2-1016784*uv*v2+43877376*v2*v2)


assert binary_values(z, 1) == (L, A, B, C, V, D)
classes = set(product(squarefree_products([2, 3, 5, 23]), squarefree_products([2, 3, 23])))
for prime, power, expected in [
    (2, 13, {(1, 6), (3, 46), (115, 46), (345, 6)}),
    (3, 7, {(1, 6), (115, 46)}),
]:
    modulus = prime**power
    squares = {n*n % modulus for n in range(modulus)}
    charts = chain(((u, 1) for u in range(modulus)), ((1, v) for v in range(0, modulus, prime)))
    surviving = set()
    for u, v in charts:
        *first, fifth, disc = binary_values(u, v)
        if disc % modulus not in squares:
            continue
        surviving.update((d, k) for d, k in classes
                         if k*u*v % modulus in squares and 3*k*fifth % modulus in squares
                         and all(d*val % modulus in squares for val in first))
    assert surviving == expected
    classes = surviving

t = s.symbols('t')
q = t**2-Q(131, 8)*t+3381
assert s.cancel(L.subs(z, 1656/t)*t**2/1656-L.subs(z, t)) == 0
assert s.cancel(A.subs(z, 1656/t)*t**2/39744-q) == 0
assert s.sqrt(Q(115*1656, 10)).is_Rational and s.sqrt(Q(115*39744, 15)).is_Rational
assert (aa.subs(z, Q(207, 2)), bb.subs(z, Q(207, 2)), cc.subs(z, Q(207, 2))) == (Q(81, 4), Q(25, 4), Q(25, 4))
conditions = [A/L, B/L, C/L, 3*z*V, D]
assert [s.sqrt(f.subs(z, Q(207, 2))) for f in conditions] == [9, 5, 5, Q(5175, 2), Q(388125, 4)]
assert (A/L).subs(z, 0) == 24
print("Verified the conic parameter, all exceptional fibers, and the two surviving squareclass pairs.")
