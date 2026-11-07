#!/usr/bin/env python3
"""Exact elliptic-factor finite local terms and a certified 13-adic height.

Uses Python 3 and SymPy only. See PROOF.md, Section 5.17.
This does not compute genus-two heights or a genus-three QC function.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import sympy as s

A, B = -560, 72900
P = (F(180), F(900))
TWO_TORSION = (F(0), F(0))
N, p = 10, 13
MOD = p**N


def val(q, ell):
    q = F(q)
    if q == 0:
        return float("inf")
    a, b, v = q.numerator, q.denominator, 0
    while a % ell == 0:
        a //= ell
        v += 1
    while b % ell == 0:
        b //= ell
        v -= 1
    return v


def add(R, S):
    if R is None:
        return S
    if S is None:
        return R
    x, y = R
    u, v = S
    if x == u and y == -v:
        return None
    slope = ((v-y)/(u-x) if x != u
             else (3*x*x + 2*A*x + B)/(2*y))
    z = slope*slope - A - x - u
    ans = (z, -y + slope*(x-z))
    assert ans[1]**2 == ans[0]**3 + A*ans[0]**2 + B*ans[0]
    return ans


def denominator(R):
    d = isqrt(R[0].denominator)
    assert d*d == R[0].denominator
    return d


def nonsingular(R, ell):
    if R is None or val(R[0], ell) < 0:
        return True
    x, y = R
    assert val(x, ell) >= 0 and val(y, ell) >= 0
    xx = x.numerator * pow(x.denominator, -1, ell) % ell
    yy = y.numerator * pow(y.denominator, -1, ell) % ell
    return ((2*yy) % ell != 0
            or (3*xx*xx + 2*A*xx + B) % ell != 0)


points = [None]
for n in range(1, 25):
    points.append(add(points[-1], P))
assert add(TWO_TORSION, TWO_TORSION) is None


@lru_cache(None)
def psi(n):
    """Division polynomial evaluated at P, including its sign."""
    x, y = P
    if n == 0:
        return F(0)
    if n == 1:
        return F(1)
    if n == 2:
        return 2*y
    if n == 3:
        return 3*x**4 + 4*A*x**3 + 6*B*x*x - B*B
    if n == 4:
        return 4*y*(x**6 + 2*A*x**5 + 5*B*x**4
                    - 5*B*B*x*x - 2*A*B*B*x - B**3)
    if n % 2:
        m = (n-1)//2
        return psi(m+2)*psi(m)**3 - psi(m-1)*psi(m+1)**3
    m = n//2
    return psi(m)/(2*y)*(psi(m+2)*psi(m-1)**2
                        - psi(m-2)*psi(m+1)**2)


for n in range(2, 13):
    assert P[0] - psi(n-1)*psi(n+1)/psi(n)**2 == points[n][0]

# Exact minimal-model invariants. The 2-adic contradiction is also proved
# symbolically in the note; this exhausts all relevant residue classes.
disc = 16*B*B*(A*A-4*B)
c4 = 16*(A*A-3*B)
assert disc == 2**12 * 3**12 * 5**7 * 11
assert c4 == 1518400
assert all(val(c4, ell) < 4 for ell in [3, 5, 11])
assert not any(
    (A+3*r-ss*ss) % 4 == 0
    and (B+2*A*r+3*r*r-2*ss*tt) % 16 == 0
    for r in range(16) for ss in range(8) for tt in [0, 4, 8, 12]
)

bad = [2, 3, 5, 11]
assert all(nonsingular(points[12], ell) for ell in bad)
assert all(nonsingular(points[24], ell) for ell in bad)
assert val(-points[12][0]/points[12][1], 13) == 1
d12 = denominator(points[12])
assert d12 == 434588679391916602332
assert [val(d12, ell) for ell in bad] == [2, 3, 0, 0]
assert [val(psi(12), ell) for ell in bad] == [74, 123, 90, 0]
local_P = {
    ell: F(val(d12, ell)-val(psi(12), ell), 144) for ell in bad
}
assert local_P == {2: F(-1, 2), 3: F(-5, 6),
                   5: F(-5, 8), 11: F(0)}

# All rational component classes, using E(Q)=Z P + <T> as an existing input.
# kappa is the positive correction subtracted from v_ell(d).
tables = {
    2: [F(0), F(1, 2)],
    3: [F(0), F(5, 6), F(4, 3), F(3, 2), F(4, 3), F(5, 6)],
    5: [F(0), F(5, 8), F(1, 2), F(5, 8)],
}
for ell, table in tables.items():
    order = len(table)
    assert nonsingular(points[order], ell)
    assert all(not nonsingular(points[n], ell) for n in range(1, order))
    for n in range(1, order+1):
        local_n = n*n*local_P[ell] + val(psi(n), ell)
        correction = val(denominator(points[n]), ell) - local_n
        assert correction == table[n % order]
assert nonsingular(add(P, TWO_TORSION), 2)
assert nonsingular(add(points[3], TWO_TORSION), 3)
assert nonsingular(add(points[2], TWO_TORSION), 5)
assert nonsingular(P, 11) and nonsingular(TWO_TORSION, 11)


def residue(q, precision=N):
    q = F(q)
    assert q.denominator % p
    modulus = p**precision
    return q.numerator*pow(q.denominator, -1, modulus) % modulus


def log13(q, precision=N):
    """Standard branch log13(13)=0, with a proved tail bound in the note."""
    q = F(q)
    v = val(q, p)
    q = q/F(p)**v
    assert val(q, p) == 0
    working = p**(precision+2)
    z = (pow(residue(q, precision+2), 12, working)-1) % working
    assert z % p == 0
    total = sum((F((-1)**(k+1)*z**k, 12*k)
                 for k in range(1, precision+3)), F(0))
    return residue(total, precision)


awayP = sum(residue(local_P[ell])*log13(ell) for ell in bad) % MOD
assert awayP == 62652731752
assert awayP % 169 == 65  # In particular, this term is NOT zero.

# Formal sigma for W=span((x-280)omega), not the unit-root splitting.
t = s.symbols("t")
cutoff = 20
u = s.series((1-A*t*t+s.sqrt((1-A*t*t)**2-4*B*t**4))/2,
             t, 0, cutoff+5).removeO()
x_series = u/t**2
w = s.series(1-t*s.diff(u, t)/(2*u), t, 0, cutoff+3).removeO()
rhs = s.series((-x_series+280)*w, t, 0, cutoff+1).removeO()
H = s.integrate(rhs, t)
assert not H.has(s.log)
g = s.integrate(s.series(H*w-1/t, t, 0, cutoff+1).removeO(), t).expand()
expected10 = (-140*t**2 + s.Rational(169525, 3)*t**4
              - s.Rational(310492000, 9)*t**6
              + s.Rational(522868628750, 21)*t**8
              - s.Rational(177457835600000, 9)*t**10)
assert s.series(g, t, 0, 12).removeO() == expected10
gcoeff = [(i, F(g.coeff(t, i))) for i in range(2, cutoff+1)]
L = s.integrate(w, t).expand()
lcoeff = [(i, F(L.coeff(t, i))) for i in range(1, cutoff+1)]


def evaluate(coefficients, q, degree=cutoff):
    return sum((c*q**i for i, c in coefficients if i <= degree), F(0))


def sigma_log(R, degree=cutoff):
    q = -R[0]/R[1]
    assert val(q, 13) >= 1
    return (log13(q)+residue(evaluate(gcoeff, q, degree))) % MOD


def half_height_from_multiple(n):
    R = points[n]
    assert all(nonsingular(R, ell) for ell in bad)
    return (log13(denominator(R))-sigma_log(R))*pow(n*n, -1, MOD) % MOD


half_height = half_height_from_multiple(12)
assert half_height == 29553437420
assert half_height_from_multiple(24) == half_height
assert sigma_log(points[12], 10) == sigma_log(points[12], 20)
height = 2*half_height % MOD
assert height == 59106874840
assert height % 169 == 91

# An independent denominator/correction check with 6P: at 5 it has a
# nontrivial component correction 1/2, which must NOT be dropped.
R6 = points[6]
away6 = (log13(denominator(R6))-residue(F(1, 2))*log13(5)) % MOD
assert (away6-sigma_log(R6))*pow(36, -1, MOD) % MOD == half_height

# Formal elliptic logarithm, useful for the rank-one quadratic coefficient.
q12 = -points[12][0]/points[12][1]
ellP = residue(evaluate(lcoeff, q12)/12)
assert val(ellP, 13) == 1
alpha_times_13 = (height//13)*pow((ellP//13)**2, -1, 13**9) % (13**9)

print("Minimal model and all rational component correction tables checked.")
print("lambda_away(P) = -1/2 log13(2) - 5/6 log13(3) - 5/8 log13(5)")
print(f"lambda_away(P) = {awayP} + O(13^10)")
print(f"h_W(P,P) = {height} + O(13^10)")
print(f"log_omega(P) = {ellP} + O(13^10)")
print(f"13*h_W(P,P)/log_omega(P)^2 = {alpha_times_13} + O(13^9)")
print("Checked the height independently through 6P, 12P, and 24P.")
print("NOT COMPUTED: genus-two heights, full genus-three QC local function or roots.")
