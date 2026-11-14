#!/usr/bin/env python3
"""Exact de Rham/correspondence input certificate for the two genus-3 targets.

Requires Python 3 and SymPy. All calculations are over Q or Q(delta).
This does NOT compute a p-adic logarithm, height, Frobenius matrix, or QC root.
See PROOF.md Section 5.15 for the geometric deductions.
"""
import sympy as s

t, z, T, delta = s.symbols("t z T delta")
Q = s.Rational


def zero(expr):
    assert s.cancel(expr) == 0, expr


def laurent_primitive(expr):
    """Termwise primitive of a finite Laurent polynomial, constant zero."""
    expr = s.expand(expr)
    assert expr.coeff(t, -1) == 0
    out = 0
    for term in s.Add.make_args(expr):
        coefficient, exponent = term.as_coeff_exponent(t)
        assert exponent.is_Integer
        assert exponent != -1
        out += coefficient * t ** (exponent + 1) / (exponent + 1)
    return s.expand(out)


def cup_matrix(coefficients, sheet_scale):
    """<a,b> = sum at infinity Res(F_a b), where dF_a=a.

    The two-sheet scale restores any constant removed from the differentials.
    All our differentials have their only poles at the two infinities.
    """
    n = len(coefficients)
    return s.Matrix(n, n, lambda i, j: s.expand(
        laurent_primitive(coefficients[i]) * coefficients[j]
    ).coeff(t, -1) * sheet_scale)


fE = (z*z + 5) * (z*z + 275)
fC = fE * (4*z*z + 30*z + 80)
reciprocal_C = s.expand(t**6 * fC.subs(z, 1/t) / 4)
SC = s.series(reciprocal_C ** Q(-1, 2), t, 0, 9).removeO()
assert [SC.coeff(t, i) for i in range(4)] == [
    1, Q(-15, 4), Q(-4125, 32), Q(64725, 128)
]
assert s.series(reciprocal_C * SC**2 - 1, t, 0, 9).removeO() == 0

# Divide the C differentials at one infinity by sqrt(delta)/2.
polynomials_C = [s.Integer(1), z, z**3 + Q(15, 4)*z**2,
                 z**4 + Q(4125, 32)*z**2]
coefficients_C = [s.series(-t*P.subs(z, 1/t)*SC, t, 0, 6)
                  .removeO().expand() for P in polynomials_C]
assert [a.coeff(t, -1) for a in coefficients_C] == [0, 0, 0, 0]
CC = cup_matrix(coefficients_C, delta/2)
CC_expected = s.Matrix([
    [0, 0, 0, delta/4],
    [0, 0, delta/2, -45*delta/16],
    [0, -delta/2, 0, 525*delta/2],
    [-delta/4, 45*delta/16, -525*delta/2, 0],
])
assert CC == CC_expected
assert CC + CC.T == s.zeros(4)
zero(CC.det() - delta**4/64)

reciprocal_E = s.expand(t**4 * fE.subs(z, 1/t))
SE = s.series(reciprocal_E ** Q(-1, 2), t, 0, 9).removeO()
coefficients_E = [s.series(-P.subs(z, 1/t)*SE, t, 0, 6)
                  .removeO().expand() for P in [s.Integer(1), z**2]]
assert [a.coeff(t, -1) for a in coefficients_E] == [0, 0]
CE = cup_matrix(coefficients_E, 2)
assert CE == s.Matrix([[0, 2], [-2, 0]])

# Both quotient maps have degree two. Their cross-cup products vanish by
# opposite eigenvalues of the elliptic-quotient involution (proof in note).
CH = s.diag(2*CE, 2*CC)
zero(CH.det() - 4*delta**4)
sigma_H = s.diag(1, 1, -1, -1, -1, -1)
Z_H = s.eye(6) + 3*sigma_H
assert Z_H == s.diag(4, 4, -2, -2, -2, -2)
assert sigma_H*CH*sigma_H.T == CH
assert Z_H*CH == CH*Z_H.T

# Rows express a holomorphic-first symplectic basis in (e0,e1,c0,c1,c2,c3).
B = s.Matrix([
    [1, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, Q(1, 4), 0, 0, 0, 0],
    [0, 0, 0, 1050/delta, 45/(4*delta), 2/delta],
    [0, 0, 0, 0, 1/delta, 0],
])
standard_cup = s.zeros(6)
for i in range(3):
    standard_cup[i, i+3] = 1
    standard_cup[i+3, i] = -1
assert s.simplify(B*CH*B.T) == standard_cup
Z_symplectic = s.diag(4, -2, -2, 4, -2, -2)
assert B*Z_H == Z_symplectic*B
assert sum(Z_symplectic[i, i] for i in range(3)) == 0
assert Z_symplectic*standard_cup == standard_cup*Z_symplectic.T

parameters = {
    1: (80-t*t, 4*t-30, 2*t*t-30*t+160, -2,
        15*t-160, 2*t-15),
    5: (25*t*t-120*t+80, 10*t*t-8, -60*t*t+100*t-48, -20,
        24-25*t, 25-30*t),
}
expected_holomorphic_actions = {
    1: s.Matrix([[Q(45,19), Q(-12,19), Q(4,95)],
                 [Q(480,19), Q(-109,19), Q(6,19)],
                 [Q(5120,19), Q(-960,19), Q(45,19)]]),
    5: s.Matrix([[Q(125,19), Q(-300,19), Q(180,19)],
                 [Q(120,19), Q(-269,19), Q(150,19)],
                 [Q(576,95), Q(-240,19), Q(125,19)]]),
}
expected_discriminants = {
    1: 2**52 * 3**33 * 5**30 * 11**2 * 19**12 * 53,
    5: 2**108 * 3**33 * 5**58 * 11**2 * 19**12 * 53,
}

for d, (N, D, M, lam, sigma_num, sigma_den) in parameters.items():
    R = s.expand(s.diff(N, t)*D - N*s.diff(D, t))
    H = s.expand((N*N + 5*D*D)*(N*N + 275*D*D))
    zero(R - lam*M)
    zero(d*M*M - (4*N*N + 30*N*D + 80*D*D))
    assert s.degree(H, t) == 8
    assert s.discriminant(H, t) == expected_discriminants[d]
    assert s.gcd(s.Poly(H, t, modulus=13),
                 s.Poly(s.diff(H, t), t, modulus=13)).degree() == 0

    u = sigma_num/sigma_den
    v = Q(9025)/sigma_den**4  # sigma(W)=v*W
    zero(u.subs(t, u, simultaneous=True) - t)
    zero(v.subs(t, u, simultaneous=True)*v - 1)
    zero(H.subs(t, u, simultaneous=True) - v*v*H)
    zero((N/D).subs(t, u, simultaneous=True) - N/D)
    zero((M/D).subs(t, u, simultaneous=True) + M/D)
    zero(v/D.subs(t, u, simultaneous=True)**2 - 1/D**2)

    pullback_coefficients = [
        R, R*N*N/D**2,
        lam*D, lam*N,
        lam*(N**3/D**2 + Q(15, 4)*N*N/D),
        lam*(N**4/D**3 + Q(4125, 32)*N*N/D),
    ]
    direct_pullbacks = [R, R*(N/D)**2] + [
        s.diff(N/D, t)*D**3/M*P.subs(z, N/D)
        for P in polynomials_C
    ]
    for i, (h, direct) in enumerate(zip(pullback_coefficients, direct_pullbacks)):
        zero(h - direct)
        sign = 1 if i < 2 else -1
        zero(h.subs(t, u, simultaneous=True)*s.diff(u, t)/v - sign*h)

    # e0,c0,c1 are visibly holomorphic and span H^0(X,Omega^1).
    holomorphic_rows = s.Matrix([
        [s.expand(h).coeff(t, j) for j in range(3)]
        for h in [R, lam*D, lam*N]
    ])
    assert holomorphic_rows.det() != 0
    S = s.Matrix([
        [s.cancel(u**i*s.diff(u, t)/v).expand().coeff(t, j)
         for j in range(3)] for i in range(3)
    ])
    assert S == expected_holomorphic_actions[d]
    assert S*S == s.eye(3)
    assert S.eigenvals() == {s.Integer(1): 1, s.Integer(-1): 2}
    for coefficient in B.subs(delta, d):
        assert s.denom(coefficient) % 13 != 0

    affine_points = [(x, w) for x in range(13) for w in range(13)
                     if (w*w - int(H.subs(t, x))) % 13 == 0]
    assert len(affine_points) == 10
    assert all(w != 0 for _, w in affine_points)
    leading = s.Poly(H, t).LC()
    assert len([w for w in range(13) if (w*w-leading) % 13 == 0]) == 2
    print(f"delta={d}: pullbacks, involution, holomorphic action, good reduction and 12 residue discs checked")

# This checks the exact polynomial projector, NOT a Frobenius matrix.
PE = T*T - 2*T + 13
PC = T**4 + 6*T*T + 169
PH = s.expand(PE*PC)
eE = (2*T+9)*PC/Q(2704)
Zpoly = 6*eE - 2
assert s.rem(eE, PE, T) == 1
assert s.rem(eE, PC, T) == 0
assert s.rem(eE*eE-eE, PH, T) == 0
assert s.rem(Zpoly-4, PE, T) == 0
assert s.rem(Zpoly+2, PC, T) == 0
assert s.rem((Zpoly-4)*(Zpoly+2), PH, T) == 0
assert 2704 == 16*13**2

print("Quotient and genus-3 cup matrices checked; det(cup_X)=4*delta^4")
print("Holomorphic-first symplectic basis and exact trace-zero Z checked")
print("Exact polynomial Frobenius-projector identities checked")
print("NOT COMPUTED: Frobenius matrix, p-adic logarithms, heights, QC roots, or Mordell-Weil sieve")
