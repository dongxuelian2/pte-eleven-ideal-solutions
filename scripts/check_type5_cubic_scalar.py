"""Exact representative-change counterexample in PROOF.md Section 5.41.

This terminates the fixed cubic-twist shortcut for an even-degree fake
two-cover class; it does not classify the rational points of the curve.
"""
from math import isqrt
import sympy as s

x, lam = s.symbols('x lam')
Q = s.Rational
t = 80*x*x-840*x+589
da = 80*x*x-40*x+4869
db = 80*x*x-1064*x+5125
f = t*da*db
c = Q(49, 4)
assert f.subs(x, c) == 393216**2
assert s.discriminant(f, x) != 0
assert s.factorint(s.discriminant(t, x)) == {2: 10, 5: 1, 101: 1}
r = (105+4*s.sqrt(505))/20
assert s.simplify(t.subs(x, r)) == 0
g = (x-r)*db/80
gamma = s.simplify(g.subs(x, c))
assert s.simplify(gamma-Q(256, 5)*(c-r)) == 0
assert gamma != 0
# Resultants give the relative cubic norm over K=Q(sqrt(505)).
assert s.simplify(s.resultant(g, c-x, x)-gamma) == 0
assert s.simplify(s.resultant(g, 2*(c-x), x)-8*gamma) == 0
assert s.simplify((8*gamma*g.subs(x, c))/(2*gamma)**2) == 2


def rational_square(value):
    value = Q(value)
    return value >= 0 and isqrt(value.p)**2 == value.p and isqrt(value.q)**2 == value.q


assert not rational_square(505)
assert not rational_square(2) and not rational_square(Q(2, 505))
# Odd-degree norm retains the scalar square class, whereas degree four
# makes any rational scalar contribution a square.
assert s.expand(lam**3/(lam**2)) == lam
assert s.expand(lam**4-(lam**2)**2) == 0
print('PASS: the same fake class yields cubic twists with different fibres at c=49/4.')
print('No complete fake-branch exclusion follows from the fixed cubic twist.')
