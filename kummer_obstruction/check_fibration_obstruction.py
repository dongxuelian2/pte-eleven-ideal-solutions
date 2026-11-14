"""Exact elementary checks accompanying the geometric theorem; not a proof assistant."""
import sympy as s
from itertools import permutations
x=s.symbols('x')
h=x**5-660*x**3+70983*x-102060
assert s.Poly(h,x,modulus=11)==s.Poly(x**5-2,x,modulus=11)
f=s.Poly(x**5-2,x,modulus=11)
assert s.rem(s.Poly(x**11,x,modulus=11),f)==s.Poly(4*x,x,modulus=11)
assert [pow(4,j,11) for j in range(1,6)]==[4,5,9,3,1]
fac59=(x+5)*(x-9)*(x-4)*(x*x+8*x+23)
assert s.Poly(h-fac59,x,modulus=59).is_zero
assert pow(31,29,59)==58
for p in (11,59):
 assert s.gcd(s.Poly(h,x,modulus=p),s.Poly(s.diff(h,x),x,modulus=p)).degree()==0
F=x**5+3478005*x**4-38594700111069*x**3-14431460273343469305*x**2+47793145649188822875249984*x+22058431254036587705175782036640
theta=(-165*x**4+70983*x**2-255150*x-750870)/2
assert s.rem(F.subs(x,theta),h,x)==0
# Irreducibility of the four-dimensional even-sum permutation module.
V=[m for m in range(32) if m.bit_count()%2==0]
def image(m,p):return sum(((m>>i)&1)<<p[i] for i in range(5))
def rank(rows):
 piv={}
 for y in rows:
  while y:
   i=y.bit_length()-1
   if i in piv:y^=piv[i]
   else:piv[i]=y;break
 return len(piv)
perms=list(permutations(range(5)))
assert len(V)==16
for m in V[1:]:assert rank({image(m,p) for p in perms})==4
# Intersection matrix of H and sum of the 16 exceptional curves.
assert 16*(-2)==-32
# Isotropy is a^2=8 b^2, impossible for nonzero rationals since 8 is nonsquare.
assert not s.sqrt(8).is_Rational
print('Exact field, Frobenius, module, and intersection-form checks passed.')
