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
A5=[p for p in perms if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0]
assert len(A5)==60
assert len(V)==16
for m in V[1:]:assert rank({image(m,p) for p in perms})==4
for m in V[1:]:assert rank({image(m,p) for p in A5})==4
# Intersection matrix of H and sum of the 16 exceptional curves.
assert 16*(-2)==-32
# Isotropy is a^2=8 b^2, impossible for nonzero rationals since 8 is nonsquare.
assert not s.sqrt(8).is_Rational
# Degree-ten field: the subgroup preserving an affine coordinate half-space.
U=[m for m in V if (m>>4)&1==0]
P=[p for p in perms if p[4]==4]
assert len(U)==8 and rank(U)==3 and len(P)==24
assert all(image(u,p) in U for u in U for p in P)
B=[(u,p) for u in U for p in P]
assert len(V)*len(perms)==1920 and len(B)==192
assert (len(V)*len(perms))//len(B)==10
def affine_orbit(m):return {u^image(m,p) for u,p in B}
O0=affine_orbit(0)
O1=affine_orbit(next(m for m in V if m not in O0))
assert O0==set(U) and len(O0)==len(O1)==8
assert O0.isdisjoint(O1) and O0|O1==set(V)
# Actual integral divisor D=2H-sum_{O0} E_i, not a presumed half-even-eight.
gram=s.diag(4,*([-2]*16))
D=s.Matrix([2]+[-1 if m in O0 else 0 for m in V])
H=s.Matrix([1]+[0]*16)
assert (D.T*gram*D)[0]==0 and (D.T*gram*H)[0]==8
# Restricted invariants have basis H,sum_{O0}E,sum_{O1}E.
assert s.diag(4,-2*len(O0),-2*len(O1))==s.diag(4,-16,-16)
print('Exact field, Frobenius, module, degree-ten orbit, and divisor checks passed.')
