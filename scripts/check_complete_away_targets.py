#!/usr/bin/env python3
"""Algebraic certificate; COMPLETE_AWAY_TARGETS.md supplies geometric inputs."""
import sympy as s
from fractions import Fraction
t,u,y=s.symbols('t u y')
N=80-t*t;D=4*t-30
F=s.expand((N*N+5*D*D)*(N*N+275*D*D))
def val(q,p=3):
 q=Fraction(q)
 if not q:return 999
 a,b=q.numerator,q.denominator;r=0
 while a%p==0:a//=p;r+=1
 while b%p==0:b//=p;r-=1
 return r
def unit(q):
 q=Fraction(q)/Fraction(3)**val(q)
 return q.numerator*pow(q.denominator,-1,3)%3
def square(q):return q==0 or (val(q)%2==0 and unit(q)==1)
def zero_mod(e,m,*v):
 P=s.Poly(s.cancel(e),*v,domain=s.QQ)
 assert all(c.q==1 for c in P.coeffs())
 return all(int(c)%m==0 for c in P.coeffs())
A=4*(20-u*u)**2+5*(4*u-15)**2
B=4*(20-u*u)**2+275*(4*u-15)**2
assert s.expand(F.subs(t,2*u)/16-A*B)==0
assert zero_mod(A*B-3,4,u)
Ni=80*u*u-1;Di=4*u-30*u*u;Si=s.expand(Ni*Ni+140*Di*Di)
Fi=s.expand((Ni*Ni+5*Di*Di)*(Ni*Ni+275*Di*Di))
rhs=-s.Rational(18225,4)*Di**4
model=y*y-Si*y-rhs
assert s.cancel(u**8*F.subs(t,1/u)-Fi)==0
assert s.expand(Fi-Si*Si-4*rhs)==0
assert zero_mod(Si-1,2,u) and zero_mod(rhs,2,u)
assert zero_mod(model-(y*y-y),2,u,y)
assert zero_mod(s.diff(model,y)-1,2,u,y)
assert s.expand(model.subs(y,Si-y)-model)==0
assert all((i*i-i-1)%2 for i in (0,1))
assert val(s.discriminant(F,t))==33
mod=3**8
zs=[next(r for r in range(1,mod,3) if (r*r+n)%mod==0) for n in (5,275)]
assert zs==[5731,3733]
q=lambda z:4*z*z+30*z+80
assert [val(q(z)) for z in zs]==[3,4]
assert [unit(q(z)) for z in zs]==[1,2]
assert [val(q(-z)) for z in zs]==[1,1]
assert [unit(q(-z)) for z in zs]==[1,1]
assert val(zs[0]-zs[1])==3
G=s.Poly(F.subs(t,4+9*u),u)
assert min(val(c) for c in G.all_coeffs())==7
assert zero_mod(G.as_expr()/3**7-(2*u*u+2*u+1),3,u)
assert all((2*i*i+2*i+1)%3 for i in range(3))
assert all(square(F.subs(t,i)) for i in [-1,0,1])
for z,want in [(0,'central'),(16,'bridge'),(17,'A'),(331,'loop')]:
 aa=Fraction(q(z),5);bb=(z*z+5)*(z*z+275)
 assert square(aa) and square(bb)
 va=val(aa)//2;da=min(va,val(-2*z+1));db=min(va,val(-2*z-4))
 got='A' if da>=1 else 'loop' if db>=2 else 'bridge' if db==1 else 'central'
 assert got==want
C=s.Matrix([[1,1,0],[-1,1,0],[0,-2,0],[0,-2,0],[0,0,1]])
gram=C.T*s.diag(4,4,2,2,1)*C
assert gram==s.diag(8,24,1)
dens=C*s.diag(-2,4,-2)*gram.inv()*C.T
assert [dens[i,i] for i in range(5)]==[-s.Rational(1,12)]*2+[s.Rational(2,3)]*2+[-2]
bridge=-s.Rational(1,12)+u/3-u*u/3
long=-s.Rational(1,12)+u*(u-4)/24
loop=-s.Rational(3,4)+u*u-u
assert [bridge.subs(u,i) for i in [0,s.Rational(1,2),s.Rational(3,2),2]]==[-s.Rational(1,12),0,-s.Rational(1,3),-s.Rational(3,4)]
assert loop.subs(u,s.Rational(1,2))==-1
assert 2*s.diff(long,u).subs(u,0)+s.diff(bridge,u).subs(u,0)==0
assert -2*s.diff(bridge,u).subs(u,2)+s.diff(loop,u).subs(u,0)-s.diff(loop,u).subs(u,1)==0
L=s.zeros(5)
for i,j,r in [(0,1,s.Rational(2)),(0,3,s.Rational(1,2)),(3,2,s.Rational(3,2)),(1,4,s.Rational(1,2)),(4,2,s.Rational(3,2))]:
 L[i,i]+=1/r;L[j,j]+=1/r;L[i,j]-=1/r;L[j,i]-=1/r
d=s.Matrix([1,0,0,-1,0]);e=s.Matrix([3,-1,0,1,-3])
phi=s.Matrix([s.Rational(5,24),s.Rational(1,24),-s.Rational(1,8),-s.Rational(1,4),0])
assert L*phi==d and (phi.T*e)[0]==s.Rational(1,3)
assert 2-s.Rational(3,2)-(phi.T*e)[0]==s.Rational(1,6)
assert val(F.subs(t,-1))==2 and val(Fraction(192,17))==1
# log13(3)=log(3^12)/12, evaluated with exact rationals.
a=3**12-1
log=sum((-1)**(n+1)*s.Rational(a**n,n) for n in range(1,7))/12
def residue(q):
 q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,2197)%2197
assert residue(log)==2093
S1=[s.Rational(0),s.Rational(1,6),s.Rational(2,3)]
S5=[s.Rational(0),-s.Rational(1,2),-s.Rational(2,3),s.Rational(4,3)]
assert [residue(-c*log) for c in S1]==[0,1482,1534]
assert [residue(-c*log) for c in S5]==[0,2145,663,871]
print('PASS: dyadic algebra, cluster inputs, local witnesses, graph, calibration, targets')
