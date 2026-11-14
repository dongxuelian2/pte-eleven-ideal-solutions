"""Exact square-center cover and rank-one elliptic projection; requires SymPy."""
import sympy as s
from itertools import combinations
from math import prod
t,a,b,z=s.symbols('t a b z')
P=6*t*t-5*t+5;Q=3*t*t+26*t+10;R=3*t*t-10*t+100;D=3*t*t-10*t+10
A=2*P/D;B=Q/(4*D);C=R/(4*D)
H=24*a*a-80*a*b-100*a+80*b*b+40*b+81
oldD=8*a*a-48*a*b+4*a+48*b*b+24*b-9
N=48*a*a*b-8*a*a-192*a*b*b-168*a*b+50*a+160*b**3+352*b*b+10*b-27
assert s.cancel(H.subs({a:A,b:B}))==0
assert s.cancel((N/(2*oldD)).subs({a:A,b:B})-C)==0
assert s.cancel(C-(-6*A+20*B+21)/8)==0
assert s.expand(4*N-oldD*(-6*a+20*b+21)-(2*a-4*b+1)*H)==0
assert s.cancel(oldD.subs({a:A,b:B})-288*t*t*(t-4)*(3*t-4)/D**2)==0
assert s.cancel(1+t*(B-s.Rational(1,4))-A)==0
for f in (P,Q,R,D):
 assert s.degree(f,t)==2 and s.degree(s.gcd(f,s.diff(f,t)),t)==0
for f,g in combinations((P,Q,R,D),2):assert s.degree(s.gcd(f,g),t)==0
assert s.expand(D.subs(t,(z+5)/3)-(z*z+5)/3)==0
assert s.expand(R.subs(t,(z+5)/3)-(z*z+275)/3)==0
assert s.degree(s.gcd((z*z+5)*(z*z+275),s.diff((z*z+5)*(z*z+275),z)),z)==0
assert s.Rational(225,4)**2==(s.Rational(5,2)**2+5)*(s.Rational(5,2)**2+275)
assert s.Rational(1449,16)**2==(s.Rational(19,4)**2+5)*(s.Rational(19,4)**2+275)

def classes(B):
 ps=list(s.factorint(abs(B)))
 return {sg*prod(I) for k in range(len(ps)+1) for I in combinations(ps,k) for sg in (-1,1)}
def squareclass(n):return (-1 if n<0 else 1)*prod(p for p,k in s.factorint(abs(n)).items() if k%2)
def local(d,A,B,p,k):
 M=p**k;sq={i*i%M for i in range(M)};e=B//d
 return any((d*i**4+A*i*i+e)%M in sq for i in range(M)) or any((d+A*i*i+e*i**4)%M in sq for i in range(0,M,p))
models=[(-560,72900,{1,5},[(180,900)]),(1120,22000,{1,-5,-11,55},[(-1100,0),(-500,12000)])]
counts=[]
for AA,BB,expected,pts in models:
 upper={d for d in classes(BB) if not(d<0 and BB//d<0 and AA<=0) and local(d,AA,BB,2,8) and local(d,AA,BB,3,1)}
 assert upper==expected
 lower={1,squareclass(BB)}
 for x,y in pts:
  assert y*y==x*(x*x+AA*x+BB)
  lower.add(squareclass(x))
 assert lower==upper;counts.append(len(upper))
assert (-2*models[0][0],models[0][0]**2-4*models[0][1])==models[1][:2]
assert counts[0]*counts[1]==8
assert -2*280==-560 and 280**2-4*1375==72900
print('Parametrization, eight branch points, square-center lift and exact elliptic rank one verified.')

q=s.symbols('q');X=2*q+2*z*z+280;Y=2*z*X
assert s.rem(s.expand(Y*Y-X*(X*X-560*X+72900)),s.Poly(q*q-z**4-280*z*z-1375,q),q)==0
print('Birational elliptic map and all exceptional parameters verified.')
