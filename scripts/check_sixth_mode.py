"""Exact sixth-mode identities and a nonexceptional rejection; needs SymPy."""
import sympy as s
from itertools import combinations
a,b,c,u,v=s.symbols('a b c u v')
def offsets(L):return [s.Rational(2*j-L+1,2) for j in range(L)]
def moment(L,z,q,k):
 P=s.Poly(sum((z+d)**(2*k) for d in offsets(L)),z)
 assert all(mo[0]%2==0 for mo,co in P.terms())
 return sum(co*q**(mo[0]//2) for mo,co in P.terms())
p=[None]+[sum(d**(2*k) for d in offsets(4))/2+moment(2,u,a,k)+c**k-moment(3,v,b,k) for k in range(1,5)]
e1=p[1];e2=(p[1]**2-p[2])/2
F=s.expand(p[3]-e1*p[2]+e2*p[1]);G=s.expand(p[4]-e1*p[3]+e2*p[2])
D=8*a*a-48*a*b+4*a+48*b*b+24*b-9
N=48*a*a*b-8*a*a-192*a*b*b-168*a*b+50*a+160*b**3+352*b*b+10*b-27
J=16*a*a-32*a*b-8*a+16*b*b-8*b+1
H=24*a*a-80*a*b-100*a+80*b*b+40*b+81
assert s.expand(F+3*(N-2*D*c)/16)==0
assert s.expand(s.resultant(N,2*D,a)+128*(4*b-25)*(4*b-9)**2*(4*b-1)**3)==0
assert s.cancel(G.subs(c,N/(2*D))-(4*b-9)*(4*b-1)*J*H/(64*D))==0
assert s.expand(J.subs({a:u*u,b:v*v})-16*(u*u-(v+s.Rational(1,2))**2)*(u*u-(v-s.Rational(1,2))**2))==0
M=s.Matrix([[48,-80,-100],[-80,160,40],[-100,40,162]])
assert [M.extract(I,I).det() for k in (1,2,3) for I in combinations(range(3),k)]==[48,160,162,1280,-2224,24320,-829440]
point={a:s.Rational(64,61)**2,b:s.Rational(1,122)**2}
assert H.subs(point)==0
assert D.subs(point)==s.Rational(70543200,13845841)
assert (N/(2*D)).subs(point)==s.Rational(26785,14884)
assert 26785%5==0 and 26785%25!=0 and 14884%5!=0
print('Sixth-mode identities, smoothness minors and nonexceptional nonsquare rejection verified.')
