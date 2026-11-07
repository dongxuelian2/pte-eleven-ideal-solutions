"""Exact identities for four symmetric size-ten shift-chain modes.
Requires SymPy. This does not enumerate rational points or find a PTE witness.
"""
from itertools import combinations_with_replacement
import sympy as s
b,a,x,z=s.symbols('b a x z')
def offsets(L):return [s.Rational(2*j-L+1,2) for j in range(L)]
def hm(L,k,w):
    return s.expand(sum(s.binomial(2*k,2*j)*w**(k-j)*sum(d**(2*j) for d in offsets(L)) for j in range(k+1)))
def direct(L,k,w):
    p=s.Poly(sum((z+d)**(2*k) for d in offsets(L)),z)
    assert all(mon[0]%2==0 for mon,c in p.terms())
    return s.expand(sum(c*w**(mon[0]//2) for mon,c in p.terms()))
patterns=[]
for m in range(1,6):
 for na in range(6):
  for aa in combinations_with_replacement(range(1,6),na):
   for bb in combinations_with_replacement(range(1,6),5-na):
    if sum(aa)==5-m and sum(bb)==5:patterns.append((m,aa,bb))
assert len(patterns)==18
F=[-4*(a-b)*(a-5*b-4)*(a*a-2*a*b-2*a+b*b-2*b+1),
(b-1)*(4*a-10*b-9)*(16*a*a-32*a*b-8*a+16*b*b-8*b+1)/8,
-4*(a-b)*(b-4)*(b-1)*(3*a-5*b-3),
4*(b-9)*(b-4)*(b-1)*(4*a-5*b-1)]
h=[800*b**6-2625*b**5+12050*b**4-529*b**3-54448*b**2-53744*b-14400,
675*b**5-7950*b**4-26125*b**3-35704*b**2-23248*b-6144,
900*b**5-7100*b**4-21675*b**3+71658*b**2+3321*b-93312,
50*b**6-1875*b**5+20600*b**4-50624*b**3+102912*b**2-2392064*b-3686400]
relations=[5*b+4,(10*b+9)/4,(5*b+3)/3,(5*b+1)/4]
expected_D=[-4*h[0],-b*h[1]/4,-4*b*h[2]/243,-h[3]/64]
curves=[-h[0],-h[1],-3*h[2],-h[3]]
for m in range(1,5):
 for L in (2,5-m):
  for k in range(1,5):assert s.expand(hm(L,k,a)-direct(L,k,a))==0
 p=[None]+[s.expand(sum(d**(2*k) for d in offsets(2*m))/2+hm(5-m,k,a)-hm(2,k,b)) for k in range(1,5)]
 e1=p[1];e2=(p[1]**2-p[2])/2;e3=(p[3]-e1*p[2]+e2*p[1])/3
 assert s.expand(p[4]-e1*p[3]+e2*p[2]-e3*p[1]-F[m-1])==0
 K=s.expand((x**3-e1*x*x+e2*x-e3).subs(a,relations[m-1]))
 assert s.expand(s.discriminant(K,x)-expected_D[m-1])==0
 f=curves[m-1]
 assert s.degree(s.gcd(f,s.diff(f,b)),b)==0
 assert s.degree(f,b)==(6,5,5,6)[m-1]
 roots=s.polys.polytools.ground_roots(K.subs(b,m*m),x)
 assert len(roots)==3 and all(s.sqrt(r).is_Rational for r in roots)
 u=s.Rational(m+5,2);v=s.Integer(m)
 A=set(offsets(2*m))|{sg*u+d for sg in (-1,1) for d in offsets(5-m)}
 B={sg*s.sqrt(r) for sg in (-1,1) for r in roots}|{sg*v+d for sg in (-1,1) for d in offsets(2)}
 assert A==B==set(offsets(10))
 print('mode',m,': identities, squarefree genus-two model, and trivial control verified')
print('18 chain types; no rational-point classification or size-eleven witness claimed.')
