import sympy as s
a,b,c,u,v=s.symbols('a b c u v')
def moments(z,q,k):
 p=s.Poly(sum((z+d)**(2*k) for d in (-1,0,1)),z)
 assert all(mo[0]%2==0 for mo,co in p.terms())
 return sum(co*q**(mo[0]//2) for mo,co in p.terms())
p=[None]+[s.Rational(1,4)**k+moments(u,a,k)+c**k-moments(v,b,k) for k in range(1,5)]
e1=p[1];e2=(p[1]**2-p[2])/2
F=s.factor(p[3]-e1*p[2]+e2*p[1]);G=s.factor(p[4]-e1*p[3]+e2*p[2])
f=s.cancel(F/(a-b));g=s.cancel(G/(a-b));D=4*a-8*b-7
H=8*a*a-40*a*b-58*a+80*b*b+50*b+95
expected=3*(4*a*a-32*a*b+12*a*c-29*a+40*b*b-24*b*c+106*b-21*c+34)/4
assert s.expand(f-expected)==0
assert s.expand(f.subs(a,2*b+s.Rational(7,4))+3*(4*b-9)*(4*b-1)/8)==0
C=-(4*a*a-32*a*b-29*a+40*b*b+106*b+34)/(3*D)
assert s.cancel(g.subs(c,C)-2*((a-b-1)**2-4*b)*H/(3*D))==0
M=s.Matrix([[16,-40,-58],[-40,160,50],[-58,50,190]])
from itertools import combinations
minors=[M.extract(I,I).det() for k in (1,2,3) for I in combinations(range(3),k)]
assert minors==[16,160,190,960,-324,27900,-163840]
assert s.expand((8*a-20*b-29)**2-(-240*b*b+760*b+81)-8*H)==0
for uu,vv in [(s.Rational(3,2),s.Rational(1,2)),(s.Rational(5,2),s.Rational(3,2))]:
 assert H.subs({a:uu*uu,b:vv*vv})==D.subs({a:uu*uu,b:vv*vv})==0
q=-240*v**4+760*v**2+81
assert s.degree(s.gcd(q,s.diff(q,v)),v)==0
assert q.subs(v,0)==9**2
print('Fifth mode elimination, exceptional branch, smoothness minors and quotient identity verified.')
