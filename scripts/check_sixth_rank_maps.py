"""Exact maps and rational points; the isogeny/rank argument is in the note."""
import sympy as s
u,v,x=s.symbols('u v x');a=u*u;b=v*v
H=24*a*a-80*a*b-100*a+80*b*b+40*b+81
fv=-80*x**4+760*x*x+139
fu=-5*(x*x-1)*(x*x-19)
fr=139*x**4+1120*x*x-1520
Y=12*a-20*b-25;Z=10*b-5*a+s.Rational(5,2)
R=u/v;Q=(81/v**2-50*R**2+20)/2
assert s.expand(Y*Y-fv.subs(x,v)-6*H)==0
assert s.expand(Z*Z-fu.subs(x,u)-s.Rational(5,4)*H)==0
assert s.cancel(Q*Q-fr.subs(x,R)-s.Rational(81,4)*H/v**4)==0
for f in (fv,fu,fr):assert s.degree(s.gcd(f,s.diff(f,x)),x)==0
Cv=[(s.Rational(1,2),s.Integer(18)),(s.Rational(3,2),s.Integer(38)),(s.Rational(5,2),s.Integer(42)),(s.Rational(41,14),s.Rational(1362,49)),(s.Rational(13,86),s.Rational(23118,1849))]
Cpoints=[(s.Integer(1),s.Rational(1,2)),(s.Integer(2),s.Rational(1,2)),(s.Integer(3),s.Rational(3,2)),(s.Integer(3),s.Rational(5,2)),(s.Integer(4),s.Rational(5,2)),(s.Rational(64,61),s.Rational(1,122))]
Cu=[];Cr=[]
for uu,vv in Cpoints:
 assert H.subs({u:uu,v:vv})==0
 Cu.append((uu,Z.subs({u:uu,v:vv})))
 Cr.append((uu/vv,Q.subs({u:uu,v:vv})))
for name,f,points in [('Ev',fv,Cv),('Eu',fu,Cu),('Er',fr,Cr)]:
 signed={(sg*t,sh*y) for t,y in points for sg in (-1,1) for sh in (-1,1)}
 assert all(y*y==f.subs(x,t) for t,y in signed)
 assert len(signed)=={'Ev':20,'Eu':18,'Er':20}[name]
 print(name,'distinct rational affine points:',len(signed))
print('All three smooth genus-one quotients exceed the rational torsion maximum of 16.')

A,B,C,t=s.symbols('A B C t')
X=A*t*t
assert s.expand(A*A*t*t*(A*t**4+B*t*t+C)-X*(X*X+B*X+A*C))==0
X=x+A+B/x
assert s.cancel(x*(x*x+A*x+B)*(1-B/x**2)**2-X*(X*X-2*A*X+A*A-4*B))==0
print('Even-quartic cover and standard 2-isogeny identities verified.')
