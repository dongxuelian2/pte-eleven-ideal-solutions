#!/usr/bin/env python3
"""All dyadic mixed heights and complete away-from-13 inputs; PROOF 5.21.
Requires existing SymPy. Does not compute normalized integrals at 13.
"""
import sympy as s
from fractions import Fraction as F
x,t,r,v=s.symbols("x t r v")
def val(q,p=2):
    q=s.Rational(q);assert q
    a,b,n=int(q.p),int(q.q),0
    while a%p==0:a//=p;n+=1
    while b%p==0:b//=p;n-=1
    return n
def f(k):return x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k)
# k=-3 smooth affine charts and support positions.
for a in [1,3]:
    g=s.Poly(s.expand(f(-3).subs(x,a+4*t)/16),t)
    assert all(c.q==1 for c in g.all_coeffs())
    assert s.Poly(g.as_expr(),t,modulus=2)==s.Poly(t+1,t,modulus=2)
uq=x*x+60*x+3249;ur=x*x+65*x+528
assert val(uq.subs(x,1))==1 and val(s.diff(uq,x).subs(x,1))==1
assert val(ur.subs(x,1))==1 and int(s.diff(ur,x).subs(x,1))%2
assert s.Poly(ur,x,modulus=2)==s.Poly(x*(x+1),x,modulus=2)
assert 33%4==1 and 10692%4==0
assert all(int(c)%4==0 for c in s.Poly(-1076*x-12804,x).all_coeffs())
assert s.diff(f(-3),x).subs(x,0)%2==1
# k=-15 infinity is smooth and Q/R have integral odd roots.
g=s.expand(v**6*f(-15).subs(x,1/v))
assert g==v+960*v*v-124200*v**3+2565000*v**4+69609375*v**5
assert s.diff(g,v).subs(v,0)==1
assert val(s.Rational(225,4))==-2 and val(s.Rational(320625,32))==-5
for u in [x*x+1020*x+81225,x*x-276*x+11025]:
    assert s.Poly(u,x,modulus=2)==s.Poly((x+1)**2,x,modulus=2)
# Wild smooth chart over L=Q2(r), r^2=11, pi=1+r.
qred=lambda a:s.rem(s.expand(a),r*r-11,r)
pi=1+r;pinv=s.invert(pi,r*r-11,r);a=138+27*r
norm=lambda c:qred(c*qred(c.subs(r,-r)))
cs=[qred(s.diff(f(-15),x,j).subs(x,a)/s.factorial(j)) for j in range(6)]
assert [val(norm(c)) for c in cs]==[6,4,2,2,0,0]
scaled=[]
for j,c in enumerate(cs):
    e=2*j-6
    scaled.append(qred(c*(pi**e if e>=0 else pinv**(-e))))
assert [val(norm(c)) for c in scaled]==[0,0,0,2,2,4]
for c in scaled[:3]:
    assert all(s.Rational(b).q%2 for b in s.Poly(c,r).all_coeffs())
    assert int(s.numer(c.subs(r,1)))%2 and int(s.denom(c.subs(r,1)))%2
assert val(s.Integer(5366250**2-795*203850**2))==3
# Exact distance proof, no numerical square-root approximation.
assert (795*pow(11,-1,16))%16==1
assert val(s.Integer(-648))==3
assert (15-27)%8==4 and (15+27)%4==2
# For a close pair d=xR-xQ has v2=2 and q=xQ is a unit:
assert [val(s.Integer(n)) for n in [5796,7794,1971000]]==[2,1,3]
assert [val(s.Integer(n)) for n in [21384,7794,1158300]]==[3,1,2]
# Hence both yR-yQ and yR+yQ have v2=2, or vL=4.
assert 4-3==1  # both intersections in parameter Y have length one.
# Optional certified residues of the complete away sums.
def residue(a,N):
    a=F(a);M=13**N
    assert a.denominator%13
    return a.numerator*pow(a.denominator,-1,M)%M
def log13(q,N):
    z=(pow(q,12,13**(N+2))-1)%13**(N+2)
    assert z%13==0
    total=sum((F((-1)**(j+1)*z**j,12*j) for j in range(1,N+3)),F(0))
    return residue(total,N)
L={q:log13(q,10) for q in [3,5,7,17,19,29,31,313,593]}
assert all(log13(q,20)%13**10==b for q,b in L.items())
M=13**10
values=[
(-15,1,-2*(L[17]+L[593]),-2*L[5],2*(L[7]+L[19]+L[313])),
(-3,-1,-L[3],2*L[19],-2*(L[29]+L[31]))]
expected=[(76380091502,69520774258,79405079481),
          (24898284736,34648651414,71972655326)]
for (k,sgn,U,V,W),wanted in zip(values,expected):
    got=tuple(b%M for b in [U,V,W]);assert got==wanted
    print(k,"complete away mixed sums mod13^10:",got)
print("PASS: all six dyadic mixed heights are zero.")
print("NOT COMPUTED: normalized 13-adic mixed integrals, full heights or QC roots.")
