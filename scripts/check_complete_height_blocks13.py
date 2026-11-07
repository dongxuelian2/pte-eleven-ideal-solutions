#!/usr/bin/env python3
"""Completed global height blocks and logarithmic quadratic forms; PROOF 5.26--5.27.
Requires adjacent check_coleman_basis13.py and existing SymPy.
Does not compute local QC functions or roots.
"""
import sympy as s
from check_coleman_basis13 import Ring,res,p,N,MOD,scale,x

for k,pa,pb,original in [(-3,s.Integer(33),s.Integer(10692),[(x*x-s.Rational(14122312959396,446872195225)*x+s.Rational(119757412617129,446872195225),(s.Integer(94388884819291490706)*x-s.Integer(1395182795689585717794))/s.Integer(298727359424984125))])]:
 odd=x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k)
 f=s.Poly(s.cancel(x**6*odd.subs(x,pa+1/x)/(pb*pb)),x,domain=s.QQ)
 assert f.degree()==6 and f.LC()==1
 ff=f.clear_denoms()[1].set_modulus(13);assert s.gcd(ff,ff.diff()).degree()==0
 pairs=[]
 for U,V in original:
  UT=s.Poly(s.cancel(x*x*U.subs(x,pa+1/x)),x,domain=s.QQ).monic()
  VT=s.Poly(s.cancel(V.subs(x,pa+1/x)*x**3/pb),x,domain=s.QQ).rem(UT)
  assert (f-VT*VT).rem(UT).is_zero
  pairs.append((UT.as_expr(),VT.as_expr()))
 rings=[Ring(u,v) for u,v in pairs]
 fp=f.diff();inv=s.invert(fp,f);R=s.Poly(f.as_expr().subs(x,x**p),x)-f**p
 def reduce(N,m):
  gps=[0]*len(rings)
  while m:
   B=(N*inv).rem(f);A=(N-B*fp).exquo(f)
   for j,rg in enumerate(rings):gps[j]-=rg.gt(B,2*m-1,-(2*m-1))
   N=A+B.diff().mul_ground(s.Rational(2,2*m-1));m-=1
  while N.degree()>=5:
   e=N.degree();j=e-5;lead=N.LC();c=lead*s.Rational(2,2*e-4)
   for h,rg in enumerate(rings):gps[h]+=rg.gt(s.Poly(lead*x**j,x),2*e-4,1)
   P=s.Poly(x**j,x)*fp.mul_ground(s.Rational(1,2))
   if j:P+=s.Poly(j*x**(j-1),x)*f
   N-=P.mul_ground(c)
  return N,gps
 cols=[];gs=[[] for _ in rings]
 for i in range(5):
  C=s.Poly(0,x,domain=s.QQ);gps=[0]*len(rings)
  for n in range(4):
   aa=(s.Poly(p*x**(p*(i+1)-1),x,domain=s.QQ)*R**n).mul_ground(s.binomial(-s.Rational(1,2),n))
   rr,gg=reduce(aa,6+13*n);C+=rr;gps=[a+b for a,b in zip(gps,gg)]
  cols.append([res(C.nth(j))%p**N for j in range(5)])
  for h,g in enumerate(gps):
   g%=MOD;assert g%scale==0;gs[h].append(g//scale)
 M=s.Matrix.hstack(*[s.Matrix(c) for c in cols]);invI=(s.eye(5)-M.T).inv_mod(p**N)
 T=s.symbols("T")
 assert all(int(c)%13**3==0 for c in s.Poly(M.charpoly(T).as_expr()-(T-13)*(T**4+6*T*T+169),T).all_coeffs())
 for label,rg,(_,vv),gg in zip(("Q","R") if k==-15 else ("R",),rings,pairs,gs):
  tiny=rg.tiny(f,s.Poly(vv,x));I=(invI*s.Matrix([2*(a-b) for a,b in zip(gg,tiny)])).applyfunc(lambda a:a%p**N)
  assert list(I)==[273,1859,143,1464,1266]
  assert int(-2*I[2])%2197==1911
  assert (1911-2021*1521-1085*1560)%2197==663
  print('T EVEN INTEGRALS',list(I),'raw',int(-2*I[2])%2197,flush=True)

# Exact Cantor identities and all finite-prime guards for the moved divisor.
def valuation(q,l):
 q=s.Rational(q);a,b=int(q.p),int(q.q)
 if a==0:return 10**6
 n=0
 while a%l==0:a//=l;n+=1
 while b%l==0:b//=l;n-=1
 return n
f=s.Poly(x*(x+3)*(x-15)*(x*x+204*x-2475),x,domain=s.QQ)
u=s.Poly(x*x+65*x+528,x,domain=s.QQ)
v=s.Poly(-1076*x-12804,x,domain=s.QQ)
for j in range(3):
 V=v+u.mul_ground((10692-v.eval(33))/u.eval(33))
 un=(f-V*V).exquo(s.Poly(x-33,x)*u).monic()
 v=(-V).rem(un);u=un
 assert (f-v*v).rem(u).is_zero
assert u.as_expr()==x*x-s.Rational(14122312959396,446872195225)*x+s.Rational(119757412617129,446872195225)
assert v.as_expr()==(s.Integer(94388884819291490706)*x-s.Integer(1395182795689585717794))/s.Integer(298727359424984125)
assert valuation(s.discriminant(u.as_expr(),x),13)==0
assert valuation(s.resultant(u.as_expr(),v.as_expr(),x),13)==0
assert valuation(u.eval(33),13)==0
assert all(valuation(c,13)>=0 for c in u.all_coeffs()+v.all_coeffs())
rr=u.eval(33)
assert s.factorint(abs(int(s.numer(rr))))=={2:1,3:5,7:1,61:1,569:1,1188727:1}
assert s.factorint(int(s.denom(rr)))=={5:2,133697:2}
assert s.factorint(298727359424984125)=={5:3,133697:3}
for l,sign in [(7,1),(61,1),(569,-1),(1188727,-1)]:
 assert s.isprime(l) and valuation(rr,l)==1 and 10692%l
 assert valuation(v.eval(33)-sign*10692,l)>=1
 assert valuation(v.eval(33)+sign*10692,l)==0
for l in [11,19,53]:
 assert valuation(rr,l)==0
d,up=u.clear_denoms();up=up.primitive()[1]
for l,want in [(5,4),(133697,41333)]:
 assert valuation(up.nth(2),l)==2
 assert valuation(up.nth(1),l)==valuation(up.nth(0),l)==0
 root=(-int(up.nth(0))*pow(int(up.nth(1)),-1,l))%l
 assert root==want and root!=33%l
shift=s.Poly(up.as_expr().subs(x,x+1),x)
assert [valuation(shift.nth(i),2) for i in range(3)]==[1,1,0]
t=s.symbols('t')
G=s.Poly(up.as_expr().subs(x,-3-9*t)/81,t)
assert G.as_expr()==446872195225*t*t+1867060681194*t+2051187669782
assert valuation(G.eval(-4),3)==1 and valuation(G.diff().eval(-4),3)==0
assert s.Poly(G,t,modulus=3).as_expr()==t*t-1
VT=s.cancel(v.as_expr().subs(x,-3-9*t)/243)
assert s.cancel(VT+38*(91996963761492681*t+181757575281293068)/298727359424984125)==0
# Reduce rational ordinates modulo3 directly, not through the13-adic res helper.
def red(q,l):
 q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,l)%l
assert red(VT.subs(t,2),3)==1 and 44%3==2
def lg(q):
 zz=pow(q,12,13**6)-1
 return res(sum((s.Rational((-1)**(j+1)*zz**j,12*j) for j in range(1,7)),s.Rational(0)))%2197
away=(lg(3)-2*lg(7)-2*lg(61)+2*lg(569)+2*lg(1188727))%2197
assert away==1560
assert (663+away)%2197==26
assert (-2*lg(29)-2*lg(31))%2197==871
assert (2*lg(7)+2*lg(19)+2*lg(313))%2197==390
print("Exact divisor move and complete finite correction PASS")
# Global reconstruction and coordinate/precision checks.
for k,IP,IQ,H,Aexpected in [
 (-15,[1768,533],[195,117],[[923,1586],[1586,689]],[[37,131],[131,135]]),
 (-3,[1079,533],[598,572],[[2093,676],[676,1027]],[[90,114],[114,36]])]:
 HH=s.Matrix(H)
 if k==-15:
  assert ((1248-1950)*pow(4,-1,2197))%2197==H[0][0]
  assert 1950*pow(4,-1,2197)%2197==H[0][1]
  assert (312-1950)*pow(4,-1,2197)%2197==H[1][1]
 else:
  assert (26-1274)*pow(12,-1,2197)%2197==H[0][0]
  assert (26-4*1274)*pow(12,-1,2197)%2197==H[0][1]
  assert (H[0][1]-793*pow(4,-1,2197))%2197==H[1][1]
 L=s.Matrix([[IP[0]//13,IQ[0]//13],[IP[1]//13,IQ[1]//13]])*pow(2,-1,169)
 L=L.applyfunc(lambda z:z%169);assert int(L.det())%13==6
 H0=HH/13;li=L.inv_mod(169);A=(li.T*H0*li).applyfunc(lambda z:z%169)
 assert A==s.Matrix(Aexpected)
 assert (13*L.T*A*L-HH).applyfunc(lambda z:z%2197)==s.zeros(2,2)
 assert int(H0.det())%13==({-15:7,-3:5}[k])
 print(k,"global height matrix",HH.tolist(),"quadratic coefficients",A.tolist())
print("Both global blocks certified modulo13^3; local QC functions and roots NOT computed.")
