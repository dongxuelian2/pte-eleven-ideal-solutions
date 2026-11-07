#!/usr/bin/env python3
"""Certified Coleman basis integrals and normalization; PROOF 5.24.
Requires existing SymPy. Shared etale arithmetic is imported by the mixed-height checker.
"""
import sympy as s
from fractions import Fraction as Q
x=s.symbols("x");p=13;N=3;G=4;MOD=p**(N+G);scale=p**G
def res(q):
 q=Q(q);assert q.denominator%p
 return q.numerator*pow(q.denominator,-1,MOD)%MOD
class Ring:
 def __init__(self,u,v):
  self.d=s.degree(u,x);self.u=[res(s.Poly(u,x).nth(i)) for i in range(self.d)]
  self.a=(0,1) if self.d==2 else((-self.u[0])%MOD,0)
  up=s.Poly(u,x,domain=s.QQ).clear_denoms()[1].set_modulus(13)
  assert s.gcd(up,up.diff()).degree()==0
  self.y=self.poly(s.Poly(v,x),self.a)
  self.inv(self.y)
 def add(self,a,b):return tuple((a[i]+b[i])%MOD for i in range(2))
 def neg(self,a):return tuple(-v%MOD for v in a)
 def mul(self,a,b):
  c=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[1]*b[1]]
  if self.d==2:return((c[0]-self.u[0]*c[2])%MOD,(c[1]-self.u[1]*c[2])%MOD)
  assert c[1]==c[2]==0
  return(c[0]%MOD,0)
 def pow(self,a,n):
  if n<0:return self.pow(self.inv(a),-n)
  r=(1,0)
  while n:
   if n&1:r=self.mul(r,a)
   a=self.mul(a,a);n//=2
  return r
 def inv(self,a):
  if self.d==1:return(pow(a[0],-1,MOD),0)
  u0,u1=self.u;norm=(a[0]**2-u1*a[0]*a[1]+u0*a[1]**2)%MOD
  b=pow(norm,-1,MOD)
  return((a[0]-u1*a[1])*b%MOD,-a[1]*b%MOD)
 def smul(self,q,a):return self.mul((res(q),0),a)
 def poly(self,P,a):
  r=(0,0)
  for c in P.all_coeffs():r=self.add(self.mul(r,a),(res(c),0))
  return r
 def tr(self,a):return (a[0] if self.d==1 else 2*a[0]-self.u[1]*a[1])%MOD
 def gt(self,B,den,yexp):
  # trace((scale*B/den)(a)*y^yexp)
  BP=B.set_domain(s.QQ).mul_ground(s.Rational(scale,den))
  return self.tr(self.mul(self.poly(BP,self.a),self.pow(self.y,yexp)))
 def tiny(self,f,vpoly):
  a13=self.pow(self.a,13);sigma=a13
  U=s.Poly(x**self.d+sum(s.Integer(b)*x**i for i,b in enumerate(self.u)),x)
  for j in range(4):
   sigma=self.add(sigma,self.neg(self.mul(self.poly(U,sigma),self.inv(self.poly(U.diff(),sigma)))))
  assert self.poly(U,sigma)==(0,0)
  y=self.poly(vpoly,sigma)
  assert all(c%13==0 for c in self.add(y,self.neg(self.pow(self.y,13))))
  h=self.add(a13,self.neg(sigma));assert all(c%13==0 for c in h)
  vals=[]
  for j in range(f.degree()-1):
   w=self.smul(s.Rational(1,2),self.mul(self.pow(sigma,j),self.inv(y)))
   dw=self.smul(s.Rational(j,2),self.mul(self.pow(sigma,j-1),self.inv(y))) if j else(0,0)
   correction=self.smul(s.Rational(1,4),self.mul(self.mul(self.pow(sigma,j),self.poly(f.diff(),sigma)),self.pow(y,-3)))
   dw=self.add(dw,self.neg(correction))
   val=self.add(self.mul(w,h),self.smul(s.Rational(1,2),self.mul(dw,self.mul(h,h))))
   vals.append(self.tr(val)%p**N)
  return vals

def main():
 for k,us,vs in [(-15,[x-s.Rational(225,4),x*x+1020*x+81225,x*x-276*x+11025],[s.Rational(320625,32),-13590*x-1564650,406350-7794*x]),(-3,[x-33,x*x+60*x+3249,x*x+65*x+528],[s.Integer(10692),-1026*x-40014,-1076*x-12804])]:
  f=s.Poly(x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k),x,domain=s.QQ);fp=f.diff();inv=s.invert(fp,f);R=s.Poly(f.as_expr().subs(x,x**p),x)-f**p
  rings=[Ring(u,v) for u,v in zip(us,vs)]
  for rg in rings:assert rg.poly(f,rg.a)==rg.mul(rg.y,rg.y)
  def reduce(N,m):
   gps=[0]*3
   while m:
    B=(N*inv).rem(f);A=(N-B*fp).exquo(f)
    for j,rg in enumerate(rings):gps[j]-=rg.gt(B,2*m-1,-(2*m-1))
    N=A+B.diff().mul_ground(s.Rational(2,2*m-1));m-=1
   while N.degree()>=4:
    e=N.degree();j=e-4;lead=N.LC();c=lead*s.Rational(2,2*e-3)
    for h,rg in enumerate(rings):gps[h]+=rg.gt(s.Poly(lead*x**j,x),2*e-3,1)
    P=s.Poly(x**j,x)*fp.mul_ground(s.Rational(1,2))
    if j:P+=s.Poly(j*x**(j-1),x)*f
    N-=P.mul_ground(c)
   return N,gps
  cols=[];gs=[[],[],[]]
  for i in range(4):
   C=s.Poly(0,x,domain=s.QQ);gps=[0]*3
   for n in range(4):
    a=(s.Poly(p*x**(p*(i+1)-1),x,domain=s.QQ)*R**n).mul_ground(s.binomial(-s.Rational(1,2),n))
    rr,gg=reduce(a,6+13*n);C+=rr;gps=[a+b for a,b in zip(gps,gg)]
   cols.append([res(C.nth(j))%p**N for j in range(4)])
   for h,g in enumerate(gps):
    g%=MOD;assert g%scale==0;gs[h].append(g//scale)
  M=s.Matrix.hstack(*[s.Matrix(c) for c in cols]);invI=(s.eye(4)-M.T).inv_mod(p**N)
  for label,rg,v,gg in zip("PQR",rings,vs,gs):
   tiny=rg.tiny(f,s.Poly(v,x));I=(invI*s.Matrix([2*(a-b) for a,b in zip(gg,tiny)])).applyfunc(lambda a:a%p**N)
   expected={-15:{'P':[1768,533,2062,1190],'Q':[195,117,182,1343],'R':[1963,650,14,22]},-3:{'P':[1079,533,852,402],'Q':[598,572,1772,1289],'R':[481,2158,1297,214]}}
   assert list(I)==expected[k][label]
   lam=[-s.Rational(1100*k**3,3)*I[0]+552*k*k*I[1]+128*k*I[2]-3*I[3],-s.Rational(16*k,3)*I[1]-I[2]]
   want={-15:{"P":[1545,1032],"Q":[1301,390],"R":[1238,1455]},-3:{"P":[2021,1085],"Q":[1409,789],"R":[614,276]}}
   assert [int(q)%2197 for q in lam]==want[k][label]
   print(k,label,'I',list(I),'normalization', [int(q)%2197 for q in lam],flush=True)
 

def reciprocity():
 t=s.symbols("t")
 for k,V,wanted in [(-15,s.Rational(33,2)*x*x+3240*x-s.Rational(448875,2),[0,0,33,314]),(-3,-10*x*x+426*x+7524,[0,0,2177,1096])]:
  f=x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k)
  H=s.series(s.expand(t**10*f.subs(x,t**-2))**s.Rational(-1,2),t,0,6).removeO()
  q=s.series(V.subs(x,t**-2)*t**5*H,t,0,4).removeO()
  lg=s.series(-2*q-2*q**3/3,t,0,4).removeO()
  residues=[res(s.expand(lg*(-t**(2-2*j)*H)).coeff(t,-1))%2197 for j in range(4)]
  assert residues==wanted
 print("Independent principal-divisor reciprocity PASS")

if __name__=='__main__':
 main()
 reciprocity()
