#!/usr/bin/env python3
"""Independent finite-field 13-saturation certificate. See PROOF.md Section 5.16.
Requires the adjacent check_genus_two_richelot_certificate.py; no packages.
Polynomial coefficients are in ascending order. No rational-point completeness claim.
"""
import sys,math
from contextlib import redirect_stdout
from io import StringIO
with redirect_stdout(StringIO()):
 import check_genus_two_richelot_certificate as m
from fractions import Fraction as F
class Jac:
 def __init__(self,p,k):
  self.p=p; self.f=[1]
  for h in [[0,1],[-k,1],[5*k,1],[-275*k*k,-68*k,1]]: self.f=self.mul(self.f,h)
  self.O=([1],[0])
 def trim(self,a):
  a=[x%self.p for x in a]
  while len(a)>1 and not a[-1]:a.pop()
  return a
 def addp(self,a,b):return self.trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
 def sc(self,a,c):return self.trim([x*c for x in a])
 def mul(self,a,b):
  c=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):c[i+j]+=x*y
  return self.trim(c)
 def div(self,a,b):
  a=self.trim(a);b=self.trim(b);q=[0]*max(1,len(a)-len(b)+1)
  while a!=[0] and len(a)>=len(b):
   c=a[-1]*pow(b[-1],-1,self.p)%self.p;j=len(a)-len(b);q[j]=c
   a=self.addp(a,self.sc([0]*j+b,-c))
  return self.trim(q),a
 def quot(self,a,b):
  q,r=self.div(a,b);assert r==[0];return q
 def egcd(self,a,b):
  r0,r1=a,b;s0,s1=[1],[0];t0,t1=[0],[1]
  while r1!=[0]:
   q,r=self.div(r0,r1);r0,r1=r1,r;s0,s1=s1,self.addp(s0,self.sc(self.mul(q,s1),-1));t0,t1=t1,self.addp(t0,self.sc(self.mul(q,t1),-1))
  c=pow(r0[-1],-1,self.p)
  return self.sc(r0,c),self.sc(s0,c),self.sc(t0,c)
 def check(self,D):
  u,v=D;assert u[-1]==1 and len(u)<=3 and (v==[0] or len(v)<len(u))
  assert self.div(self.addp(self.f,self.sc(self.mul(v,v),-1)),u)[1]==[0]
 def add(self,A,B):
  u1,v1=A;u2,v2=B
  d1,e1,e2=self.egcd(u1,u2)
  d,c1,c2=self.egcd(d1,self.addp(v1,v2))
  s1=self.mul(c1,e1);s2=self.mul(c1,e2)
  u=self.quot(self.mul(u1,u2),self.mul(d,d))
  v=self.quot(self.addp(self.addp(self.mul(self.mul(s1,u1),v2),self.mul(self.mul(s2,u2),v1)),self.mul(c2,self.addp(self.mul(v1,v2),self.f))),d)
  v=self.div(v,u)[1]
  while len(u)>3:
   u=self.quot(self.addp(self.f,self.sc(self.mul(v,v),-1)),u);u=self.sc(u,pow(u[-1],-1,self.p))
   v=self.div(self.sc(v,-1),u)[1]
  u=self.sc(u,pow(u[-1],-1,self.p));out=(u,v);self.check(out);return out
 def times(self,n,D):
  R=self.O
  while n:
   if n&1:R=self.add(R,D)
   D=self.add(D,D);n//=2
  return R
 def frac(self,x):
  x=F(x);return x.numerator*pow(x.denominator,-1,self.p)%self.p
 def gens(self,k):
  z,Y,A,B,a,b=(F(-15,4),F(1425,32),-68,361,906,-6954) if k==-15 else (F(-11),F(1188),-20,361,342,-4446)
  P=([self.frac(-k*z),1],[self.frac(k*k*Y)])
  Q=(self.trim([B*k*k,A*k,1]),self.trim([b*k*k,a*k]))
  self.check(P);self.check(Q);return P,Q
base=[1]
for q in [[0,1],[-1,1],[5,1],[-275,-68,1]]:base=m.mul(base,q)
found={-15:[],-3:[]}
for p in [181,199]:
 if p in [11,19,53] or any(p%d==0 for d in range(2,math.isqrt(p)+1)):continue
 n1,n2,a1,a2,_=m.counts(base,p)
 assert (n1,n2,a1,a2)==({181:(192,32970,-10,154),199:(200,39854,0,126)}[p])
 for k in [-15,-3]:
  if len(found[k])>=2:continue
  n=p*p+1-m.legendre(k,p)*a1*(p+1)+a2
  if n%13 or n%169==0:continue
  J=Jac(p,k);P,Q=J.gens(k)
  assert J.times(n,P)==J.O and J.times(n,Q)==J.O
  P1,Q1=J.times(n//13,P),J.times(n//13,Q)
  if P1==J.O and Q1==J.O:continue
  if P1==J.O:row=(0,1)
  else:
   rat=next(a for a in range(13) if J.times(a,P1)==Q1);row=(1,rat)
  print('WITNESS',k,p,n,(n1,n2,a1,a2),'row',row,'P1',P1,'Q1',Q1,flush=True)
  if not found[k] or (found[k][0][0]*row[1]-found[k][0][1]*row[0])%13:found[k].append(row)
 if all(len(v)>=2 for v in found.values()):break
assert found == {-15:[(1,10),(0,1)],-3:[(1,7),(1,3)]}
print('Both rank-two subgroups are 13-saturated by independent finite-field tests.')
