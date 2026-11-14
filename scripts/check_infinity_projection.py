#!/usr/bin/env python3
"""Certified infinity projection for sixth-mode curves; see INFINITY_PROJECTION.md.
Requires existing SymPy. No installation or network actions.
"""
import sympy as s
p=13; N=25; S=8; mod=p**N; x=s.symbols('x')
def trim(A):
 while len(A)>1 and not A[-1]:A.pop()
 return A
def add(A,B):return trim([((A[i] if i<len(A) else 0)+(B[i] if i<len(B) else 0))%mod for i in range(max(len(A),len(B)))])
def scale(A,c):return trim([a*c%mod for a in A])
def sub(A,B):return add(A,scale(B,-1))
def mul(A,B):
 C=[0]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]=(C[i+j]+a*b)%mod
 return trim(C)
def divrem(A,B):
 A=A[:];Q=[0]*max(1,len(A)-len(B)+1);ib=pow(B[-1],-1,mod)
 while len(A)>=len(B) and A!=[0]:
  j=len(A)-len(B);c=A[-1]*ib%mod;Q[j]=c
  for i,b in enumerate(B):A[i+j]=(A[i+j]-c*b)%mod
  trim(A)
 return trim(Q),A
def rem(A,B):return divrem(A,B)[1]
def deriv(A):return trim([i*A[i]%mod for i in range(1,len(A))] or [0])
def pw(A,n):
 B=[1]
 while n:
  if n%2:B=mul(B,A)
  A=mul(A,A);n//=2
 return B
def pd(A,d):
 t=0
 while d%p==0:t+=1;d//=p
 q=p**t
 assert all(a%q==0 for a in A),(d,t)
 return trim([a//q*pow(d,-1,mod)%mod for a in A])
def rat(q):
 q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,mod)%mod
for delta,Npoly,Dpoly,c in [(1,80-x*x,4*x-30,1),(5,25*x*x-120*x+80,10*x*x-8,5625)]:
 fq=s.Poly(s.expand((Npoly*Npoly+5*Dpoly*Dpoly)*(Npoly*Npoly+275*Dpoly*Dpoly))/s.Integer(c*c),x)
 f=[rat(fq.nth(i)) for i in range(9)];fp=deriv(f)
 fi=s.invert(s.Poly(sum(z*x**i for i,z in enumerate(fp)),x,modulus=p),s.Poly(sum(z*x**i for i,z in enumerate(f)),x,modulus=p))
 inv=[int(fi.nth(i))%mod for i in range(8)]
 for j in range(6):inv=rem(mul(inv,sub([2],rem(mul(fp,inv),f))),f)
 assert rem(mul(fp,inv),f)==[1]
 fx=[0]*(8*p+1)
 for i,z in enumerate(f):fx[i*p]=z
 R=sub(fx,pw(f,p));assert all(z%p==0 for z in R)
 def reduce(A,m):
  while m:
   B=rem(mul(A,inv),f);C,rr=divrem(sub(A,mul(B,fp)),f);assert rr==[0]
   Bdiv=pd(B,2*m-1)
   A=add(C,scale(deriv(Bdiv),2));m-=1
  while len(A)>7:
   e=len(A)-1;j=e-7;c0=pd([A[-1]],2*e-6)[0]
   B=[0]*j+[c0]
   A=sub(A,add(mul(B,fp),scale(mul(deriv(B),f),2)))
  return A
 mats=[]
 for i in range(7):
  A=[0]
  for n in range(5):
   term=scale([0]*(p*(i+1)-1)+pw(R,n),p**(S+1)*rat(s.binomial(-s.Rational(1,2),n)))
   A=add(A,reduce(term,6+13*n))
  assert all(z%p**S==0 for z in A)
  mats.append([(A[j] if j<len(A) else 0)//p**S%p**4 for j in range(7)])
 MM=s.Matrix.hstack(*map(s.Matrix,mats));print('DELTA',delta,'M',MM.tolist(),flush=True)
 u=s.symbols('u')
 rev=sum(fq.nth(8-j)*u**j for j in range(9))
 sr=s.series(rev**(-s.Rational(1,2)),u,0,4).removeO()
 rr=[sr.coeff(u,j) for j in range(4)]
 BB=s.zeros(7)
 for i in range(3):BB[i,i]=1
 for i in range(3):BB[i+4,i+3]=1;BB[3,i+3]=-rr[i+1]
 BB[3,6]=1
 def rsmall(q,n=4):
  q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,p**n)%p**n
 TT=(BB.inv()*MM*BB).applyfunc(lambda z:rsmall(z))
 assert list(TT[6,:])==[0]*6+[13]
 FF=TT[:6,:6];vv=TT[:6,6];AA=FF-13*s.eye(6)
 ai=AA.inv();sol=ai*vv
 def vp(q):
  q=s.Rational(q);nn,dd=int(q.p),int(q.q);z=0
  if nn==0:return 999
  while nn%13==0:nn//=13;z+=1
  while dd%13==0:dd//=13;z-=1
  return z
 print('rho equation','v',list(vv),'inverse minval',min(vp(z) for z in ai),'sol minval',min(vp(z) for z in sol),flush=True)
 avec=[rsmall(z,3) for z in sol]
 print('PSI_RHO_eta',delta,avec,flush=True)
 t=s.symbols('T')
 expected=(t-13)*(t*t-2*t+13)*(t**4+6*t*t+169)
 assert all(int(z)%p**4==0 for z in s.Poly(MM.charpoly(t).as_expr()-expected,t).all_coeffs())
 assert (AA*s.Matrix(avec)-vv).applyfunc(lambda z:z%p**3)==s.zeros(6,1)

import sympy as s
x,u=s.symbols('x u')
for delta,N,D,M,c,avec in [
 (1,80-x*x,4*x-30,2*x*x-30*x+160,1,[1433,392,1539,1586,728,1248]),
 (5,25*x*x-120*x+80,10*x*x-8,-60*x*x+100*x-48,5625,[1194,1844,368,1690,533,1898])]:
 f=s.Poly(s.expand((N*N+5*D*D)*(N*N+275*D*D))/s.Integer(c*c),x,domain=s.QQ);dp=s.Poly(D,x,domain=s.QQ);fp=f.diff()
 K=s.cancel((s.diff(N,x)*D-N*s.diff(D,x))/M)
 forms=[(K*M*N*N/(D*D),2),(K*(N**3+s.Rational(15,4)*N*N*D)/D**2,2),(K*(N**4+s.Rational(4125,32)*N*N*D*D+525*N*D**3)/D**3,3)]
 def reduce(R,m):
  P=s.Poly(s.cancel(R*D**m/c),x,domain=s.QQ)
  while m>1:
   A=(-P*s.invert(dp.diff()*f*(m-1),dp)).rem(dp)
   exact=A.diff()*dp*f-A*dp.diff()*f*(m-1)+A*dp*fp/2
   P=(P-exact).exquo(dp);m-=1
  P=P.exquo(dp)
  while P.degree()>=7:
   e=P.degree();j=e-7;lead=P.LC();A=s.Poly(x**j,x)
   exact=A.diff()*f+A*fp/2
   P-=exact*(lead/s.Integer(e-3))
  return s.Matrix([2*P.nth(i) for i in range(7)])
 rev=sum(f.nth(8-j)*u**j for j in range(9));rs=s.series(rev**(-s.Rational(1,2)),u,0,4).removeO()
 B=s.zeros(7)
 for i in range(3):B[i,i]=1
 for i in range(3):B[i+4,i+3]=1;B[3,i+3]=-rs.coeff(u,i+1)
 B[3,6]=1
 W=[B.inv()*reduce(R,m) for R,m in forms]
 assert all(v[6]==0 for v in W)
 S=s.Matrix.hstack(s.eye(6)[:,:3],*[v[:6,0] for v in W])
 sol=S.inv()*s.Matrix(avec)
 def r(q):
  q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,2197)%2197
 print('delta',delta,'K',K,'r',[rs.coeff(u,i) for i in range(4)],'splitdet',s.factor(S.det()),flush=True)
 print('W_eta',s.Matrix.hstack(*W).tolist(),flush=True)
 print('lambda_rho_omega012',[r(z) for z in sol[:3]],'lambda_nuinfplus',[r(-z) for z in sol[:3]],flush=True)
