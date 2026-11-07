#!/usr/bin/env python3
"""Certified remaining Q/R integrals; PROOF 5.26.
Existing SymPy only. Does not compute the local QC functions or roots.
"""
import sympy as s
p=13; N=20; S=6; mod=p**N; x=s.symbols('x')
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
def run(a):
 a=int(a)
 for j in range(N):a=(a-(a*a+60*a+3249)*pow(2*a+60,-1,mod))%mod
 assert (a*a+60*a+3249)%mod==0
 b=(-1026*a-40014)%mod
 f0=[0]+[int(s.Poly(s.expand((x+3)*(x-15)*(x*x+204*x-2475)),x).nth(i)) for i in range(5)]
 f=[0]*7
 for j,c in enumerate(f0):
  for l in range(j+1):f[6-l]=(f[6-l]+c*int(s.binomial(j,l))*pow(a,j-l,mod)*pow(b,-2,mod))%mod
 assert f[6]==1
 fp=deriv(f)
 fi=s.invert(s.Poly(sum(c*x**i for i,c in enumerate(fp)),x,modulus=p),s.Poly(sum(c*x**i for i,c in enumerate(f)),x,modulus=p))
 inv=[int(fi.nth(i))%mod for i in range(6)]
 for j in range(6):inv=rem(mul(inv,sub([2],rem(mul(fp,inv),f))),f)
 assert rem(mul(fp,inv),f)==[1]
 fx=[0]*(6*p+1)
 for i,c in enumerate(f):fx[i*p]=c
 R=sub(fx,pw(f,p));assert all(c%p==0 for c in R)
 def reduce(A,m):
  ps=[]
  while m:
   B=rem(mul(A,inv),f);C,rr=divrem(sub(A,mul(B,fp)),f);assert rr==[0]
   Bdiv=pd(B,2*m-1);ps.append((scale(Bdiv,-1),1-2*m))
   A=add(C,scale(deriv(Bdiv),2));m-=1
  while len(A)>5:
   e=len(A)-1;j=e-5;c=pd([A[-1]*pow(f[-1],-1,mod)%mod],2*e-4)[0]
   B=[0]*j+[c];G=add(mul(B,fp),scale(mul(deriv(B),f),2))
   A=sub(A,G);ps.append((B,1))
  return A,ps
 mats=[];prims=[]
 for i in range(5):
  A=[0];ps=[]
  for n in range(4):
   term=scale([0]*(p*(i+1)-1)+pw(R,n),p**(S+1)*rat(s.binomial(-s.Rational(1,2),n)))
   out,pp=reduce(term,6+13*n);A=add(A,out);ps+=pp
  assert all(c%p**S==0 for c in A)
  mats.append([(A[j] if j<len(A) else 0)//p**S%p**3 for j in range(5)]);prims.append(ps)
 print('a',a%p**3,'M columns',mats,flush=True)
 U1=65;U0=528;zero=(0,0);one=(1,0)
 ad=lambda A,B:((A[0]+B[0])%mod,(A[1]+B[1])%mod)
 sc=lambda A,c:(A[0]*c%mod,A[1]*c%mod)
 su=lambda A,B:ad(A,sc(B,-1))
 def mu(A,B):return ((A[0]*B[0]-U0*A[1]*B[1])%mod,(A[0]*B[1]+A[1]*B[0]-U1*A[1]*B[1])%mod)
 co=lambda A:((A[0]-U1*A[1])%mod,-A[1]%mod)
 def iv(A):
  C=co(A);d=mu(A,C);assert d[1]==0 and d[0]%p
  return sc(C,pow(d[0],-1,mod))
 def po(A,n):
  if n<0:return po(iv(A),-n)
  B=one
  while n:
   if n%2:B=mu(B,A)
   A=mu(A,A);n//=2
  return B
 def ev(A,z):
  B=zero
  for c in A[::-1]:B=ad(mu(B,z),(c,0))
  return B
 U=iv((-a%mod,1));V=sc(mu((-12804%mod,-1076%mod),po(U,3)),pow(b,-1,mod))
 assert mu(V,V)==ev(f,U)
 up=po(co(U),p);h=su(up,U);assert all(c%p==0 for c in h)
 zz=su(mu(ev(f,up),po(co(V),-2*p)),one);assert all(c%p==0 for c in zz)
 sq=zero
 for n in range(N):sq=ad(sq,sc(po(zz,n),rat(s.binomial(s.Rational(1,2),n))))
 vn=mu(po(co(V),p),sq);assert all(c%p==0 for c in su(vn,V))
 tiny=[]
 for i in range(5):
  c0=sc(mu(po(U,i),iv(V)),rat(s.Rational(1,2)))
  c1=sc(mu(po(U,i),mu(ev(fp,U),po(V,-3))),rat(-s.Rational(1,4)))
  if i:c1=ad(c1,sc(mu(po(U,i-1),iv(V)),rat(s.Rational(i,2))))
  tiny.append(ad(mu(c0,h),sc(mu(c1,mu(h,h)),rat(s.Rational(1,2)))))
 gs=[]
 for ps in prims:
  G=zero
  for B,e in ps:G=ad(G,mu(ev(B,U),po(V,e)))
  gs.append(G)
 tr=lambda A:(2*A[0]-U1*A[1])%mod
 rhs=[2*(tr(gs[i])-p**S*tr(tiny[i]))%mod for i in range(5)]
 assert all(z%p**S==0 for z in rhs)
 rhs=s.Matrix([z//p**S%p**3 for z in rhs])
 M=s.Matrix.hstack(*map(s.Matrix,mats))
 VEC=(s.eye(5)-M.T).inv()*rhs
 ans=[rat(c)%p**3 for c in VEC]
 assert ans==({2035:[520,1872,377,308,1352],102:[663,962,1339,953,1491]}[a%2197])
 assert ans[0]==(b*(a*481-2158))%2197 and ans[1]==(-b*481)%2197
 print('I',ans,'raw',(-2*ans[2])%p**3,flush=True)
for a in range(p):
 if (a*a+60*a+3249)%p==0:run(a)

# Unramified quadratic coefficient calculation for k=-15.
import sympy as s
p=13; N=22; S=6; mod=p**N; x=s.symbols('x')
def trim(A):
 while len(A)>1 and not A[-1]:A.pop()
 return A
def add(A,B):return trim([((A[i] if i<len(A) else 0)+(B[i] if i<len(B) else 0))%mod for i in range(max(len(A),len(B)))])
def scale(A,c):return trim([a*c%mod for a in A])
def sub(A,B):return add(A,scale(B,-1))
def mul(A,B):
 C=[K(0)]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]=(C[i+j]+a*b)%mod
 return trim(C)
def divrem(A,B):
 A=A[:];Q=[K(0)]*max(1,len(A)-len(B)+1);ib=pow(B[-1],-1,mod)
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
 q=s.Rational(q);return K(int(q.p)*pow(int(q.q),-1,mod)%mod)
class K:
 def __init__(self,a=0,b=0):
  if isinstance(a,K):self.a,self.b=a.a,a.b
  else:self.a,self.b=int(a)%mod,int(b)%mod
 def __add__(self,z):
  z=K(z);return K(self.a+z.a,self.b+z.b)
 __radd__=__add__
 def __neg__(self):return K(-self.a,-self.b)
 def __sub__(self,z):return self+-K(z)
 def __rsub__(self,z):return K(z)+-self
 def __mul__(self,z):
  z=K(z);return K(self.a*z.a-81225*self.b*z.b,self.a*z.b+self.b*z.a-1020*self.b*z.b)
 __rmul__=__mul__
 def co(self):return K(self.a-1020*self.b,-self.b)
 def __pow__(self,n,unused=None):
  if n<0:
   C=self.co();d=self*C;assert d.b==0 and d.a%p
   return (C*pow(d.a,-1,mod))**(-n)
  z=K(1);a=self
  while n:
   if n%2:z=z*a
   a=a*a;n//=2
  return z
 def __mod__(self,n):return K(self.a%n,self.b%n)
 def __floordiv__(self,n):return K(self.a//n,self.b//n)
 def __eq__(self,z):
  z=K(z);return self.a==z.a and self.b==z.b
 def __bool__(self):return self.a!=0 or self.b!=0
 def __repr__(self):return str((self.a,self.b))
def solve(A,B):
 A=[list(row)+[b] for row,b in zip(A,B)];n=len(B)
 for j in range(n):
  ix=next(i for i in range(j,n) if (A[i][j].a%p or A[i][j].b%p))
  A[j],A[ix]=A[ix],A[j];iv=A[j][j]**-1;A[j]=[c*iv for c in A[j]]
  for i in range(n):
   if i!=j:
    z=A[i][j];A[i]=[a-z*b for a,b in zip(A[i],A[j])]
 return [row[-1] for row in A]
def run():
 a=K(0,1);b=-13590*a-1564650
 f0=[K(0)]+[K(s.Poly(s.expand((x+15)*(x-75)*(x*x+1020*x-61875)),x).nth(i)) for i in range(5)]
 f=[K(0)]*7
 for j,c in enumerate(f0):
  for l in range(j+1):f[6-l]+=c*int(s.binomial(j,l))*a**(j-l)*b**-2
 assert f[6]==1
 fp=deriv(f)
 cols=[rem(mul([K(0)]*j+[K(1)],fp),f) for j in range(6)]
 A=[[cols[j][i] if i<len(cols[j]) else K(0) for j in range(6)] for i in range(6)]
 inv=solve(A,[K(1)]+[K(0)]*5)
 assert rem(mul(fp,inv),f)==[1]
 fs=[c.co() for c in f];fps=deriv(fs)
 fx=[K(0)]*(6*p+1)
 for i,c in enumerate(fs):fx[i*p]=c
 R=sub(fx,pw(f,p));assert all(c%p==0 for c in R)
 def reduce(A,m):
  ps=[]
  while m:
   B=rem(mul(A,inv),f);C,rr=divrem(sub(A,mul(B,fp)),f);assert rr==[0]
   Bdiv=pd(B,2*m-1);ps.append((scale(Bdiv,-1),1-2*m))
   A=add(C,scale(deriv(Bdiv),2));m-=1
  while len(A)>5:
   e=len(A)-1;j=e-5;c=pd([A[-1]*pow(f[-1],-1,mod)%mod],2*e-4)[0]
   B=[0]*j+[c];G=add(mul(B,fp),scale(mul(deriv(B),f),2))
   A=sub(A,G);ps.append((B,1))
  return A,ps
 mats=[];prims=[]
 for i in range(5):
  A=[K(0)];ps=[]
  for n in range(4):
   term=scale([K(0)]*(p*(i+1)-1)+pw(R,n),p**(S+1)*rat(s.binomial(-s.Rational(1,2),n)))
   out,pp=reduce(term,6+13*n);A=add(A,out);ps+=pp
  assert all(c%p**S==0 for c in A)
  mats.append([(A[j] if j<len(A) else K(0))//p**S%p**3 for j in range(5)]);prims.append(ps)
 print('M columns',mats,flush=True)

 # Independent affine-Frobenius structural audit.
 rr=[K(0),K(0),rat(s.Rational(-1,2)),f[5]*rat(s.Rational(1,4)),f[4]*rat(s.Rational(1,4))-f[5]**2*rat(s.Rational(3,16))]
 assert all((sum(rr[j]*mats[i][j] for j in range(5))-13*rr[i].co())%2197==0 for i in range(5))
 B2=[[sum(mats[i][k].co()*mats[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
 power=[[K(int(i==j)) for j in range(5)] for i in range(5)]
 traces=[]
 for nn in range(1,6):
  power=[[sum(power[i][k]*B2[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
  traces.append(sum(power[i][i] for i in range(5)))
 coeff=[K(1)]
 for nn in range(1,6):
  coeff.append(-sum(coeff[nn-i]*traces[i-1] for i in range(1,nn+1))*pow(nn,-1,mod))
 TT=s.symbols('TT');want=s.Poly((TT-169)*(TT**2+6*TT+169)**2,TT).all_coeffs()
 assert all((z-int(w))%2197==0 for z,w in zip(coeff,want))
 print('Independent Frobenius-square characteristic polynomial and residue checks PASS',flush=True)

 def ev(A,z):
  v=K(0)
  for c in A[::-1]:v=v*z+c
  return v
 roots=[K(i,j) for i in range(p) for j in range(p) if (K(i,j)**2-276*K(i,j)+11025)%p==0]
 assert len(roots)==2
 r=roots[0]
 for j in range(N):r=r-(r*r-276*r+11025)*(2*r-276)**-1
 assert r*r-276*r+11025==0
 print('rootR',r%p**3,flush=True)
 gs=[K(0)]*5;ts=[K(0)]*5
 for r in [r,276-r]:
  U=(r-a)**-1;V=(-7794*r+406350)*U**3*b**-1
  assert V*V==ev(f,U)
  up=U**p;h=up-U.co();assert h%p==0
  zz=ev(fs,up)*V**(-2*p)-1;assert zz%p==0
  sq=K(0)
  for n in range(N):sq+=rat(s.binomial(s.Rational(1,2),n))*zz**n
  vn=V**p*sq;assert (vn-V.co())%p==0
  us=U.co();vs=V.co()
  for i in range(5):
   c0=us**i*vs**-1*rat(s.Rational(1,2))
   c1=us**i*ev(fps,us)*vs**-3*rat(-s.Rational(1,4))
   if i:c1+=us**(i-1)*vs**-1*rat(s.Rational(i,2))
   ts[i]+=2*(c0*h+c1*h*h*rat(s.Rational(1,2)))
   for B,e in prims[i]:gs[i]+=2*ev(B,U)*V**e
 rhs=[gs[i]-p**S*ts[i] for i in range(5)]
 assert all(z%p**S==0 for z in rhs)
 rhs=[z//p**S%p**3 for z in rhs]
 A=[[K(int(i==j))-sum(mats[i][k].co()*mats[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
 B=[sum(mats[i][k].co()*rhs[k] for k in range(5))+rhs[i].co() for i in range(5)]
 ans=solve(A,B);ans=[z%p**3 for z in ans]
 assert all((ans[i].co()-sum(mats[i][j]*ans[j] for j in range(5))-rhs[i])%p**3==0 for i in range(5))
 assert ans[0]==(b*(a*1963-650))%p**3
 assert ans[1]==(-b*1963)%p**3
 raw=-2*ans[2];tr=raw+raw.co()
 assert tr%2197==K(1716)
 print('I',ans,'raw one source',raw%p**3,'trace',tr%p**3,flush=True)
run()

assert (962-1409*481-789*2158)%2197==2119
assert (1716-1301*1963-390*650)%2197==2119
assert (2119+871)%2197==793
assert (2119+390)%2197==312
print("Both remaining Q/R local and global mixed heights PASS")
