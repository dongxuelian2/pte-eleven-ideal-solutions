#!/usr/bin/env python3
"""Three certified normalized mixed heights; PROOF 5.25.
Requires adjacent check_coleman_basis13.py and existing SymPy.
Does not complete the height matrices or compute QC roots.
"""
import sympy as s
from check_coleman_basis13 import Ring,res,p,N,MOD,scale,x
for k,pa,pb,original in [
(-15,s.Rational(225,4),s.Rational(320625,32),[(x*x+1020*x+81225,-13590*x-1564650),(x*x-276*x+11025,406350-7794*x)]),
(-3,s.Integer(33),s.Integer(10692),[(x*x+65*x+528,-1076*x-12804)])]:
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
  expected={(-15,"Q"):[416,2002,1222,323,714],(-15,"R"):[468,234,910,1142,762],(-3,"R"):[2015,325,1079,1213,1872]}
  assert list(I)==expected[k,label]
  old={(-15,"Q"):[195,117,182,1343],(-15,"R"):[1963,650,14,22],(-3,"R"):[481,2158,1297,214]}[k,label]
  assert (int(I[0])-res(pb*(pa*old[0]-old[1])))%2197==0
  assert (int(I[1])+res(pb*old[0]))%2197==0
  lam={-15:[1545,1032],-3:[2021,1085]}[k]
  raw=int(-2*I[2])%2197
  local=(raw-lam[0]*old[0]-lam[1]*old[1])%2197
  assert local=={(-15,"Q"):1755,(-15,"R"):884,(-3,"R"):1781}[k,label]
  def log(q):
   zz=(pow(q,12,13**5)-1)%13**5
   return res(sum((s.Rational((-1)**(j+1)*zz**j,12*j) for j in range(1,6)),s.Rational(0)))%2197
  away=(-2*(log(17)+log(593)) if label=="Q" else (-2*log(5) if k==-15 else 2*log(19)))%2197
  glob=(local+away)%2197
  assert glob=={(-15,"Q"):1950,(-15,"R"):1248,(-3,"R"):1274}[k,label]
  print(k,label,"even integrals",list(I),"raw",raw,"local height",local,"global mixed",glob,flush=True)
 coeff=f.all_coeffs()
 rr=s.Matrix([[0,0,-s.Rational(1,2),coeff[1]/4,coeff[2]/4-3*coeff[1]**2/16]])
 assert (rr*M-13*rr).applyfunc(lambda a:res(a)%2197)==s.zeros(1,5)
 print("Even Frobenius and infinity-residue checks PASS",k,flush=True)

assert (1248-1950)*pow(4,-1,2197)%2197==923
assert 1950*pow(4,-1,2197)%2197==1586
assert 1274*pow(4,-1,2197)%2197==1417
print("Partial global entries certified; complete matrices and QC roots NOT computed.")
