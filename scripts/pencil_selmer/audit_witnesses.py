import json
from pathlib import Path
import sympy as s
x=s.symbols('x');F=x**5+3478005*x**4-38594700111069*x**3-14431460273343469305*x**2+47793145649188822875249984*x+22058431254036587705175782036640
out=[]
for p,c in [(61,s.Integer(2657)),(773,s.Integer(205278)),(3087467,s.Integer(48537585605194255596273016)),(2,s.Rational(5,4)),(2,s.Integer(2)),(2,s.Integer(885))]:
 q=F.subs(x,c);v=int(s.ntheory.factor_.multiplicity(p,abs(q.p))-s.ntheory.factor_.multiplicity(p,q.q));unit=q/(s.Integer(p)**v);m=8 if p==2 else p;res=int(unit.p)*pow(int(unit.q),-1,m)%m
 assert v%2==0 and (res==1 if p==2 else pow(res,(p-1)//2,p)==1)
 out.append({'p':p,'c':str(c),'valuation':v,'normalized_residue':res,'modulus':m})
for p,h in [(35498020282877,x*x+401107556428946897390988312*x+895917242855804582442647861843282088065835134709300681),(5,x*x+85*x+23),(7,x*x+22*x+78271),(3087467,x*x+2190851906626*x+19753949860961491482589978)]:
 factors=s.factor_list(h,modulus=p)[1];assert all(mult==1 for fac,mult in factors)
 r=s.rem(F,h,x);entry={'p':p,'h':str(h),'F_remainder':str(r),'components':[]}
 if len(factors)==1:
  assert s.degree(factors[0][0],x)==2
  v=min(int(s.ntheory.factor_.multiplicity(p,abs(int(c)))) for c in s.Poly(r,x).all_coeffs() if c)
  assert v%2==0;rr=s.Poly(r/(s.Integer(p)**v),x,modulus=p);hh=s.Poly(h,x,modulus=p)
  n=(p*p-1)//2;acc=s.Poly(1,x,modulus=p)
  while n:
   if n&1:acc=(acc*rr).rem(hh)
   rr=(rr*rr).rem(hh);n//=2
  assert acc==s.Poly(1,x,modulus=p)
  entry['components'].append({'degree':2,'valuation':v,'euler_test':'1'})
 else:
  assert len(factors)==2
  for fac,mult in factors:
   cc=s.Poly(fac,x).all_coeffs();root=(-int(cc[1])*pow(int(cc[0]),-1,p))%p;res=int(r.subs(x,root))%p
   assert res!=0 and pow(res,(p-1)//2,p)==1
   entry['components'].append({'degree':1,'root_mod_p':root,'valuation':0,'unit_residue':res,'euler_test':'1'})
 out.append(entry)
print(json.dumps(out,indent=2))

print('All rational and quadratic divisor witnesses are independently locally valid.')
