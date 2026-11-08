#!/usr/bin/env python3
"""Sixth-mode QC root certificate; see SIXTH_MODE_QC_EXCLUSION.md.
Uses the explicitly listed previously certified PR inputs. Existing SymPy only.
Does not classify all rational curve points or solve general size-eleven PTE.
"""
import sympy as s
import json
x,a,h,u=s.symbols('x a h u'); Q=s.Rational
rows=[]
computed_qc={}
def mod(q):
 q=Q(q); assert q.q%13; return int(q.p)*pow(int(q.q),-1,2197)%2197
for delta,N,D,M,c,A,L,b0,W0,lamrho in [
 (1,80-x*x,4*x-30,2*x*x-30*x+160,1,15*x-160,2*x-15,Q(15,2),Q(9025,16),[289,925,590]),
 (5,25*x*x-120*x+80,10*x*x-8,-60*x*x+100*x-48,5625,24-25*x,25-30*x,Q(5,6),Q(9025,144),[427,1168,628])]:
 f=s.Poly(s.expand((N*N+5*D*D)*(N*N+275*D*D))/c**2,x,domain=s.QQ); fp=f.diff(); fe=f.as_expr(); dp=s.Poly(D,x,domain=s.QQ)
 K=s.cancel((s.diff(N,x)*D-N*s.diff(D,x))/M)
 forms=[(K*M*N*N/D**2,2),(K*(N**3+Q(15,4)*N*N*D)/D**2,2),(K*(N**4+Q(4125,32)*N*N*D*D+525*N*D**3)/D**3,3)]
 def reduce(R,m):
  P=s.Poly(s.cancel(R*D**m/c),x,domain=s.QQ)
  while m>1:
   C=(-P*s.invert(dp.diff()*f*(m-1),dp)).rem(dp)
   exact=C.diff()*dp*f-C*dp.diff()*f*(m-1)+C*dp*fp/2
   P=(P-exact).exquo(dp);m-=1
  P=P.exquo(dp)
  while P.degree()>=7:
   e=P.degree();C=s.Poly(x**(e-7),x)
   P-=(C.diff()*f+C*fp/2)*P.LC()/s.Integer(e-3)
  return s.Matrix([2*P.nth(i) for i in range(7)])
 rev=sum(f.nth(8-j)*h**j for j in range(9))
 rs=s.series(rev**(-Q(1,2)),h,0,4).removeO()
 B=s.zeros(7)
 for i in range(3): B[i,i]=1; B[i+4,i+3]=1; B[3,i+3]=-rs.coeff(h,i+1)
 B[3,6]=1
 W=[B.inv()*reduce(R,m) for R,m in forms]
 assert all(v[6]==0 for v in W)
 S=s.Matrix.hstack(s.eye(6)[:,:3],*[v[:6,0] for v in W])
 E=S.inv()[:3,:]*B.inv()[:6,:]
 assert all(Q(z).q%13 for z in E)
 def project(poly):
  vec=s.Matrix([s.expand(poly).coeff(x,i) for i in range(7)])
  assert s.expand((B.inv()*vec)[6])==0
  return (E*vec).applyfunc(s.expand)
 # Source derivative modulo exact differential, with y(a)^2=f(a).
 R=s.cancel(((fp.as_expr().subs(x,a)+fp.as_expr())*(x-a)-2*(fe-fe.subs(x,a)))/(x-a)**2)
 assert s.denom(R)==1 or not s.denom(R).has(x,a)
 assert s.degree(R,x)<=6
 # Verify raw source derivative minus derivative_x[-1/(2(x-a))-Y/(2y(a)(x-a))].
 identity=(fp.as_expr().subs(x,a)/(x-a)+2*fe.subs(x,a)/(x-a)**2 + fp.as_expr()/(x-a)-2*fe/(x-a)**2-R)
 assert s.cancel(identity)==0
 PR=project(R)
 y0=W0/c; f0=fe.subs(x,b0); assert y0*y0==f0
 first=PR.subs(a,b0)/(2*y0)
 second=(PR.diff(a).subs(a,b0)-PR.subs(a,b0)*fp.as_expr().subs(x,b0)/(2*f0))/(4*y0)
 # Obtain the finite-basepoint projection from sigma pullback of nu_infinity+.
 sig=A/L; scale=Q(9025)/L**4
 assert s.cancel((N/D).subs(x,sig)-N/D)==0
 assert s.cancel(scale/D.subs(x,sig)**2-1/D**2)==0
 # sigma fixes (z,b); iota changes b's sign, so all qi(z) are unchanged.
 sh=[s.cancel(sig**i*s.diff(sig,x)/scale) for i in range(3)]
 H=s.Matrix([[s.expand(q).coeff(x,j) for j in range(3)] for q in sh])
 difference=s.cancel(y0/(x-b0)+sig**3*s.diff(sig,x)/scale)
 assert s.denom(difference).is_number and s.degree(difference,x)<=2
 l0=-H.T*s.Matrix(lamrho)+s.Matrix([s.expand(difference).coeff(x,j) for j in range(3)])
 # Infinity-regular source derivative; v(h)^2=rev(h).
 num=-s.diff(rev,h)*x**3*(1-h*x)-2*rev*x**4+h**3*fp.as_expr()*(1-h*x)+2*h**4*fe
 RI=s.cancel(num/(1-h*x)**2)
 assert s.denom(RI).is_number and s.degree(RI,x)<=6
 PI=project(RI)
 ifirst=PI.subs(h,0)/2
 isecond=(PI.diff(h).subs(h,0)-PI.subs(h,0)*s.diff(rev,h).subs(h,0)/2)/4
 # Cup pairing in residue-zero eta basis; explicit Laurent residue convention.
 sr=s.series(rev**(-Q(1,2)),h,0,9).removeO()
 etas=[sum(B[j,i]*x**j for j in range(7)) for i in range(6)]
 la=[s.series(-h**2*pp.subs(x,1/h)*sr/2,h,0,7).removeO().expand() for pp in etas]
 def prim(pp):
  ans=0
  for term in s.Add.make_args(pp):
   co,ex=term.as_coeff_exponent(h)
   assert ex!=-1
   ans+=co*h**(ex+1)/(ex+1)
  return ans
 cup=s.Matrix(6,6,lambda i,j:2*s.expand(prim(la[i])*la[j]).coeff(h,-1))
 assert cup==-cup.T
 vR=s.Matrix([s.expand(R).coeff(x,i) for i in range(7)])
 dd=B.inv()*vR
 # Cup(Psi dnu,hol) = -hol(source), fixing the Abel-Jacobi sign.
 test=(dd[:6,0].T*cup[:,:3]).applyfunc(s.simplify)
 print("cup derivative check",test); assert test==s.Matrix([[-1,-a,-a*a]])
 psirho=s.Matrix([1433,392,1539,1586,728,1248] if delta==1 else [1194,1844,368,1690,533,1898])
 # log centered B = -cup(Psi nu_B,hol); hol correction pairs to zero.
 LB=(psirho.T*cup[:,:3]*H.T).T
 assert all(mod(z)%13==0 for z in LB)
 k=-15 if delta==1 else -3
 rawhol=[(s.diff(N,x)*D-N*s.diff(D,x))/c,2*K*D/(k*c),2*K*N/c]
 G=s.Matrix([[s.expand(pp).coeff(x,j) for j in range(3)] for pp in rawhol])
 LC=G*LB
 vb=s.Matrix([b0**i/(2*y0) for i in range(3)])
 vbp=s.Matrix([(i*b0**(i-1) if i else 0)/(2*y0)-b0**i*fp.as_expr().subs(x,b0)/(4*y0*f0) for i in range(3)])
 d1=G*vb;d2=G*vbp/2
 Aglob=s.Matrix([[37,131],[131,135]] if delta==1 else [[90,114],[114,36]])
 J=s.diag(144,-Aglob/2)
 glob1=2*(LC.T*J*d1)[0]/13
 glob2=((d1.T*J*d1)[0]+2*(LC.T*J*d2)[0])/13
 glob3=2*(d1.T*J*d2)[0]/13
 print(delta,'LOGCENTER',[mod(z) for z in LC])
 # Actual normalized local-height tiny integral on the basepoint disc.
 # F(r,q) is the normalized kernel integrated from r=0 to r=q.
 f1=fp.as_expr().subs(x,b0); f2=s.diff(fe,x,2).subs(x,b0)
 m0=M.subs(x,b0); m1=s.diff(M,x).subs(x,b0); m2=s.diff(M,x,2).subs(x,b0)
 dg=f1/(2*f0)-3*m1/m0
 dgr=f2/(2*f0)-f1*f1/(2*f0*f0)-3*(m2/m0-m1*m1/m0**2)
 ni=-b0**3/(2*y0)
 nir=-3*b0*b0/(2*y0)+b0**3*f1/(4*y0*f0)
 no=f1/(4*f0)
 nor=f2/(8*f0)-3*f1*f1/(16*f0*f0)
 noq=f2/(8*f0)-f1*f1/(16*f0*f0)
 hp=s.diff(L,x).subs(x,b0)/A.subs(x,b0)
 viq=s.diff(rev,h).subs(h,0)*hp/2
 niq=-hp/2-(viq*b0**3+hp*b0**4)/(2*y0)
 basis=s.Matrix([b0**i/(2*y0) for i in range(3)])
 basisr=s.Matrix([(i*b0**(i-1) if i else 0)/(2*y0)-b0**i*f1/(4*y0*f0) for i in range(3)])
 lt=-6*s.Matrix(lamrho)+2*l0
 ltq=3*ifirst*hp+first
 hv1=dg+6*ni-2*no-(lt.T*basis)[0]
 hv2=(dgr+6*nir-2*nor-(lt.T*basisr)[0])/2+3*niq-noq-(ltq.T*basis)[0]
 print(delta,'LOCAL13_HEIGHT_AT_B linear,quadratic',[mod(hv1),mod(hv2)])
 computed_qc[delta]=[mod(hv1-glob1)%169,mod(13*(hv2-glob2))%169,mod(-169*glob3)%169]
 print(delta,'RHO_ON_u_DIV13',computed_qc[delta])
 rows.append({'delta':delta,'local_height_basepoint_disc_linear_quadratic_mod_2197':[mod(hv1),mod(hv2)]})
 # Exact rational validation of the infinity derivative reduction.
 # Multiply difference by 4v; d[-x/(2(1-hx))-h^3Y/(2v(1-hx))].
 # Independent raw derivative minus exact primitive check, split rational/Y.
 vv,YY=s.symbols('vv YY', nonzero=True)
 kernel=-h/(2*(1-h*x))-vv*x**3/(2*YY*(1-h*x))
 raw=s.diff(kernel,h)+s.diff(kernel,vv)*s.diff(rev,h)/(2*vv)
 primitive=-x/(2*(1-h*x))-h**3*YY/(2*vv*(1-h*x))
 exact=s.diff(primitive,x)+s.diff(primitive,YY)*fp.as_expr()/(2*YY)
 guard=s.together((raw-exact-RI/(4*vv*YY))*4*vv*YY)
 gn,gd=s.fraction(guard)
 gn=s.Poly(gn,YY).rem(s.Poly(YY**2-fe,YY)).as_expr()
 gn=s.Poly(gn,vv).rem(s.Poly(vv**2-rev,vv)).as_expr()
 assert s.cancel(gn)==0

 # Verify residue zero directly at infinity independently of coordinate conversion.
 for poly in [R,RI]:
  assert s.expand(poly).coeff(x,3)+sum(s.expand(poly).coeff(x,j)*rs.coeff(h,j-3) for j in range(4,7))==0
 for chart,const,one,two in [('finite_B',l0,first,second),('infinity',-s.Matrix(lamrho),ifirst,isecond)]:
  coefficients=[[mod(z) for z in col] for col in [const,one,two]]
  rows.append({'delta':delta,'positive_sheet_chart':chart,'constant_linear_quadratic_vectors_mod_2197':coefficients})
  print(delta,chart,coefficients)
 # Integrality for all derivatives: polynomial coefficients, inverse square root,
 # projection E and basepoint unit square root are 13-integral.
 for poly in [R,RI]:
  assert all(Q(z).q%13 for z in s.Poly(poly,x,a,h).coeffs())
 assert mod(y0)%13 and mod(fe.subs(x,b0))%13
print('PASS: exact derivative reductions, sigma basepoint transport, residue checks, integral projection.')
print('Negative-sheet coefficients are the negatives of the displayed vectors.')
print('For source parameter u in 13 Z_13, omitted integrated powers are O(13^3).')

# Exhaustive residue-root certificate for the EXACT analytic QC functions.
# See SIXTH_MODE_QC_EXCLUSION.md for the integral analytic-tail/Hensel proof.
root_rows=[]
for delta,coeffs,targets,wanted,base_n in [
 (1,[149,136,0],[Q(0),-Q(1,6),-Q(2,3)],[[0,12],[2,10],[]],0),
 (5,[129,110,117],[Q(0),Q(1,2),Q(2,3),-Q(4,3)],[[0,11],[2,9],[3,8],[5,6]],12)]:
 assert computed_qc[delta]==coeffs
 for target,expected in zip(targets,wanted):
  target169=mod(target*161)%169
  def value(u):return sum(c*u**(i+1) for i,c in enumerate(coeffs))
  roots=[u for u in range(13) if (value(u)-target169)%13==0]
  assert roots==expected
  derivatives=[sum((i+1)*c*u**i for i,c in enumerate(coeffs))%13 for u in roots]
  assert all(derivatives)
  lifts=[u for u in range(169) if (value(u)-target169)%169==0]
  assert len(lifts)==len(roots) and sorted(u%13 for u in lifts)==roots
  surviving=[u for u in roots if (base_n+11*u)%13==base_n]
  assert surviving==([0] if target==0 else [])
  root_rows.append({'delta':delta,'target_multiple_log3':str(target),'roots_mod13':roots,'derivatives_mod13':derivatives,'roots_mod169':lifts,'sieve_survivors_mod13':surviving})
  print('ROOTS',root_rows[-1])
assert sum(len(r['roots_mod13']) for r in root_rows)==12
assert sum(len(r['sieve_survivors_mod13']) for r in root_rows)==2
print('PASS: 12 roots on two B discs, transported to48 roots; elliptic sieve leaves only the8 known orbit points.')
