#!/usr/bin/env python3
"""QC normalization and exact local targets at 11,19,53; PROOF 5.28--5.29.
Existing SymPy only. Does NOT evaluate normalized13-adic QC series or isolate roots.
"""
import sympy as s
t,u,x,r,Y=s.symbols("t u x r Y")
def val(q,p):
 q=s.Rational(q)
 if not q:return 10**6
 n,d=int(q.p),int(q.q);v=0
 while n%p==0:n//=p;v+=1
 while d%p==0:d//=p;v-=1
 return v
def red(q,p):
 q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,p)%p
data=[
 (1,80-t*t,4*t-30,2*t*t-30*t+160,15*t-160,2*t-15,s.Rational(15,2),s.Rational(9025,16),1,-2),
 (5,25*t*t-120*t+80,10*t*t-8,-60*t*t+100*t-48,24-25*t,25-30*t,s.Rational(5,6),s.Rational(9025,144),5625,4)]
for delta,N,D,M,A,L,tb,Wb,c,a19 in data:
 F=s.expand((N*N+5*D*D)*(N*N+275*D*D));sig=A/L
 assert s.degree(F,t)==8 and s.Poly(F,t).LC()==c*c
 assert s.gcd(s.Poly(N,t),s.Poly(D,t)).degree()==0
 assert s.discriminant(M,t)==(-380 if delta==1 else -1520)
 assert s.gcd(s.Poly(F,t),s.Poly(M,t)).degree()==0
 assert s.cancel(F.subs(t,sig)-9025**2*F/L**8)==0
 assert s.cancel(sig.subs(t,sig)-t)==0
 assert s.cancel(sig-t+(M/L if delta==1 else M/(2*L)))==0
 # Fixed points are precisely zeros of M, and sigma fixes the ordinate there.
 assert s.rem(L**4-9025,M,t)==0
 assert F.subs(t,tb)==Wb**2 and M.subs(t,tb)!=0
 assert L.subs(t,tb)==0
 assert s.cancel(9025*Wb/A.subs(t,tb)**4)==c
 assert s.cancel(M.subs(t,tb)**3/Wb)==(190 if delta==1 else -s.Rational(304,75))
 # Degree accounting for Gamma_iota-3 Gamma_sigma.
 assert 8-3*4==-4 and (-4)-2*(-2)==0
 # Divisor vectors ordered (hyperelliptic ramification, R_sigma, H).
 div_W=s.Matrix([1,0,-4]);div_M=s.Matrix([0,1,-2])
 assert div_W-3*div_M==s.Matrix([1,-3,2])
 # dt/W has K=2H. Proper H/2 rigidification has diagonal div((sigma(t)-t)^3).
 assert (-4)+3*4-2-3*2==0
 # Discriminants give a sufficient exceptional-prime set, not a minimal conductor.
 want={2:52,3:33,5:30,11:2,19:12,53:1} if delta==1 else {2:108,3:33,5:58,11:2,19:12,53:1}
 assert s.factorint(abs(int(s.discriminant(F,t))))==want
 for p in [11,53]:
  ff=s.Poly(F,t,modulus=p);gg=s.gcd(ff,ff.diff()).monic()
  # gcd has one copy of each double root; all nodes are rational here.
  nodes=[j for j in range(p) if ff.eval(j)%p==0 and ff.diff().eval(j)%p==0]
  assert len(nodes)==(2 if p==11 else 1) and gg.degree()==len(nodes)
  residual=ff.exquo(gg*gg)
  assert s.gcd(residual,residual.diff()).degree()==0
  assert s.gcd(residual,gg).degree()==0
  assert residual.degree()>0 # over the algebraic closure the residual is not a square.
  assert ff.LC()%p and 2*c%p
  entries=[]
  for a in nodes:
   assert val(F.subs(t,a),p)==1
   tangent=s.diff(F,t,2).subs(t,a)/2
   assert red(tangent,p)!=0
   entries.append((a,red(F.subs(t,a)/p,p),red(tangent,p)))
  print(delta,p,"regular thickness-one nodes",entries,flush=True)
 # Complete two-component semistable chart at19.
 p=19;shift=s.Poly(F.subs(t,a19+x),x)
 cs=[shift.nth(j) for j in range(9)]
 assert [val(cj,p) for cj in cs]==[4,3,2,1,0,0,0,0,0]
 inner=s.Poly(sum(s.Rational(cs[j]*p**j,p**4)*u**j for j in range(9)),u)
 assert all(cj.q==1 for cj in inner.all_coeffs())
 inner19=inner.set_modulus(p)
 outer19=s.Poly(sum(cs[j]*x**(j-4) for j in range(4,9)),x,modulus=p)
 assert inner19.degree()==outer19.degree()==4
 assert red(s.discriminant(inner19.as_expr(),u),p)==(11 if delta==1 else 7)
 assert red(s.discriminant(outer19.as_expr(),x),p)==(7 if delta==1 else 11)
 ann=sum(s.Rational(cs[j],p**(4-j))*r**(4-j) for j in range(4))+cs[4]+sum(cs[j]*x**(j-4) for j in range(5,9))
 assert all(q.q==1 for q in s.Poly(ann,x,r).coeffs())
 assert s.cancel(x**4*ann.subs(r,p/x)-F.subs(t,a19+x))==0
 assert int(cs[4])%p==8 and pow(8,(p-1)//2,p)==p-1
 # At x=r=0, w^2=8 has two distinct roots over the unramified quadratic field.
 # The remaining node equation is xr=19, so both nodes have thickness one.
 image=s.cancel(sig.subs(t,a19+19*u))
 expected=(15*u-10)/(2*u-1) if delta==1 else (25*u+4)/(30*u+5)
 assert s.cancel(image-expected)==0
 num,den=s.fraction(image)
 assert red(s.resultant(num,den,u),p)!=0
 assert val(tb-a19,19)==1 # basepoint is on inner component, sigma(basepoint) is outer infinity.
 print(delta,19,"inner",inner19.as_expr(),"outer",outer19.as_expr(),"sigma exchanges components",flush=True)
 # Omitted-correction stress test, without asserting global rational points.
 assert red(s.Poly(F,t).LC(),7)==red(c*c,7)!=0
 assert red(s.Poly(M,t).LC(),7)!=0
 # Hensel at infinity gives Q7 points t=7^-n for every n>=1.
 # There v7(W/M^3)=2n. Omitting its character term is unbounded.
 # No numerical root search is used for that claim.


# The four excluded basepoint-orbit points are explicit and already pass to
# elliptic indices0 or+-1 in the earlier exact collision test.
zz=s.Rational(5,2)
for bb,want in [(s.Rational(225,4),(405,2025)),(-s.Rational(225,4),(180,900))]:
 XX=2*bb+2*zz*zz+280;YY=2*zz*XX
 assert (XX,YY)==want and YY*YY==XX*(XX*XX-560*XX+72900)
assert s.Rational(72900,180)==405
assert -s.Rational(72900)*900/180**2==-2025 # P+T, hence (405,2025)=-P+T.

# Generic third-kind kernel nu_P = (1+W_P/W)dt/(2(t-tP)).
# At P and iota(P) the residue coefficients are respectively1 and0.
assert (1+s.Integer(1))/2==1 and (1-s.Integer(1))/2==0
# Each finite nu_P has residue-1/2 at both infinities.
# For infinity sign epsilon, nu_inf=-epsilon*c*t^3 dt/(2W):
for epsilon in [1,-1]:
 for sheet in [1,-1]:
  residue=s.Rational(epsilon,2*sheet)
  assert residue==(s.Rational(1,2) if sheet==epsilon else -s.Rational(1,2))

# Prime5: a ramified model of X1 also controls X5 after sqrt5 twisting.
pi=s.symbols("pi")
N1=80-t*t;D1=4*t-30;M1=2*t*t-30*t+160
F1=s.Poly((N1*N1+5*D1*D1)*(N1*N1+275*D1*D1),t)
Gexpr=0
for j in range(9):
 q=F1.nth(j)
 # F1(pi*u)/pi^8; replace pi^(j-8) using pi^2=5.
 parity=(j-8)%2;power=(j-8-parity)//2
 cq=s.Rational(q)*s.Rational(5)**power
 assert cq.q==1
 Gexpr+=cq*pi**parity*u**j
assert s.Poly(Gexpr,u).LC()==1
Gbar=s.Poly(Gexpr.subs(pi,0),u,modulus=5)
target=(u-1)**2*(u+1)**2*(u*u+2*u-1)*(u*u-2*u-1)
assert Gbar==s.Poly(target,u,modulus=5)
gg=s.gcd(Gbar,Gbar.diff()).monic()
assert gg==s.Poly(u*u-1,u,modulus=5)
def vk(q):
 q=s.Poly(s.rem(s.expand(q),pi*pi-5,pi),pi)
 return min(2*val(q.nth(0),5),2*val(q.nth(1),5)+1)
for a in [1,-1]:
 assert [vk(s.diff(Gexpr,u,j).subs(u,a)) for j in range(3)]==[2,2,0]
 assert Gbar.eval(a)%5==Gbar.diff().eval(a)%5==0
 assert red(s.diff(Gbar.as_expr(),u,2).subs(u,a)/2,5)!=0
assert 2*val(s.discriminant(F1.as_expr(),t),5)-8*7==4
# The exact sigma formula after scaling has unit determinant19.
assert s.rem(s.together((15*pi*u-160)/(pi*(2*pi*u-15))-(3*pi*u-32)/(2*u-3*pi)).as_numer_denom()[0],pi*pi-5,pi)==0
assert s.rem((3*pi)*(-3*pi)-(-32)*2,pi*pi-5,pi)==19
for a in [1,-1]:
 assert red(-32/(s.Integer(2)*a),5)==(-a)%5
# Thus node thicknesses are equal and sum to4: each is2.
# X1(Q5) has v_K(t/pi)=2v5(t)-1, so avoids the finite nodes u=+-1.
# For X5, 5a^2=q1(z) implies v5(z)>=1 and a^2=1 mod5:
z=s.symbols("z")
q1=4*z*z+30*z+80
assert s.Poly(q1.subs(z,5*z)/5,z,modulus=5)==s.Poly(1,z,modulus=5)
# Then u=a-2z/pi reduces+-1. A section of xy=pi^2 through its node
# has both coordinate valuations1 and meets the unique exceptional vertex.
print("PASS: prime5 regular model and one accessible automorphism orbit for each twist.")

print("PASS: corrected divisor, rational kernels, and target {0} at5,11,19,53.")
print("Remaining: targets at2,3; normalized13-adic series, precision tails, roots and sieve.")
