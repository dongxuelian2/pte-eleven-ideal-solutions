"""Exact polynomial certificate for all eighteen normalized shift-chain types.
See PROOF.md, Theorem 14. The finite arithmetic checks verify the explicitly
proved reductions; Faltings and finite fibers give finiteness. No point search,
network access, generated output, GRH or conjecture about Sha is used.
"""
from pathlib import Path
import runpy
_script_dir = Path(__file__).resolve().parent
for _certificate in ("check_shift_chains.py", "check_fifth_mode.py", "check_sixth_mode.py"):
    runpy.run_path(str(_script_dir / _certificate))
print("Previously proved types 1,2,9,10,14,17 verified.")

import sympy as s
from itertools import combinations_with_replacement
x,z,b,c,S,T,U,V=s.symbols('x z b c S T U V')
Q=s.Rational

def eq(a,b):
    assert s.cancel(a-b)==0

def offsets(L):
    return [Q(2*j-L+1,2) for j in range(L)]

def moment(L,k,v):
    return s.expand(sum(s.binomial(2*k,2*j)*v**(k-j)*sum(d**(2*j) for d in offsets(L)) for j in range(k+1)))

def central(m,k):
    return sum(d**(2*k) for d in offsets(2*m))/2

def pair_power(k,sm,pr):
    ps=[s.Integer(2),sm]
    for j in range(2,k+1):
        ps.append(s.expand(sm*ps[-1]-pr*ps[-2]))
    return ps[k]

def pair_moment(L,k,sm,pr):
    poly=s.Poly(moment(L,k,x),x)
    return s.expand(sum(co*pair_power(mon[0],sm,pr) for mon,co in poly.terms()))

def elementary(p):
    e=[s.Integer(1)]
    for k in range(1,len(p)):
        e.append(s.expand(sum((-1)**(i-1)*e[k-i]*p[i] for i in range(1,k+1))/k))
    return e

def root_poly(p):
    e=elementary(p);n=len(e)-1
    return s.expand(sum((-1)**i*e[i]*x**(n-i) for i in range(n+1)))

def fourth_residual(p):
    e=elementary(p[:4])
    return s.factor(p[4]-e[1]*p[3]+e[2]*p[2]-e[3]*p[1])

def curve(label,f,var,degree,genus):
    p=s.Poly(f,var,domain=s.QQ)
    assert p.degree()==degree
    assert s.gcd(p,p.diff()).degree()==0
    assert (degree-1)//2==genus>=2
    print('type',label,': squarefree degree',degree,'genus',genus)

def central_collision(m,L,center):
    assert set(offsets(2*m)) & {center+d for d in offsets(L)}

# The formulas below use one-sided chain moments; direct expansion checks
# the convention independently of the binomial implementation.
for L in range(1,6):
    for k in range(1,5):
        direct=s.Poly(sum((z+d)**(2*k) for d in offsets(L)),z)
        assert all(mon[0]%2==0 for mon,co in direct.terms())
        converted=sum(co*x**(mon[0]//2) for mon,co in direct.terms())
        eq(moment(L,k,x),converted)
patterns=[]
for m in range(1,6):
    for n in range(6):
        for aa in combinations_with_replacement(range(1,6),n):
            for bb in combinations_with_replacement(range(1,6),5-n):
                if sum(aa)==5-m and sum(bb)==5:
                    patterns.append((m,aa,bb))
assert len(patterns)==18
for i,p in {5:(1,(2,2),(1,2,2)),8:(1,(1,1,1,1),(5,)),12:(2,(1,1,1),(1,4)),13:(2,(1,1,1),(2,3)),15:(3,(1,1),(1,1,3)),16:(3,(1,1),(1,2,2))}.items():
    assert patterns[i-1]==p

# Type 5.
E=[s.expand(central(1,k)+pair_moment(2,k,S,T)-pair_moment(2,k,U,V)-c**k) for k in range(1,5)]
uv={U:(8*S-4*c+1)/8,V:(-64*S*c+16*S+128*T+48*c*c-56*c+11)/128}
for i in [0,1]:eq(E[i].subs(uv,simultaneous=True),0)
e3=s.factor(E[2].subs(uv,simultaneous=True));e4=s.factor(E[3].subs(uv,simultaneous=True))
eq(e3,-3*(4*c-1)*(-96*S*c+88*S+128*T+80*c*c-200*c+173)/512)
ts=(96*S*c-88*S-80*c*c+200*c-173)/128
eq(e4.subs(T,ts),(4*c-1)*(4*c-9)*(4*c-25)*(16*S-20*c+85)/4096)
eq((x*x-S*x+ts).subs(c,Q(9,4)),(x-1)*(x-S+1))
eq((x*x-S*x+ts).subs(c,Q(25,4)),(x-4)*(x-S+4))
central_collision(1,1,Q(1,2));central_collision(1,2,s.Integer(1))
assert Q(5,2) in {s.Integer(2)+d for d in offsets(2)}
vals={S:(20*c-85)/16,T:(80*c*c-840*c+589)/256,U:(12*c-83)/16,V:(4*c-49)*(4*c-9)/256}
for e in E:eq(e.subs(vals,simultaneous=True),0)
f5=c*(80*c*c-840*c+589)*(4*c-49)*(4*c-9)
eq(256**2*c*vals[T]*vals[V],f5)
curve(5,f5,c,5,2)

# Type 8.
p=[None]+[moment(5,k,b)-central(1,k) for k in range(1,5)]
K=root_poly(p)
h8=26214400*b**10-904396800*b**9+13285785600*b**8-130282291200*b**7+1026147491840*b**6-6126693150720*b**5+21780357808640*b**4-27633876692480*b**3+12983200818768*b**2+1848981232152*b+6457704345
eq(s.discriminant(K,x),5*(4*b-1)**2*h8/2**24)
central_collision(1,5,Q(1,2))
curve(8,5*h8,b,10,4)
product8=5*(256*b**4-1152*b**3+1824*b**2-388*b+405)
eq(256*K.subs(x,0),product8)
v8=s.symbols('v8')
B8=s.prod(x-(v8+j)**2 for j in range(-2,3))
eq((x-Q(1,4))*K.subs(b,v8**2), B8-B8.subs(x,Q(1,4)))


# Types 12 and 13.
p=[None]+[moment(4,k,b)+c**k-central(2,k) for k in range(1,5)]
eq(fourth_residual(p),-4*b*(b-4)*(b-1)*(b+4*c-10))
for center in [0,1,2]:central_collision(2,4,s.Integer(center))
K=root_poly(p[:4]).subs(c,(10-b)/4)
h12=50*b**6-1650*b**5+26675*b**4-187600*b**3+49540*b*b-34390*b-81
eq(s.discriminant(K,x),-h12/64)
curve(12,-h12,b,6,2)
p=[None]+[moment(2,k,b)+moment(3,k,c)-central(2,k) for k in range(1,5)]
eq(fourth_residual(p),-(b-1)*(4*c-9)*(4*c-1)*(6*b+4*c-25)/8)
central_collision(2,2,s.Integer(1))
for center in [Q(1,2),Q(3,2)]:central_collision(2,3,center)
K=root_poly(p[:4]).subs(b,(25-4*c)/6)
h13_product=-240*c**3+1600*c*c-260*c+475
eq(-K.subs(x,0),h13_product/144)
h13=921600*c**6-24832000*c**5+192608000*c**4-644627200*c**3+936509680*c*c-438571880*c-82872167
eq(s.discriminant(K,x),-h13/62208)
assert 432**2==3*62208
curve(13,-3*h13,c,6,2)

# Type 15.
E=[s.expand(central(3,k)+pair_power(k,S,T)-pair_power(k,U,V)-moment(3,k,c)) for k in range(1,5)]
uv={U:S-3*c+Q(27,4),V:(-48*S*c+108*S+16*T+96*c*c-228*c+27)/16}
for i in [0,1]:eq(E[i].subs(uv,simultaneous=True),0)
e3=s.factor(E[2].subs(uv,simultaneous=True));e4=s.factor(E[3].subs(uv,simultaneous=True))
eq(e3,-3*(4*c-9)*(-96*S*c+12*S+48*T+160*c*c+160*c-47)/64)
ts=(96*S*c-12*S-160*c*c-160*c+47)/48
eq(e4.subs(T,ts),(4*c-25)*(4*c-9)*(4*c-1)*(-8*S+20*c+15)/64)
for center in [Q(1,2),Q(3,2),Q(5,2)]:central_collision(3,3,center)
vals={S:(20*c+15)/8,T:(160*c*c-20*c+49)/96,U:(69-4*c)/8,V:(2*c-23)*(4*c-31)/48}
for e in E:eq(e.subs(vals,simultaneous=True),0)
f15=2*c*(2*c-23)*(4*c-31)*(160*c*c-20*c+49)
eq(96**2*c*vals[T]*vals[V],f15)
curve(15,f15,c,5,2)

# Type 16.
E=[s.expand(central(3,k)+pair_power(k,S,T)-pair_moment(2,k,U,V)-c**k) for k in range(1,5)]
uv={U:(4*S-4*c+31)/8,V:(-16*S*S-32*S*c+296*S+64*T+48*c*c-296*c-73)/128}
for i in [0,1]:eq(E[i].subs(uv,simultaneous=True),0)
e3=s.factor(E[2].subs(uv,simultaneous=True));e4=s.factor(E[3].subs(uv,simultaneous=True))
D=4*S+4*c-51
ts=(64*S**3+64*S*S*c-944*S*S+192*S*c*c-1376*S*c+1996*S-320*c**3+2320*c*c-1996*c-441)/(64*D)
eq(e3.subs(T,ts),0)
eq(e4.subs(T,ts),(4*c-25)*(4*c-1)*(-4*S+4*c+1)*(-4*S+4*c+9)*(-4*S+4*c+49)*(12*S+20*c-195)/(1024*D**2))
eq(e3.subs(S,Q(51,4)-c),-3*(4*c-25)*(4*c-21)*(4*c-1)/64)
for center in [Q(1,2),Q(5,2)]:central_collision(3,1,center)
assert not s.sqrt(Q(21,4)).is_Rational
for j in [1,9,49]:
    sk=c+Q(j,4)
    eq(x*x-sk*x+ts.subs(S,sk),(x-c)*(x-Q(j,4)))
vals={S:(195-20*c)/12,T:(160*c*c-1335*c+2304)/36,U:4*(9-c)/3,V:2*(c-6)*(4*c-9)/3}
for e in E:eq(e.subs(vals,simultaneous=True),0)
f16=6*c*(c-6)*(4*c-9)*(160*c*c-1335*c+2304)
eq(18**2*c*vals[T]*vals[V],f16)
curve(16,f16,c,5,2)
print('PASS: exact moment eliminations, all exceptional branches, discriminants, square products and finite-fiber reconstruction identities for types 5,8,12,13,15,16.')


# Independent remaining-family elimination and genus certificates.
from itertools import combinations
import sympy as s
A,B,C,S,T,X,Z=s.symbols('A B C S T X Z')

def moment(length,k,w):
    offsets=[s.Rational(2*j-length+1,2) for j in range(length)]
    return s.expand(sum(s.binomial(2*k,2*j)*w**(k-j)*sum(x**(2*j) for x in offsets) for j in range(k+1)))
def central(m,k):
    return sum(s.Rational(2*j-2*m+1,2)**(2*k) for j in range(2*m))/2
def pair(k,u,v):
    return [2,u,u*u-2*v,u**3-3*u*v,u**4-4*u*u*v+2*v*v][k]
def pairmoment(length,k,u,v):
    p=s.Poly(moment(length,k,X),X)
    return s.expand(sum(c*pair(mon[0],u,v) for mon,c in p.terms()))
def num(p):return s.factor(s.cancel(p).as_numer_denom()[0])
def proportional(p,q,variables):
    ratio=s.cancel(p/q)
    assert not any(ratio.has(v) for v in variables) and ratio!=0
    return ratio

def eliminate_single_pair(m,left,right,firstvar,secondvar):
    eq=[s.expand(central(m,k)+left(k)-right(k)) for k in range(1,5)]
    r1=s.solve(eq[0],firstvar)[0]
    r2=s.solve(eq[1].subs(firstvar,r1),secondvar)[0]
    subs={firstvar:r1,secondvar:r2}
    rem=[num(e.subs(subs,simultaneous=True)) for e in eq[2:]]
    return r1,r2,rem

# Type 3: A1=A,A3=B; B1=C, two B2 squares have sum S/product T.
a,t,eq=eliminate_single_pair(1,lambda k:A**k+moment(3,k,B),lambda k:C**k+pairmoment(2,k,S,T),A,T)
H3=1024*B*B-1520*B*C+636*B+640*C*C-1860*C+2025
R3=(4*B-25)*(4*B-9)*(4*B-1)*(4*C-1)**3*((B-C-1)**2-4*C)**2*H3
proportional(s.resultant(*eq,S),R3,(B,C))
assert s.Poly(eq[0],S).degree()==3 and not s.Poly(eq[0],S).LC().has(B,C)
# Type 4: A2 pair S,T; B1 pair A,B; B3=C.
a,b,eq=eliminate_single_pair(1,lambda k:pairmoment(2,k,S,T),lambda k:pair(k,A,B)+moment(3,k,C),A,B)
R4=(4*C-9)*(4*C-1)*(-6*S+20*C+25)*(-2*S+4*C+1)*(4*S*S-16*S*C-20*S+16*C*C+24*C+25)
proportional(s.resultant(*eq,T),R4,(S,C))
s0=(20*C+25)/6;t0=(240*C*C+40*C+167)/48
u0=(44*C+91)/12;v0=(256*C*C+97*C+466)/36
assert all(s.cancel(q.subs({S:s0,T:t0},simultaneous=True))==0 for q in eq)
assert s.cancel(a.subs(S,s0)-u0)==0 and s.cancel(b.subs({S:s0,T:t0},simultaneous=True)-v0)==0
P4=3*C*(240*C*C+40*C+167)*(256*C*C+97*C+466)
assert s.cancel(P4/(C*t0*v0))==72**2
# Every extra branch forces A2/B3 overlap. Verify squared A2 centers.
v=s.symbols('v')
for sign in (-1,1):
    sx=2*v*v+s.Rational(5,2)+2*sign*v
    tx=(2*v-sign)**2*(2*v+3*sign)**2/16
    assert all(s.cancel(q.subs({S:sx,T:tx,C:v*v},simultaneous=True))==0 for q in eq)
    proportional(eq[0].subs({S:sx,C:v*v},simultaneous=True),(2*v-3*sign)*(2*v-sign)*(T-tx),(T,v))
    assert s.expand(X*X-sx*X+tx-(X-(v-sign/s.Integer(2))**2)*(X-(v+3*sign/s.Integer(2))**2))==0
sx=2*C+s.Rational(1,2);tx=(4*C-1)**2/16
assert all(s.cancel(q.subs({S:sx,T:tx},simultaneous=True))==0 for q in eq)
proportional(eq[0].subs(S,sx),(4*C+11)*(T-tx),(T,C))
assert s.expand((X*X-sx*X+tx).subs(C,v*v)-(X-(v-s.Rational(1,2))**2)*(X-(v+s.Rational(1,2))**2))==0
# T coefficient can vanish only at already invalid/non-square C on each branch.
proportional(eq[0].subs(S,s0),(4*C+11)*(T-t0),(T,C))

# Types 6/7: A1 pair S,T, A2=A; B1=B,B4=C or B2=B,B3=C.
H6=A**3-5*A*A*C-9*A*A+10*A*C*C+15*A*C+24*A-10*C**3+10*C*C-20*C-16
H7=256*A**3-1280*A*A*C-64*A*A+2560*A*C*C-2560*A*C-2976*A-960*C**3+3920*C*C+780*C+3159
E6=3*A*A-8*A*C-7*A+6*C*C+2*C+4
for typ,right,H in [(6,lambda k:B**k+moment(4,k,C),H6),(7,lambda k:moment(2,k,B)+moment(3,k,C),H7)]:
    ss,tt,eq=eliminate_single_pair(1,lambda k:pair(k,S,T)+moment(2,k,A),right,S,T)
    if typ==6:
        expected=(A-C)*((A-C-1)**2-4*C)*E6*H
        coeff=s.Poly(eq[0],B).all_coeffs()
        proportional(s.resultant(*coeff,A),C*(C-4)**2*(C-1)**3,(C,))
    else:
        J=16*A*A-32*A*C-8*A+16*C*C-8*C+1
        expected=(4*C-9)*(4*C-1)**2*(8*A-12*C-5)*J**2*H
        aa=(12*C+5)/8;bb=(20*C+19)/24
        proportional(eq[0].subs(A,aa),(4*C-9)*(4*C-1)*(24*B-20*C-19),(B,C))
        proportional(eq[1].subs({A:aa,B:bb},simultaneous=True),(4*C-9)*(4*C-1)**2*(4*C+11),(C,))
        # The linear subresultant recovers b off the exceptional fibers.
        br=(10*A-4*C+15)/6
        assert s.expand(H.subs(A,(4*C-15)/10)+27*(4*C-25)*(144*C**2+440*C+1225)/125)==0
        u7,t7=s.symbols('u7 t7')
        a7=45*t7**3-145*t7**2+135*t7-27
        b7=-8*(45*t7**2-80*t7+21)
        c7=144*(5*t7-2)
        f7=27*t7**3-79*t7**2+57*t7-9
        assert s.expand(H.subs({A:1+u7,C:(25+(10-6*t7)*u7)/4})
                        -72*u7*(a7*u7**2+b7*u7+c7))==0
        assert s.expand(b7**2-4*a7*c7-320*f7)==0
        assert s.expand(f7.subs(t7,3*u7**2)
                        -9*(9*u7**3-5*u7**2-3*u7+1)*(9*u7**3+5*u7**2-3*u7-1))==0
        linear=s.Poly(s.subresultants(*eq,B)[-2],B)
        leading=s.rem(linear.nth(1),H,A)
        assert s.rem(linear.nth(0)+leading*br,H,A)==0
        exceptional=(4*C-9)**4*(4*C-1)**6*(4*C+11)*(144*C**2+152*C+529)*(144*C**2+440*C+1225)
        proportional(s.resultant(leading,H,A),exceptional,(C,))
        assert s.expand(ss.subs(B,br)-(16*A+20*C+81)/12)==0
        assert s.cancel(tt.subs(B,br)+(128*A**2-640*A*C+48*A+80*C**2-1380*C-1701)/144)==0
        g=s.groebner(s.Poly(eq[0],B).all_coeffs(),A,C)
        assert all(g.reduce(z)[1]==0 for z in [8*A-12*C-5,16*C*C-40*C+9])
    proportional(s.resultant(*eq,B),expected,(A,C))
    # Smooth cubic over the algebraic closure, including infinity.
    assert list(s.groebner([H,s.diff(H,A),s.diff(H,C)],A,C))==[1]
    HH=s.Poly(H,A,C).homogenize(Z).as_expr()
    infinity=s.Poly(HH.subs({C:1,Z:0}),A)
    assert infinity.degree()==3 and s.discriminant(infinity,A)!=0
    assert HH.subs({A:0,C:1,Z:0})!=0
    assert s.degree(H.subs(A,0),C)==3 and s.discriminant(H.subs(A,0),C)!=0
    print(f'Type {typ}: smooth cubic and six simple branch places verified; double cover genus 4.')

# Type 11: A1=A,A2=B; B1=C, B2 pair S,T.
a,t,eq=eliminate_single_pair(2,lambda k:A**k+moment(2,k,B),lambda k:C**k+pairmoment(2,k,S,T),A,T)
# Shared factor gives A=C, hence cross-color collision.
assert s.cancel(a.subs(S,B+1)-C)==0
eq=[s.cancel(q/(-S+B+1)) for q in eq]
H11=16*B*B+80*B*C-140*B-80*C*C+280*C-281
J=16*B*B-32*B*C-8*B+16*C*C-8*C+1
proportional(s.resultant(*eq,S),(4*C-9)*(4*C-1)*J*H11,(B,C))
assert s.Poly(eq[0],S).degree()==2 and not s.Poly(eq[0],S).LC().has(B,C)

# Type 11: recover the center squares and fix the denominator form's squareclass.
S11=(5*B-4*C+16)/3
L11=16*B*C+4*B-16*C*C+56*C-85
linear_subresultant=s.subresultants(*eq,S)[-2]
proportional(s.rem(linear_subresultant,H11,B),s.rem(L11*(S-S11),H11,B),(B,C,S))
assert s.expand(s.resultant(L11,H11,B)-16*(4*C-49)*(4*C-9)**2*(4*C-1))==0
assert s.expand(H11.subs(C,s.Rational(49,4))-8*(B-9)*(2*B+123))==0
exception={B:9,C:s.Rational(49,4)}
proportional(s.gcd(*(q.subs(exception) for q in eq),S),(S-4)*(S-5),(S,))
assert [a.subs(exception).subs(S,v0) for v0 in (4,5)]==[s.Rational(1,4),s.Rational(9,4)]
a11=a.subs(S,S11)
t11=t.subs(S,S11)
assert s.cancel(a11-(4*B-5*C+26)/3)==0
assert s.rem(t11+(4*C-9)*(5*B-4*C+4)/6,H11,B)==0

w=s.symbols('w')
d11=w*w+5*w-5
fa=25*w*w-53*w+27
fb=41*w*w-70*w+30
fc=49*w*w-37*w-5
fe=(4*w-5)*(47*w-40)
fj=4643*w**4-18384*w**3+27840*w*w-18880*w+4800
param={B:-3*fb/(2*d11),C:fc/(4*d11)}
assert s.cancel(H11.subs(param,simultaneous=True))==0
assert s.cancel((B-9-w*(C-s.Rational(49,4))).subs(param,simultaneous=True))==0
assert s.cancel(a11.subs(param,simultaneous=True)+15*fa/(4*d11))==0
assert s.cancel(t11.subs(param,simultaneous=True)-fe*(5*w-4)**2/(2*d11**2))==0
assert s.cancel((S11*S11-4*t11).subs(param,simultaneous=True)-3*fj/(4*d11**2))==0
assert s.discriminant(fb,w)==-20
assert s.resultant(d11,fa,w)==-36
assert s.resultant(d11,fb,w)==225
assert s.resultant(d11,fc,w)==-1620

modulus=37
square_residues={r*r%modulus for r in range(modulus)}
normalize=[s.Poly(f,w).homogenize(Z).as_expr() for f in (d11,fa,fb,fc,fe)]
projective_residues=[(r,1) for r in range(modulus)]+[(1,0)]
for squareclass,expected in [(-1,[(23,1),(30,1)]),(-5,[(0,1),(4,1)])]:
    factors=[squareclass,-15*squareclass,-6*squareclass,squareclass,2]
    def residue_passes(residue,count):
        subs={w:residue[0],Z:residue[1]}
        return all(int(k*f.subs(subs))%modulus in square_residues
                   for k,f in zip(factors[:count],normalize[:count]))
    assert [r for r in projective_residues if residue_passes(r,4)]==expected
    assert [r for r in projective_residues if residue_passes(r,5)]==(
        [] if squareclass==-1 else expected)
assert all(s.Poly(f,w).degree()==6 and s.degree(s.gcd(f,s.diff(f,w),w),w)==0
           for f in (-2*d11*fa*fb,3*d11*fa*fc,30*d11*fb*fc,-2*fa*fb*fc))

r11=s.symbols('r')
w11=5*(r11-8)/(4*r11-47)
dr=r11*r11-11*r11-1
ar=r11*r11-261*r11-1
br=7*r11*r11-2*r11+18
cr=9*r11*r11-69*r11-49
assert s.cancel(((47*w-40)/(4*w-5)).subs(w,w11)-r11)==0
for original,model,scale in [
        (-d11*fa*fb*fe,r11*dr*ar*br,3375),
        (15*d11*fb*fc*fe,r11*dr*br*cr,50625),
        (-fa*fb*fc*fe,r11*ar*br*cr,3375)]:
    assert s.cancel((4*r11-47)**8*original.subs(w,w11)-scale**2*model)==0
    reduction=s.Poly(model,r11,modulus=37)
    assert reduction.degree()==7 and reduction.gcd(reduction.diff()).degree()==0
    assert reduction.eval(0)==0
    assert int(reduction.eval(8))%37 in square_residues-{0}
assert [(47*w0-40)*pow((4*w0-5)%37,-1,37)%37 for w0 in [0,4]]==[8,0]
assert w11.subs(r11,0)==s.Rational(40,47) and fe.subs(w,w11.subs(r11,0))==0
assert w11.subs(r11,8)==0 and param[C].subs(w,0)==s.Rational(1,4)
print('Type 11: individual square functions and squareclass obstruction at 37 verified.')

# Even plane quartic smoothness: every principal minor is nonzero.
for name,poly,variables in [('3',H3,(B,C)),('11',H11,(B,C))]:
    h=s.Poly(poly,*variables).homogenize(Z).as_expr()
    M=s.hessian(h,(*variables,Z))
    assert all(M.extract(I,I).det()!=0 for n in range(1,4) for I in combinations(range(3),n))
    print(f'Type {name}: all principal minors nonzero; even plane quartic genus 3.')
# Type 18: doubled centered roots have squares 1,9,25,49,81.
f=s.prod(X-(2*j-1)**2 for j in range(1,6));q=s.expand(f-f.subs(X,0)-Z*Z)
P18=s.discriminant(q,X)
assert s.degree(P18,Z)==8 and s.degree(s.gcd(P18,s.diff(P18,Z)),Z)==0
assert s.expand(q-(X**5-165*X**4+8778*X**3-172810*X**2+1057221*X-Z*Z))==0
# Absolute prime cutoff in proofs/type-18-rigid-primes.md.
seed_derivative = s.diff(f, X)
critical_resultant = s.resultant(f - Z, seed_derivative, X)
assert s.degree(critical_resultant, Z) == 4
assert s.discriminant(seed_derivative, X) == 2**42 * 3**5 * 5 * 11**3 * 47 * 59
assert s.discriminant(critical_resultant, Z) == 2**202 * 3**27 * 5**11 * 11**15 * 47**3 * 59**3
assert critical_resultant.subs(Z, f.subs(X, 0)) == 3**10 * 5**6 * 7**2 * 11**5 * 59**2 * 181**2 * 1231 * 2417
assert seed_derivative.subs(X, 0) == 3**2 * 11 * 59 * 181
assert all(s.isprime(p) for p in (2, 3, 5, 7, 11, 47, 59, 181, 1231, 2417))
print('PASS: all branch-exception factors for the type-18 rigid-prime cutoff.')
for poly,var in [(P4,C),(P18,Z)]:
    assert s.degree(s.gcd(poly,s.diff(poly,var)),var)==0
print('Types 4 and 18: squarefree hyperelliptic degrees 5 and 8, genera 2 and 3.')
print('All six remaining-type reductions and exceptional-fiber checks passed exactly.')

assert len(patterns)==18
print("PASS: all eighteen types have explicit finite-fiber maps to smooth curves of genus at least two.")
print("This verifies the algebra in Theorem 14; it neither enumerates rational points nor proves all types empty.")
