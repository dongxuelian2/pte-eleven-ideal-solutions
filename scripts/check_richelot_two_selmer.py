#!/usr/bin/env python3
"""Deterministic exact certificate for the two Richelot Jacobian 2-Selmer groups.
Also certifies both Mordell-Weil ranks equal two.
No CAS, floating point, searches, or nonstandard packages are required.
See PROOF.md, Section 5.13, for the mathematical deductions.
"""
from fractions import Fraction as F
from itertools import product
D=159
PRIMES=[2,3,5,11,19,53]
QB=[-1]+PRIMES
KB=[(-1,0),(2,0),(3,0),(5,0),(11,0),(19,0),(0,1),(13,1),(12,1),(25,2),(11,1)]

def vp(x,p):
 x=F(x)
 if not x:return 10**9
 a,b=x.numerator,x.denominator; v=0
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v

def pw(p,n):return F(p**n) if n>=0 else F(1,p**(-n))

def mod(x,m):
 x=F(x);return x.numerator*pow(x.denominator,-1,m)%m

def qclass(x,p):
 if not x:raise ValueError('zero qclass')
 if p=='R':return int(x<0)
 v=vp(x,p);u=x/pw(p,v)
 if p==2:
  r=mod(u,8);return (v%2)|({1:0,3:1,5:2,7:3}[r]<<1)
 return v%2 | (int(pow(mod(u,p),(p-1)//2,p)==p-1)<<1)

def tonelli(n,p):
 n%=p
 if not n:return 0
 if pow(n,(p-1)//2,p)!=1:raise ValueError('nonsquare')
 if p%4==3:return pow(n,(p+1)//4,p)
 q=p-1;s=0
 while q%2==0:q//=2;s+=1
 z=2
 while pow(z,(p-1)//2,p)!=p-1:z+=1
 c=pow(z,q,p);x=pow(n,(q+1)//2,p);t=pow(n,q,p);m=s
 while t!=1:
  i=1;t2=t*t%p
  while t2!=1:t2=t2*t2%p;i+=1
  b=pow(c,1<<(m-i-1),p);x=x*b%p;t=t*b*b%p;c=b*b%p;m=i
 return x

def sqrtunit(x,p,N):
 m=p**N;n=mod(x,m)
 if p==2:
  assert n%8==1
  r=1
  for j in range(3,N):
   if (r*r-n)%(1<<(j+1)):r+=1<<(j-1)
  return r
 r=tonelli(n,p);mm=p
 for j in range(1,N):
  r+=(((n-r*r)//mm)*pow(2*r,-1,p)%p)*mm;mm*=p
 return r

def kmul(x,y):return (x[0]*y[0]+D*x[1]*y[1],x[0]*y[1]+x[1]*y[0])

def kdiv(x,y):
 a,b=y;n=a*a-D*b*b
 return kmul(x,(F(a,n),F(-b,n)))

def kpow(x,n):
 if n<0:return kpow(kdiv((1,0),x),-n)
 r=(F(1),F(0))
 while n:
  if n&1:r=kmul(r,x)
  x=kmul(x,x);n//=2
 return r

def signk(a,b):
 if not b:return int(a<0)
 if not a:return int(b<0)
 if a*b>0:return int(a<0)
 return int(a<0) if a*a>D*b*b else int(b<0)

def unit_mul(x,y):return tuple(t%8 for t in kmul(x,y))

USQUARE={unit_mul(x,x) for x in product(range(8),repeat=2) if sum(x)%2}

UCOS={x:0 for x in USQUARE};UB=[]

for x in product(range(8),repeat=2):
 if sum(x)%2 and x not in UCOS:
  b=1<<len(UB);UB.append(x)
  UCOS.update({unit_mul(x,y):v|b for y,v in list(UCOS.items())})

assert len(UB)==3 and len(UCOS)==32

def kclass(a,b,p):
 a,b=F(a),F(b)
 if not a and not b:raise ValueError('zero kclass')
 if p=='R':return signk(a,b)|(signk(a,-b)<<1)
 if p==2:
  v=vp(a*a-D*b*b,2)
  u=kdiv((a,b),kpow((1,1),v))
  return v%2 | (UCOS[(mod(u[0],8),mod(u[1],8))]<<1)
 if p in (3,53):
  v=min(2*vp(a,p),2*vp(b,p)+1);t=v//2
  u=(a if v%2==0 else b)/(D**t if t>=0 else F(1,D**(-t)))
  return v%2 | (int(pow(mod(u,p),(p-1)//2,p)==p-1)<<1)
 N=12
 while True:
  r=sqrtunit(F(D),p,N)
  vals=[a+b*r,a-b*r]
  if all(x and vp(x,p)<N+vp(b,p) for x in vals):
   return qclass(vals[0],p)|(qclass(vals[1],p)<<2)
  N*=2

def lc(coords,p):
 n=1 if p=='R' else (3 if p==2 else 2)
 c=0
 for i in range(3):c|=qclass(coords[i],p)<<(n*i)
 return c | (kclass(*coords[3],p)<<(3*n))

def span_add(basis,x):
 while x:
  i=x.bit_length()-1
  if i not in basis:basis[i]=x;return True
  x^=basis[i]
 return False

def residue(basis,x):
 for i in sorted(basis,reverse=True):
  if (x>>i)&1:x^=basis[i]
 return x

def fval(x,k):return x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k)

def pointcoords(x,k):return (x,x-k,x+5*k,(x-34*k,-3*k))

def torscoords(e,k):
 c=list(pointcoords(e,k));idx=[0,k,-5*k].index(e)
 d=F(1)
 for r in [0,k,-5*k]:
  if e!=r:d*=e-r
 d*=e*e-68*k*e-275*k*k
 c[idx]=d
 return c

def global_columns():
 cs=[]
 for i in range(3):
  for q in QB:
   c=[1,1,1,(1,0)];c[i]=q;cs.append(c)
 for z in KB:cs.append([1,1,1,z])
 return cs

def normclass(c):
 z=c[0]*c[1]*c[2]*(c[3][0]**2-D*c[3][1]**2)
 v=int(z<0)
 for i,p in enumerate(PRIMES,1):v|=(vp(z,p)%2)<<i
 return v

def transpose_cols(cs):
 for b in range(max(cs).bit_length()):yield sum(((c>>b)&1)<<i for i,c in enumerate(cs))

WITNESSES={-15: {'R': [['tors', 0, 20], ['tors', -15, 23]], '2': [['tors', 0, 6454], ['tors', -15, 6896], ['tors', 75, 4778], ['point', '-100', 5142], ['point', '-81', 3646]], '3': [['tors', 0, 60], ['tors', -15, 45], ['tors', 75, 145]], '5': [['tors', 0, 365], ['tors', -15, 891], ['tors', 75, 622], ['point', '-100', 108]], '11': [['tors', 0, 417], ['tors', -15, 514], ['tors', 75, 136], ['point', '-95', 40]], '19': [['tors', 0, 520], ['tors', -15, 900], ['tors', 75, 554], ['point', '-99', 130]], '53': [['tors', 0, 162], ['tors', -15, 128], ['point', '-94', 138]]}, -3: {'R': [['tors', 0, 20], ['tors', -3, 23]], '2': [['tors', 0, 5142], ['tors', -3, 6132], ['tors', 15, 7822], ['point', '-100', 6454], ['point', '-95', 1600]], '3': [['tors', 0, 148], ['tors', -3, 143], ['point', '-75', 99]], '5': [['tors', 0, 57], ['tors', -3, 554], ['tors', 15, 811], ['point', '-97', 674]], '11': [['tors', 0, 417], ['tors', -3, 514], ['tors', 15, 136], ['point', '-100', 554]], '19': [['tors', 0, 520], ['tors', -3, 900], ['tors', 15, 554], ['point', '-93', 130]], '53': [['tors', 0, 10], ['tors', -3, 34], ['point', '-100', 128]]}}

EXPECTED={-15:[25,23,18,15,11,8,5,5],-3:[25,23,18,15,11,8,6,5]}

def rank(vectors):
    b={}
    for v in vectors:span_add(b,v)
    return len(b)

def nullbasis(rows,n):
    b={}
    for r in rows:span_add(b,r)
    ans=[]
    for f in range(n):
        if f in b:continue
        v=1<<f
        for i in sorted(b):
            if bin(b[i]&v).count('1')%2:v|=1<<i
        assert all(bin(r&v).count('1')%2==0 for r in rows)
        ans.append(v)
    return ans

def poly_mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c

def poly_rem(f,u):
    f=list(f)
    assert u[-1]==1
    while len(f)>=len(u):
        t=f[-1];shift=len(f)-len(u)
        for i,c in enumerate(u):f[shift+i]-=t*c
        assert f[-1]==0
        f.pop()
    return f

def quadcoords(A,B,k):
    # Evaluate x²+A*x+B at 0,k,−5k,k*(34+3s).
    return (B,k*k+A*k+B,25*k*k-5*A*k+B,
            (2587*k*k+34*A*k+B,204*k*k+3*A*k))

def resultant(a,b):
    # Exact Sylvester determinant by fraction-preserving rational elimination.
    m,n=len(a)-1,len(b)-1
    rows=[]
    for shift in range(n):rows.append([F(0)]*shift+list(map(F,reversed(a)))+[F(0)]*(n-1-shift))
    for shift in range(m):rows.append([F(0)]*shift+list(map(F,reversed(b)))+[F(0)]*(m-1-shift))
    det=F(1)
    for i in range(m+n):
        j=next((j for j in range(i,m+n) if rows[j][i]),None)
        if j is None:return F(0)
        if j!=i:rows[i],rows[j]=rows[j],rows[i];det=-det
        piv=rows[i][i];det*=piv
        for j in range(i+1,m+n):
            t=rows[j][i]/piv
            rows[j]=[x-t*y for x,y in zip(rows[j],rows[i])]
    return det

def verify_global_setup():
    # Discriminant and all possible bad primes.
    h=[1]
    for f in ([0,1],[-1,1],[5,1],[-275,-68,1]):h=poly_mul(h,f)
    disc=resultant(h,[(i+1)*h[i+1] for i in range(5)])
    assert disc==2**8*3**13*5**8*11**2*19**2*53
    assert D%4==3 and 4*D==636 and 636<26**2
    assert pow(D%5,2,5)==1 and pow(D%7,3,7)==6
    assert pow(D%11,5,11)==1 and pow(D%19,9,19)==1
    # Minkowski bound sqrt(636)/2<13. Norm identities imply the ideal
    # class group is generated by P2, hence becomes trivial after S-inversion.
    for a,b,n in [(13,1,10),(12,1,-15),(25,2,-11),(11,1,-38)]:
        assert a*a-D*b*b==n
    for a,b in KB:
        n=abs(a*a-D*b*b)
        for p in PRIMES:
            while n%p==0:n//=p
        assert n==1
    places=['R']+PRIMES
    combined=[sum(kclass(a,b,p)<<(8*i) for i,p in enumerate(places)) for a,b in KB]
    assert rank(combined)==11 # Independent in K*/K*²; dimension theorem gives basis.
    # Exhaustive unit-mod-8 computation underlying the dyadic K square-class map.
    assert USQUARE=={(1,0),(5,4),(7,0),(3,4)}
    assert len(UB)==3 and len(UCOS)==32
    for x in UCOS:
        for y in UCOS:assert UCOS[unit_mul(x,y)]==UCOS[x]^UCOS[y]
    # Local maps are homomorphisms on the global basis (additional error check).
    for p in places:
        for a in KB:
            for b in KB:assert kclass(*kmul(a,b),p)==kclass(*a,p)^kclass(*b,p)
    return disc

def check_twist(k):
    cs=global_columns();eq={}
    assert len(cs)==32
    for r in transpose_cols([normclass(c) for c in cs]):span_add(eq,r)
    dims=[len(cs)-len(eq)]
    locals={}
    for p in ['R']+PRIMES:
        b={}
        for typ,data,expected in WITNESSES[k][str(p)]:
            if typ=='tors':
                assert data in [0,k,-5*k]
                coords=torscoords(data,k)
            else:
                assert typ=='point';x=F(data);y2=fval(x,k)
                assert y2 and (y2>0 if p=='R' else qclass(y2,p)==0)
                coords=pointcoords(x,k)
            val=lc(coords,p)
            assert val==expected
            assert span_add(b,val)
        target=2 if p=='R' else (5 if p==2 else (3 if p in [3,53] else 4))
        assert len(b)==target # Therefore the exhibited Kummer image is complete.
        locals[p]=b
        for r in transpose_cols([residue(b,lc(c,p)) for c in cs]):span_add(eq,r)
        dims.append(len(cs)-len(eq))
    assert dims==EXPECTED[k]
    assert len(eq)==27
    sb=nullbasis(list(eq.values()),32)
    assert len(sb)==5 and rank(sb)==5
    # The known non-Weierstrass point also lies in every certified local image.
    z,Y=(F(-15,4),F(1425,32)) if k==-15 else (F(-11),F(1188))
    assert Y*Y==k*z*(z-1)*(z+5)*(z*z-68*z-275)
    x=k*z
    for p,b in locals.items():assert residue(b,lc(pointcoords(x,k),p))==0
    # Its Kummer class is independent of the three rational branch classes.
    places=['R']+PRIMES
    def combined(coords):return sum(lc(coords,p)<<(20*i) for i,p in enumerate(places))
    low=[combined(torscoords(e,k)) for e in [0,k,-5*k]]+[combined(pointcoords(x,k))]
    assert rank(low)==4
    # A fixed rational degree-two divisor gives a second independent free class.
    A,B,m,n=(-68,361,906,-6954) if k==-15 else (-20,361,342,-4446)
    original=[k]
    for fac in ([0,1],[-1,1],[5,1],[-275,-68,1]):original=poly_mul(original,fac)
    v2=poly_mul([n,m],[n,m])
    for i,c in enumerate(v2):original[i]-=c
    assert poly_rem(original,[B,A,1])==[0,0]
    # Under x=k*z, u becomes x²+k*A*x+k²*B.
    coords=quadcoords(k*A,k*k*B,k)
    for p,b in locals.items():assert residue(b,lc(coords,p))==0
    low.append(combined(coords))
    assert rank(low)==5
    print('k =',k,'; successive global dimensions =',dims)
    print('2-Selmer basis bit masks:',sb)
    print('Degree-two divisor: u = z² + (%s)z + %s; v = %sz + (%s).'%(A,B,m,n))
    print('dim Sel_2 = 5; rational 2-torsion dimension = 3; exact rank = 2.')
    return dims

if __name__=='__main__':
    verify_global_setup()
    for k in (-15,-3):check_twist(k)
    print('All deterministic 2-descent certificate checks passed.')
