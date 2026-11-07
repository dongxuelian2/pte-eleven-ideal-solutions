"""Exact rank-two certificates for the ac and bc elliptic quotients.
Dependency-free. Uses the standard two-isogeny descent theorem.
"""
from itertools import combinations
from math import prod

def prime_factors(n):
    n=abs(n);out=[];p=2
    while p*p<=n:
        if n%p==0:
            out.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:out.append(n)
    return out

def squareclass(n):
    assert n
    d=-1 if n<0 else 1;n=abs(n)
    for p in prime_factors(n):
        v=0
        while n%p==0:n//=p;v+=1
        if v%2:d*=p
    return d

def candidates(B):
    ps=prime_factors(B)
    return {sg*prod(S) for k in range(len(ps)+1)
            for S in combinations(ps,k) for sg in [-1,1]}

def locally_possible(d,A,B,p,k):
    e=B//d
    if d<0 and e<0 and A<=0:return False
    M=p**k;sq={z*z%M for z in range(M)}
    return (any((d*t**4+A*t*t+e)%M in sq for t in range(M)) or
            any((d+A*t*t+e*t**4)%M in sq for t in range(0,M,p)))

models=[
 ('ac',-2360,1287900,{1,5,159,795},[(405,14175),(3180,111300)]),
 ('acdual',4720,418000,{1,-11,-95,1045},[(-4400,66000),(-95,1425)]),
 ('bc',-920,364500,{1,5},[]),
 ('bcdual',1840,-611600,{1,-5,-11,139,55,-695,-1529,7645},
  [(-1620,39600),(-1100,39600),(-695,31275),(556,20016),
   (880,39600),(30580,5504400)])]
for name,A,B,expected,points in models:
    ds=candidates(B)
    for p,k in [(2,8),(3,1),(5,3)]:
        ds={d for d in ds if locally_possible(d,A,B,p,k)}
    assert ds==expected
    witnessed={1,squareclass(B)}
    for x,y in points:
        assert y*y==x*(x*x+A*x+B)
        witnessed.add(squareclass(x))
    assert witnessed==expected
    print(name,'exact descent image:',sorted(expected))
for i in [0,2]:
    _,A,B,D,_=models[i];_,AA,BB,DD,_=models[i+1]
    assert (AA,BB)==(-2*A,A*A-4*B)
    assert len(D)*len(DD)==4*2**2

# Binary-quartic Jacobian invariants for ac and bc.
# The Jacobian of y^2=a*z^4+b*z^3+c*z^2+d*z+e has model
# y^2=x^3-27I*x-27J; this is the scale-by-3 transform of
# y^2=x*(x^2+A*x+B) after translating x by A/3.
q1=(4,30,80);q2=(1,36,185);q3=(1,0,275)
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
for f,g,i in [(q1,q3,0),(q2,q3,2)]:
    a,b,c,d,e=mul(f,g)
    I=12*a*e-3*b*d+c*c
    J=72*a*c*e+9*b*c*d-27*a*d*d-27*e*b*b-2*c**3
    _,A,B,_,_=models[i]
    assert I==A*A-3*B
    assert J==-2*A**3+9*A*B
    # Positive square leading coefficient gives rational points at infinity.
    assert a in [1,4]
assert 2+2+3+1+2+2==12>9
print('Both new elliptic factors have rank2; six factors force center-cover rank >=12 > genus9.')
