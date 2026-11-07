"""Exact 2-isogeny descent images for the sixth-mode elliptic factors.
Dependency-free. The descent theorem and isogenies are proved/cited in PROOF.
No rational-point search or rank oracle is used.
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
    sign=-1 if n<0 else 1;n=abs(n);d=sign
    for p in prime_factors(n):
        exponent=0
        while n%p==0:n//=p;exponent+=1
        if exponent%2:d*=p
    return d

def candidates(B):
    ps=prime_factors(B)
    return {sg*prod(I) for k in range(len(ps)+1)
            for I in combinations(ps,k) for sg in (-1,1)}

M=256
squares={z*z%M for z in range(M)}
def locally_possible(d,A,B):
    assert B%d==0;e=B//d
    # Both endpoint coefficients negative and A<=0: no real point.
    if d<0 and e<0 and A<=0:return False
    # A primitive (U,V) has V odd, or U odd and V even.
    # Scale the odd coordinate to 1 modulo 256; Z scales by its square.
    return (any((d*t**4+A*t*t+e)%M in squares for t in range(M)) or
            any((d+A*t*t+e*t**4)%M in squares for t in range(0,M,2)))

models=[
 ('v',-95,2430,{1,5,6,30},[(24,132),(45,90)]),
 ('vdual',190,-695,{1,-5,139,-695},[(-125,1050),(139,2502)]),
 ('u',-200,8100,{1,5,6,30},[(6,204),(20,300),(30,300)]),
 ('udual',400,7600,{1,-5,19,-95},[(-380,0),(-320,2400)]),
 ('r',-560,131220,{1,5,6,30},[(6,876),(30,1860)]),
 ('rdual',1120,-211280,{1,-5,19,139,-95,-695,2641,-13205},
  [(-1280,2880),(-695,18765),(-380,13680),(304,8208),
   (556,20016),(10564,1140912)])]
counts=[]
for name,A,B,expected,points in models:
    assert B and A*A-4*B
    surviving={d for d in candidates(B) if locally_possible(d,A,B)}
    assert surviving==expected
    witnessed={1,squareclass(B)}  # infinity and (0,0)
    for x,y in points:
        assert y*y==x*(x*x+A*x+B)
        witnessed.add(squareclass(x))
    assert witnessed==surviving
    counts.append(len(surviving))
    print(name,'exact descent image:',sorted(surviving))
for i,rank in zip((0,2,4),(2,2,3)):
    _,A,B,_,_=models[i];_,AA,BB,_,_=models[i+1]
    assert (AA,BB)==(-2*A,A*A-4*B)
    assert counts[i]*counts[i+1]==4*2**rank
# Even-quartic -> intermediate cubic -> 2-isogenous cubic normalization.
for (qA,qB,qC),scale,i in [((-80,760,139),4,0),((-5,100,-95),1,2),((139,1120,-1520),2,4)]:
    _,A,B,_,_=models[i]
    assert -2*qB==A*scale**2
    assert qB*qB-4*qA*qC==B*scale**4
print('Certified elliptic ranks 2, 2, 3; their product has rank 7.')
