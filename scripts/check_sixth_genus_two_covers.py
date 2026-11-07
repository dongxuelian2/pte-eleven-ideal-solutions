"""Exact finite covering identities; NOT a rational-point classification.
Dependency-free. Genus and global squareclass deductions are in the proof.
"""
from fractions import Fraction as F
from itertools import combinations
from math import prod

# Polynomial coefficients are listed in descending degree.
q=[(1,0,5),(4,30,80),(1,36,185),(1,0,275)]
def mul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
def at(f,x):
    r=0
    for c in f:r=r*x+c
    return r
def quad_resultant(f,g):
    a,b,c=f;d,e,h=g
    return (a*h-c*d)**2-(a*e-b*d)*(b*h-c*e)
assert [b*b-4*a*c for a,b,c in q]==[-20,-380,556,-1100]
assert [quad_resultant(q[i],q[j]) for i,j in combinations(range(4),2)]==[
    8100,38880,72900,131220,1287900,364500]
expected=[
 [4,174,3000,56280,537300,2318250,4070000],
 [1,36,465,10080,53175,49500,254375],
 [4,30,1200,8400,27900,41250,110000],
 [4,174,1920,9300,24300,42150,74000]]
for j in range(4):
    f=[1]
    for i in range(4):
        if i!=j:f=mul(f,q[i])
    assert f==expected[j]

# Check coordinate identities as polynomial identities by degree+1 values.
# Multiplying by 3 removes denominators after t=(z+5)/3.
for z in [F(-2),F(0),F(1)]:
    t=(z+5)/3
    d=3*t*t-10*t+10;p=6*t*t-5*t+5
    bnum=3*t*t+26*t+10;cnum=3*t*t-10*t+100
    assert (3*d,6*p,3*bnum,3*cnum)==tuple(at(f,z) for f in q)

# Only squarefree divisors of 270 can occur by the primitive valuation argument.
candidates={1,2,3,5,6,10,15,30}
# Primitive pairs modulo 8 capture the required parity exclusions.
def homogeneous(f,a,b):
    x,y,z=f;return x*a*a+y*a*b+z*b*b
survivors=set()
for delta in candidates:
    residues={delta*c*c%8 for c in range(8)}
    if any(all(homogeneous(f,a,b)%8 in residues for f in q)
           for a in range(8) for b in range(8) if a%2 or b%2):
        survivors.add(delta)
assert survivors=={1,5}
for delta in [1,5]:
    residues={delta*c*c%8 for c in range(8)}
    for a in range(8):
        for b in range(8):
            if not(a%2 or b%2):continue
            if all(homogeneous(f,a,b)%8 in residues for f in q):
                assert a%2==1 and b%2==0

# Formal identity: X=2s+2z^2+280, with z^2=delta*C0^2/B^2-5,
# s=delta*C0*C3/B^2, delta*(C3^2-C0^2)=270B^2.
# The difference X-delta*(C0+C3)^2/B^2 is exactly
# 270-delta*(C3^2-C0^2)/B^2, which vanishes by that relation.
assert 280-2*5==270
# Check the known finite control, including all four quotient ordinates.
z=F(5,2);delta=5;cs=[F(3,2),F(6),F(15,2),F(15,2)]
assert [at(f,z) for f in q]==[delta*c*c for c in cs]
ys=[delta*prod(cs[i] for i in range(4) if i!=j) for j in range(4)]
assert ys==[F(3375,2),F(3375,8),F(675,2),F(675,2)]
for f,y in zip(expected,ys):assert delta*y*y==at(f,z)
for sign in [-1,1]:
    s=sign*delta*cs[0]*cs[3]
    X=2*s+2*z*z+280
    assert X==delta*(cs[0]+sign*cs[3])**2
    assert X in {180,405}
# The leading coefficients are squares, giving two rational infinities
# on each Delta=1 genus-two curve; Delta=5 is nonsquare there.
assert [f[0] for f in expected]==[4,1,4,4]
assert 16*(-2)+8*8==2*17-2
assert (2*17-2)==2*(2*9-2)
print('Two squareclasses and eight genus-two quotient equations verified.')
print('No Jacobian ranks or complete rational-point classification are asserted.')
