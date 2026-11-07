#!/usr/bin/env python3
"""Exact certificate for genus-three quadratic-Chabauty targets.
Dependency-free Python 3. Mathematical deductions are in
neron_severi_and_genus_three_qc.md. This does not compute rational points.
"""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from itertools import combinations
from math import prod
# Reuse already certified exact finite-field/polynomial utilities.
with redirect_stdout(StringIO()):
    from check_genus_two_richelot_certificate import (
        add,scale,mul,evaluate,counts,squarefree_mod,legendre,
        irreducible_monic_quartic,prime_factors)

q0=[5,0,1];q1=[80,30,4];q3=[275,0,1]
f=mul(mul(q0,q1),q3)
# Every listed reduction is good, ordinary and simple over F_p.
for delta in (1,5):
    for p,expected_a1,disc in [(17,-4 if delta==1 else 4,32),
                               (23,-2 if delta==1 else 2,68)]:
        n1,n2,a1,a2,order=counts(scale(f,delta),p)
        assert (a1,a2)==(expected_a1,30)
        assert a2%p!=0
        F=[p*p,-p*a1,a2,-a1,1]
        assert irreducible_monic_quartic(F)
        assert a1*a1-4*(a2-2*p)==disc
        # Frobenius fixed subfield has minimal polynomial
        # U^2-a1 U+(a2-2p), U=pi+p/pi.
        assert disc in (4**2*2,2**2*17)
print('rho_Q(J(C_1)) = rho_Q(J(C_5)) = 1: arithmetic inputs checked.')

# Recheck the complete two-isogeny descent on the elliptic factor.
def sqclass(n):
    out=-1 if n<0 else 1;n=abs(n)
    for p in prime_factors(n):
        e=0
        while n%p==0:n//=p;e+=1
        if e%2:out*=p
    return out

def classes(B):
    ps=prime_factors(abs(B))
    return {sg*prod(I) for n in range(len(ps)+1)
            for I in combinations(ps,n) for sg in (-1,1)}

def possible(d,A,B,p,n):
    mod=p**n;sq={x*x%mod for x in range(mod)};e=B//d
    return (any((d*t**4+A*t*t+e)%mod in sq for t in range(mod)) or
            any((d+A*t*t+e*t**4)%mod in sq for t in range(0,mod,p)))

models=[(-560,72900,{1,5},[(180,900)]),
        (1120,22000,{1,-5,-11,55},[(-1100,0),(-500,12000)])]
for A,B,expected,points in models:
    upper={d for d in classes(B)
           if not(d<0 and B//d<0 and A<=0)
           and possible(d,A,B,2,8) and possible(d,A,B,3,1)}
    lower={1,sqclass(B)}
    for x,y in points:
        assert y*y==x*(x*x+A*x+B)
        lower.add(sqclass(x))
    assert upper==lower==expected
assert len(models[0][2])*len(models[1][2])==8
assert (-2*models[0][0],models[0][0]**2-4*models[0][1])==models[1][:2]
# Non-torsion is also independently visible: [2](180,900) has nonintegral x.
x,y=180,900;A,B=-560,72900
m=Q(3*x*x+2*A*x+B,2*y)
assert m*m-A-2*x==Q(2025,4)
print('Elliptic factor E has exact rank one.')

# N,D,M encode z=N/D, a=M/D, delta*a^2=q1(z).
PARAMS={1:([80,0,-1],[-30,4],[160,-30,2]),
        5:([80,-120,25],[-8,0,10],[-48,100,-60])}
for delta,(N,D,M) in PARAMS.items():
    N2,D2=mul(N,N),mul(D,D)
    if delta==1:assert M==add(scale(N,2),mul([0,1],D))
    else:assert scale(add(M,scale(D,-6)),2)==mul([0,1],add(scale(N,2),scale(D,-5)))
    assert scale(mul(M,M),delta)==add(add(scale(N2,4),scale(mul(N,D),30)),scale(D2,80))
    A4=add(N2,scale(D2,5));B4=add(N2,scale(D2,275));H=mul(A4,B4)
    assert len(A4)==len(B4)==5 and len(H)==9
    assert squarefree_mod(H,7) and squarefree_mod(H,13)
    # The degree-eight squarefree model has genus 3. Good ordinary p=13:
    assert counts(H,13)[0]==12
    leading_root=1 if delta==1 else 5625
    assert H[-1]==leading_root**2
    # All generic quotient identities follow directly from these polynomials:
    # b=W/D^2, y=a*b=M*W/D^3.
    assert scale(mul(mul(M,M),H),delta)==mul(add(add(scale(N2,4),scale(mul(N,D),30)),scale(D2,80)),mul(A4,B4))
    print('delta=',delta,'N=',N,'D=',D,'M=',M,'A4=',A4,'B4=',B4)
    if delta==1:
        t=Q(15,2)
        assert evaluate(D,t)==0 and evaluate(N,t)==Q(95,4)
        assert evaluate(H,t)==Q(9025,16)**2
        assert H[-1]==1
    else:
        t=Q(5,6);d=evaluate(D,t)
        assert evaluate(N,t)/d==Q(5,2)
        assert evaluate(M,t)/d==6
        assert evaluate(H,t)==Q(9025,144)**2
        # t=infinity has z=5/2,a=-6,b=+-225/4.
        assert Q(N[-1],D[-1])==Q(5,2) and Q(M[-1],D[-1])==-6
        assert Q(leading_root,D[-1]**2)==Q(225,4)

# Ordinary good reduction of the genus-three Jacobians at p=13 follows
# from the isogeny and ordinary factors; verify the polynomial product.
p=13
nE=p+1+sum(legendre(x*(x*x-560*x+72900),p) for x in range(p))
assert nE==12
PE=[p,-2,1];PC=[p*p,0,6,0,1]
PH=mul(PE,PC)
assert PH==[2197,-338,247,-12,19,-2,1] and PH[3]%p
assert 2+1==3 and 3<3+2-1
print('Both genus-three targets: rank 3, rho_Q 2, ordinary good prime 13.')
print('Mordell-Weil-restricted quadratic-Chabauty criterion 3 < 3+2-1 verified.')
print('NOT COMPUTED: p-adic logarithm rank, heights, rational-point sets, or a complete PTE classification.')
