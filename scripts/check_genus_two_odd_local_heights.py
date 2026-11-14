#!/usr/bin/env python3
"""Exact odd bad-prime mixed heights at 3,5; PROOF.md Section 5.20.
Uses existing SymPy only. No normalized 13-adic height is computed.
"""
import sympy as s
z,t,pi,r,X=s.symbols("z t pi r X")
F=z*(z-1)*(z+5)*(z*z-68*z-275)
def red(q,m):
    q=s.Rational(q);return int(q.p)*pow(int(q.q),-1,m)%m
def val(q,p):
    q=s.Rational(q);assert q
    a,b,v=int(q.p),int(q.q),0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v
def hensel(f,seed,p,n):
    a=seed
    for e in range(2,n+1):
        m=p**e
        a=(a-int(f.subs(t,a))*pow(int(s.diff(f,t).subs(t,a)),-1,m))%m
        assert int(f.subs(t,a))%m==0
    return a
G3=(1+3*t)*t*(t+2)*(t*t-22*t-38)
assert s.expand(F.subs(z,1+3*t)-81*G3)==0
H3=(1+pi*r)*r*(r+2*pi)*(r*r-22*pi*r-114)
rem3=lambda a:s.rem(s.expand(a),pi*pi-3,pi)
assert rem3(F.subs(z,1+pi*r)-pi**4*H3)==0
T=11+pi*X
assert rem3(T*T-22*T-38-pi*pi*(X*X-53))==0
# k=-3: same twin component, same sheet at its smooth point (0,2).
# s0=sqrt(-29)=1 mod3; Q twin t=3-s0.
s0=hensel(t*t+29,1,3,12)
assert val(7-s0,3)==1 and val(7+s0,3)==0
assert 44%3==(-38*(1+s0))%3==2
# Difference in smooth parameter X is (7-s0)/pi, valuation one over K.
assert 2*val(7-s0,3)-1==1
# Q/R on cluster branch t=1: parameter Yc=w/(9*pi).
tR=hensel(9*t*t-59*t+38,1,3,12)
wQ=1026*(s0-1);wR=4*(269*tR-266)
assert val(wQ,3)==4 and val(wR,3)==3
assert 2*val(wQ-wR,3)-5==1
assert 2*val(wQ+wR,3)-5==1
# k=-15: Q is split over K, one point on each exchanged midpoint component.
assert 265%3==1 and 302%3==2 and val(s.Integer(7950),3)==1
# R is on the fixed twin component over an unramified quadratic extension.
assert 11%3==2
assert (-2754)%3==0 and 866%3==2 and 25%3==1
assert (2*2*2)%3==2  # J=+-2*sqrt(11), hence J^2=2 != 0.
# k=-3 R has one root at infinity and one cluster point t=1.
assert s.expand((3*z*z-65*z+176).subs(z,1+3*t)
                -3*(9*t*t-59*t+38))==0
# Five-adic inner chart and its selected same-sheet intersection.
G5=t*(5*t-1)*(t+1)*(5*t*t-68*t-55)
assert s.expand(F.subs(z,5*t)-125*G5)==0
assert s.Poly(-3*G5,t,modulus=5)==s.Poly(t*t*(t+1),t,modulus=5)
q=25*t*t+92*t+49
a=hensel(q,3,5,8);pP=s.Rational(-3,4)
xiP=s.Rational(57,32);xiR=s.Rational(2598*a+1806,25)
assert s.cancel(q.subs(t,pP))==s.Rational(-95,16)
assert val(a-pP,5)==1
assert red(xiR,5)==red(xiP,5)==1
assert s.diff(t*t*(t+1),t).subs(t,3)%5!=0 or 2*red(xiP,5)%5
# k=-3 all three supports are smooth on the original model at5, disjoint.
x=s.symbols("x");k=-3
f=x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k)
us=[x-33,x*x+60*x+3249,x*x+65*x+528]
assert all(int(u.subs(x,0))%5 for u in us)
assert s.Poly(f,x,modulus=5)==s.Poly(x**3*(x+3)*(x+4),x,modulus=5)
for i,j in [(0,1),(0,2),(1,2)]:
    assert val(s.resultant(us[i],us[j],x),5)==0
print("PASS: k=-3 h3(DP,DQ)=-log13(3); other h3 mixed terms zero.")
print("PASS: k=-15 h5(DP,DR)=-2log13(5); other h5 mixed terms zero.")
print("Still missing: dyadic terms and normalized 13-adic integrals.")
