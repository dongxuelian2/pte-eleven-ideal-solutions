#!/usr/bin/env python3
"""Exact Richelot de Rham and splitting transport; PROOF Section 5.22.
Requires existing SymPy. Computes no normalized Coleman integrals.
"""
import sympy as s
x,z,t,k,delta=s.symbols("x z t k delta")
A=x*x+5;B=4*x*x+30*x+80;C=x*x+275;f=A*B*C
Q=z*z-68*z-275;F=z*(z-1)*(z+5)*Q
E=A*Q+18*B*z
Ls=[-30*Q,-540*z,30*(z-1)*(z+5)]
assert s.expand(A*Ls[0]+B*Ls[1]+C*Ls[2]-8100*(x-z)**2)==0
# Squared ordinate identity yY=k*A*Q*(x-z), with delta*k=-15.
test=s.together((f/delta)*(k*F)-(k*A*Q*(x-z))**2).subs(delta,-15/k)
assert s.rem(s.numer(s.cancel(test)),E,z)==0
# Direct trace in the quadratic correspondence algebra.
ep=s.Poly(E,z,domain=s.QQ.frac_field(x))
reduce=lambda h:s.rem(h,E,z)
inverse=lambda h:s.invert(h,E,z)
# y/Y=f/(delta*k*A*Q*(x-z))=-B*C/(15*Q*(x-z)).
# z'=-E_x/E_z, so trace coefficient below multiplies dx/y.
g=reduce(B*C*s.diff(E,x)*inverse(15*Q*(x-z)*s.diff(E,z)))
def trace(h):
    h=s.Poly(reduce(h),z)
    return s.factor(2*h.nth(0)-h.nth(1)*ep.nth(1)/ep.nth(2))
direct=[trace(z**j*g) for j in range(4)]
assert direct[:2]==[2,2*x]
basis=[1,x,x**3+s.Rational(15,4)*x*x,x**4+s.Rational(4125,32)*x*x]
T=s.Matrix([[2,0,0,0],[0,2,0,0],
 [0,-s.Rational(32,3),-s.Rational(8,15),0],
 [-s.Rational(2200,9),-s.Rational(2464,9),-s.Rational(1114,45),-s.Rational(16,45)]])
H=-s.Rational(176,9)-s.Rational(214,9)*x+s.Rational(256,45)*x*x+s.Rational(2,45)*x**3
R=[0,0,s.Rational(2,15)/A,H/A**2]
for j in range(4):
    exact=s.diff(R[j],x)*f+R[j]*s.diff(f,x)/2
    assert s.cancel(direct[j]-exact-sum(T[j,i]*basis[i] for i in range(4)))==0
assert T.det()==s.Rational(512,675)
CC=s.Matrix([[0,0,0,delta/4],[0,0,delta/2,-45*delta/16],
 [0,-delta/2,0,525*delta/2],[-delta/4,45*delta/16,-525*delta/2,0]])
# Independently compute target cups as Laurent residues at its infinity.
series=s.series(s.expand(t**10*F.subs(z,t**-2))**s.Rational(-1,2),t,0,10).removeO()
forms=[s.expand(-2*t**(2-2*j)*series) for j in range(4)]
primitives=[s.integrate(w,t) for w in forms]
CD=s.Matrix(4,4,lambda i,j:s.expand(primitives[i]*forms[j]).coeff(t,-1)/k)
assert (T*CC*T.T-2*CD).subs(delta,-15/k).applyfunc(s.simplify)==s.zeros(4)
W=s.Matrix([[0,s.Rational(16,3),1,0],
 [s.Rational(1100,9),s.Rational(392,9),0,1]])
assert W*CD*W.T==s.zeros(2)
image=W*T
assert image[:,0]==s.zeros(2,1)
assert image[:,1]==525*image[:,3]
for kk in [-15,-3]:
    scale=s.diag(*[s.Rational(kk)**(j-1)/2 for j in range(4)])
    M=scale*T
    assert M.det()==s.Rational(32*kk*kk,675)
    for entry in list(M)+list(M.inv()):
        assert s.Rational(entry).q%13
    coeffs=(s.Rational(16*kk,3),s.Rational(1100*kk**3,9),s.Rational(392*kk*kk,9))
    print(kk,"monic complement coefficients:",coeffs)
print("PASS: trace, exact primitives, cup scaling, splitting and 13-integrality.")
print("NOT COMPUTED: normalized Coleman integrals, full height blocks or QC roots.")
