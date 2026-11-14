#!/usr/bin/env python3
"""Certified Frobenius modulo 13^3; PROOF 5.23.
Exact rational arithmetic, existing SymPy; no modular working-roundoff.
The infinite-tail proof is in PROOF.md. This is not a height calculation.
"""
import sympy as s
x=s.symbols("x");p=13;MOD=p**3
expected={
-15:[[949,39,842,458],[481,1456,1678,583],[234,1248,1524,661],[1092,1391,833,465]],
-3:[[1248,871,933,2000],[1989,741,1422,1207],[741,351,673,1186],[1911,377,229,1732]]}
def residue(a):
    a=s.Rational(a);assert a.q%p
    return int(a.p)*pow(int(a.q),-1,MOD)%MOD
matrices={}
for k in [-15,-3]:
    f=s.Poly(x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k),x,domain=s.QQ)
    fp=f.diff()
    assert s.gcd(s.Poly(f,x,modulus=p),s.Poly(fp,x,modulus=p)).degree()==0
    inv=s.invert(fp,f)
    R=s.Poly(f.as_expr().subs(x,x**p),x)-f**p
    assert R.degree()<=64 and all(int(c)%p==0 for c in R.all_coeffs())
    def reduce(N,m):
        while m:
            B=(N*inv).rem(f)
            A=(N-B*fp).exquo(f)
            assert N==A*f+B*fp
            N=A+B.diff().mul_ground(s.Rational(2,2*m-1))
            m-=1
        while N.degree()>=4:
            e=N.degree();j=e-4
            c=N.LC()*s.Rational(2,2*e-3)
            G=s.Poly(x**j,x)*fp.mul_ground(s.Rational(1,2))
            if j:G+=s.Poly(j*x**(j-1),x)*f
            N-=G.mul_ground(c)
        return N
    cols=[]
    for i in range(4):
        N=s.Poly(0,x,domain=s.QQ)
        for n in range(4):
            a=(s.Poly(p*x**(p*(i+1)-1),x,domain=s.QQ)*R**n).mul_ground(s.binomial(-s.Rational(1,2),n))
            N+=reduce(a,6+13*n)
        cols.append([residue(N.nth(j)) for j in range(4)])
    M=s.Matrix.hstack(*[s.Matrix(c) for c in cols])
    assert M==s.Matrix(expected[k])
    T=s.symbols("T")
    assert all(int(c)%MOD==0 for c in s.Poly(M.charpoly(T).as_expr()-(T**4+6*T*T+169),T).all_coeffs())
    C=s.Matrix([[0,0,0,s.Rational(1,3)],[0,0,1,s.Rational(128*k,3)],
      [0,-1,0,-184*k*k],[-s.Rational(1,3),-s.Rational(128*k,3),184*k*k,0]])
    assert (M.T*C*M-p*C).applyfunc(residue)==s.zeros(4)
    assert (M-s.eye(4)).det()%p!=0
    a=(M-13*s.eye(4)).det()%MOD
    assert a%(p*p)==0 and a%MOD!=0
    # rank mod13 is exactly two: at least one unit 2-minor, all 3-minors zero.
    from itertools import combinations
    minors=lambda n:[M.extract(I,J).det()%p for I in combinations(range(4),n) for J in combinations(range(4),n)]
    assert any(minors(2)) and not any(minors(3))
    matrices[k]=M
    print(k,M)
D=s.diag(s.Rational(1,5),1,5,25)
assert (D*matrices[-15]+matrices[-3]*D).applyfunc(residue)==s.zeros(4)
print("PASS: certified matrices, characteristic polynomial, cup and twist checks.")
print("NOT COMPUTED: Frobenius primitives at endpoints, normalized mixed integrals or QC roots.")
