#!/usr/bin/env python3
"""Exact mixed-divisor inputs and all good-prime finite terms; PROOF 5.18.
Requires existing SymPy. Does NOT compute normalized 13-adic integrals,
bad-prime horizontal intersections, a full height matrix or QC roots.
"""
import sympy as s
x=s.symbols("x")
BAD={2,3,5,11,19,53}
def red(a,p):
    a=s.cancel(a)
    return int(s.numer(a))*pow(int(s.denom(a)),-1,p)%p
def poly_mod(a,p):
    a=s.Poly(a,x,domain=s.QQ)
    return s.Poly.from_dict({i:red(c,p) for i,c in a.terms()},x,modulus=p)
data=[
(-15,1,s.Rational(225,4),s.Rational(320625,32),
 x*x+1020*x+81225,-13590*x-1564650,
 x*x-276*x+11025,406350-7794*x),
(-3,-1,s.Integer(33),s.Integer(10692),
 x*x+60*x+3249,-1026*x-40014,
 x*x+65*x+528,-1076*x-12804)]
expected=[
[{17:-2,593:-2},{},{7:2,313:2}],
[{},{},{29:-2,31:-2}]]
for row,ex in zip(data,expected):
    k,sign,a,b,uq,vq,ur,vr=row
    f=x*(x-k)*(x+5*k)*(x*x-68*k*x-275*k*k)
    us=[x-a,uq,ur]; vs=[b,vq,vr]
    assert s.factorint(abs(int(s.discriminant(f,x)))).keys()==BAD
    for u,v in zip(us,vs):
        assert s.rem(f-v*v,u,x)==0
        assert s.gcd(u,s.diff(u,x))==1
        assert s.gcd(u,v)==1
    V=sign*vq+(b-sign*vq.subs(x,a))/uq.subs(x,a)*uq
    quotient,rem=s.div(f-V*V,(x-a)*uq,x)
    assert rem==0
    assert s.Poly(quotient,x).monic().as_expr()==ur
    assert s.rem(-V,ur,x)==vr
    for (i,j),wanted in zip([(0,1),(0,2),(1,2)],ex):
        res=s.resultant(us[i],us[j],x)
        assert res!=0
        # Denominators are powers of two, so all odd models are integral.
        assert set(s.factorint(int(s.denom(res))))<={2}
        factors=s.factorint(abs(int(s.numer(res))))
        got={}
        for ell,e in factors.items():
            if ell in BAD or ell==13: continue
            assert e==1
            g=s.gcd(poly_mod(us[i],ell),poly_mod(us[j],ell))
            assert g.degree()==1
            root=-int(g.nth(0))*pow(int(g.nth(1)),-1,ell)%ell
            yi=red(vs[i].subs(x,root),ell)
            yj=red(vs[j].subs(x,root),ell)
            assert yi and yj
            same=1 if yi==yj else -1
            assert (yi-same*yj)%ell==0
            got[ell]=-2*same
        assert got==wanted,(k,i,j,got)
        print(k,(i,j),"resultant",res,"good local log coefficients",got)
    nums=[s.rem(v*s.diff(u,x),u,x) for u,v in zip(us,vs)]
    print("third-kind numerators:",nums)
    for u,v,n in zip(us,vs,nums):
        assert s.rem(n-v*s.diff(u,x),u,x)==0
        assert s.degree(n,x)<s.degree(u,x)
    # One ordinary geometric node; discriminant valuation gives thickness.
    for ell,node,thickness in [(11,0,2),(19,k,2),(53,34*k,1)]:
        ff=poly_mod(f,ell)
        common=s.gcd(ff,ff.diff())
        assert common.degree()==1
        assert common.eval(node)==0
        q,r=s.div(ff,poly_mod((x-node)**2,ell))
        assert r.is_zero and q.eval(node)!=0
        disc=abs(int(s.discriminant(f,x)));v=0
        while disc%ell==0:disc//=ell;v+=1
        assert v==thickness
        print("node",ell,node%ell,"thickness",v,"tangent",int(q.eval(node))%ell)
print("PASS. Missing bad-prime horizontal terms and normalized 13-adic integrals.")
