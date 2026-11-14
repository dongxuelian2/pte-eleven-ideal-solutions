#!/usr/bin/env python3
"""Regular-model mixed heights at 11,19,53; see PROOF.md Section 5.19.
Requires adjacent check_genus_two_mixed_heights.py and existing SymPy.
Does not compute primes 2,3,5, normalized 13-adic integrals or QC roots.
"""
from check_genus_two_mixed_heights import data, x, red, poly_mod
import sympy as s
from itertools import combinations
from math import gcd
blowups={
(-15,11):[],(-15,19):[(0,18,9),(2,6,4)],
(-3,11):[(0,3,4),(2,4,10)],(-3,19):[(1,16,13),(2,10,3)]}
shared_expected={
(-15,11):[],(-15,19):[(0,2,4,0,0),(1,2,6,8,11)],(-15,53):[],
(-3,11):[(0,2,0,0,0)],(-3,19):[(0,2,14,14,5),(1,2,16,0,0)],(-3,53):[]}
lengths={
-15:[[(0,0,0),(0,0,0)],[(0,1,0),(0,1,0)],[(0,0,0),(0,1,0)]],
-3:[[(0,0,0),(0,0,0)],[(1,0,0),(1,1,0)],[(0,1,0),(0,1,0)]]}
def val(q,p):
    q=s.Rational(q);assert q
    a,b,v=int(q.p),int(q.q),0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v
for k,sign,a,b,uq,vq,ur,vr in data:
    us=[x-a,uq,ur];vs=[b,vq,vr]
    for p in [11,19,53]:
        shared=[]
        for i,j in [(0,1),(0,2),(1,2)]:
            g=s.gcd(poly_mod(us[i],p),poly_mod(us[j],p))
            assert g.degree()<=1
            for root in range(p):
                if red(us[i].subs(x,root),p) or red(us[j].subs(x,root),p):continue
                yi,yj=red(vs[i].subs(x,root),p),red(vs[j].subs(x,root),p)
                shared.append((i,j,root,yi,yj))
                if yi and yj:
                    assert yi==(-yj)%p
                    assert val(s.resultant(us[i],us[j],x),p)==1
        assert shared==shared_expected[k,p]
    for p in [11,19]:
        node=0 if p==11 else k
        got=[]
        for i,(u,v) in enumerate(zip(us,vs)):
            if red(u.subs(x,node),p):continue
            d=red(s.diff(u,x).subs(x,node),p*p);assert d%p
            h=(node-red(u.subs(x,node),p*p)*pow(d,-1,p*p))%(p*p)
            y=red(v.subs(x,h),p*p)
            assert red(u.subs(x,h),p*p)==0 and y%p==0
            got.append((i,((h-node)//p)%p,(y//p)%p))
        assert got==blowups[k,p]
        if len(got)==2:assert got[0][1]!=got[1][1]
        print(k,p,"resolved node coordinates",got)
    # Independent naive-model module lengths. At nodes these are NOT
    # automatically regular-model intersection numbers; blowups above decide.
    for pair,(i,j) in enumerate([(0,1),(0,2),(1,2)]):
        d=s.degree(us[i],x);actual=[]
        for eps in [1,-1]:
            cols=[]
            for h in [us[j],vs[i]-eps*vs[j]]:
                for t in range(d):
                    r=s.Poly(s.rem(x**t*h,us[i],x),x)
                    cols.append(s.Matrix([r.nth(n) for n in range(d)]))
            M=s.Matrix.hstack(*cols)
            minors=[M[:,list(cc)].det() for cc in combinations(range(2*d),d)]
            # All denominators are powers of two, units at the primes used.
            assert all(set(s.factorint(int(s.denom(q))))<={2} for q in minors)
            g=gcd(*(abs(int(s.numer(q))) for q in minors));assert g
            actual.append(tuple(val(g,p) for p in [11,19,53]))
        assert actual==lengths[k][pair]
print("PASS: all mixed finite contributions outside {2,3,5,13} are certified.")
print("Uncomputed: primes 2,3,5; normalized 13-adic integrals; full QC and roots.")
