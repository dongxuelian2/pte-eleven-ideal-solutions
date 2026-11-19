// Both Type 12 parity classes. Run the companion PARI certificate first.
SetColumns(0);
SetSeed(1);
Q:=Rationals();
A<T>:=PolynomialRing(Q);
K<w>:=NumberField(T^3+38*T-504);
r:=(4-w)/8;
L:=NumberField(T^6+228*T^4+12996*T^2+7077920);
disc:=-5540347111683392000;
assert Discriminant(MaximalOrder(L)) eq disc;
function ProofBound(O)
    if AbsoluteDegree(NumberField(O)) eq 6 and Discriminant(O) eq disc then
        assert IsIsomorphic(NumberField(O),L);
        return 100;
    end if;
    return MinkowskiBound(O);
end function;
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd)->Integers() | O:->Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd)->Integers() | O:->ProofBound(O)>);
Z<X>:=PolynomialRing(K);
q:=X^2+(r-3/2)*X+r^2-3*r/2+43/32;
assert 160*(X-r)*q eq 160*X^3-240*X^2+215*X+90;
assert IsIsomorphic(AbsoluteField(ext<K|q>),L);
p:=907;
F:=GF(p);
R<z>:=PolynomialRing(F);
rs:=[(4-u[1])/8:u in Roots(z^3+38*z-504)];
assert #rs eq 3 and not IsSquare(F!-1);
for anchor in [1,9] do
    centers:={};
    for b in F do
        if not IsSquare(b) or not IsSquare(10-b) then continue;end if;
        if not &and[IsSquare((b-a)*(anchor-a)):a in rs] then continue;end if;
        pol:=64*z^3-(240*b+320)*z^2+(320*b^2+280*b+436)*z
             -(160*b^3-240*b^2+215*b+90);
        roots:=Roots(pol);
        if &+[Integers()|u[2]:u in roots] eq 3 and
           &and[IsSquare(u[1]):u in roots] then Include(~centers,b);end if;
    end for;
    assert centers eq {F!anchor};
end for;

procedure CheckDisk(E,P,slope)
    aa:=aInvariants(E)[2]; bb:=aInvariants(E)[4];
    D:=pAdicField(p,32);
    R<t>:=PolynomialRing(D);
    base:=[]; tangents:=[];
    for pair in Roots(t^3+38*t-504) do
        ev:=func<a|Evaluate(Polynomial(Q,Eltseq(K!a)),pair[1])>;
        ED:=EllipticCurve([D|0,ev(aa),0,ev(bb),0]);
        EP:=ED![ev(P[1]),ev(P[2]),1];
        assert Valuation(Discriminant(ED)) eq 0;
        EF:=EllipticCurve([F|0,Integers()!ev(aa),0,Integers()!ev(bb),0]);
        n:=#EF;
        T0:=n*EP;
        assert T0 ne ED!0 and Valuation(EP[2]) eq 0;
        h:=-T0[1]/T0[2];
        assert Valuation(h) ge 1 and AbsolutePrecision(h) ge 2;
        assert AbsolutePrecision(h/p) ge 1;
        Append(~base,F!(Integers()!(h/p)));
        Append(~tangents,F!(Integers()!(n*ev(slope)/(2*ev(P[2])))));
    end for;
    assert #base eq 3;
    assert &and[base[i]*tangents[j] ne base[j]*tangents[i]:
                i,j in [1..3]|i lt j];
    print "Log and tangent rows",base,tangents;
end procedure;

// b0 = 1.
d:=-160*(1-r);
aa:=d*Evaluate(Derivative(q),10);
bb:=d^2*Evaluate(q,10);
E:=EllipticCurve([K|0,aa,0,bb,0]);
xx:=-9*d;
ok,y:=IsSquare(xx*(xx^2+aa*xx+bb));
assert ok;
P:=E![xx,y,1];
assert Order(P) eq 0 and #TwoSelmerGroup(E) eq 16 and #TwoTorsionSubgroup(E) eq 2;
Cs,maps,ms:=TwoDescent(E:RemoveTorsion:=true,RemoveGens:={P});
G:=Domain(ms);
assert Invariants(G) eq [2,2] and CasselsTatePairing(ms(G.1),ms(G.2)) eq 1;

CheckDisk(E,P,d);

// b0 = 9. Use both center squares, then descend on the 2-isogenous curve.
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd)->Integers() | O:->Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd)->Integers() | O:->MinkowskiBound(O)>);
SetSeed(3);
d:=160*(9-r);
c:=10*d*Coefficient(q,0);
F0:=(X-c/10)*(X^2+10*d*Coefficient(q,1)*X+10*d*c);
assert Evaluate(F0,c*X) eq c^2*d*(10*X-1)*(1+Coefficient(q,1)*X+Coefficient(q,0)*X^2);
h:=Evaluate(F0,X+c/10);
E:=EllipticCurve([K|0,Coefficient(h,2),0,Coefficient(h,1),0]);
xx:=c/9-c/10;
ok,y:=IsSquare(Evaluate(h,xx));
assert ok;
P:=E![xx,y,1];
E1,phi:=IsogenyFromKernel(E,X);
Cs,maps,ms:=TwoDescent(E1:RemoveTorsion:=true,RemoveGens:={phi(P)});
G:=Domain(ms);
assert Invariants(G) eq [2,2] and CasselsTatePairing(ms(G.2),ms(G.1+G.2)) eq 1;
assert Order(P) eq 0 and #TwoSelmerGroup(E1) eq 16 and #TwoTorsionSubgroup(E1) eq 2;
CheckDisk(E,P,-c/81);
print "Both rank-one covers force b to equal their colliding anchor at every rational height.";
