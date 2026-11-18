// The two global parity classes allowed by the Type 12 centers.
SetColumns(0);
SetSeed(2);
Q := Rationals();
R<x> := PolynomialRing(Q);
a := 160;
f := 160*x^3-240*x^2+215*x+90;
offsets := [Q|-3/2,-1/2,1/2,3/2];
moments := [&+[Binomial(2*k,2*j)*x^(k-j)*(&+[t^(2*j): t in offsets]):j in [0..k]]
            + ((10-x)/4)^k - (&+[t^(2*k): t in offsets])/2 : k in [1..3]];
e1 := moments[1];
e2 := (e1^2-moments[2])/2;
e3 := (moments[3]-e1*moments[2]+e2*moments[1])/3;
assert 64*e1 eq 240*x+320 and 64*e2 eq 320*x^2+280*x+436 and 64*e3 eq f;
E := EllipticCurve([Q|0,-240,0,215*a,90*a^2]);
P1 := E![160,-2400,1];
P9 := E![1440,50400,1];
gens,proved,full := Generators(E);
assert proved and full and #gens eq 2 and #TorsionSubgroup(E) eq 1;
assert &and[&or[P in {G,-G}: G in gens]: P in [P1,P9]];
K<w> := NumberField(x^3+38*x-504);
r := (4-w)/8;
assert IsIrreducible(f) and Evaluate(f,r) eq 0;
p := [pair[1]: pair in Decomposition(MaximalOrder(K),2) | pair[2] eq 2][1];
assert Valuation(w,p) eq 1 and Valuation(r,p) eq -5 and Valuation(K!a,p) eq 10;
print "Type 12: precisely the two single-generator classes can satisfy the dyadic parity condition.";

// The primary squareclass covers have full rational two-torsion.
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd) -> Integers() | O :-> Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd) -> Integers() | O :-> MinkowskiBound(O)>);
for b0 in [1,9] do
    d := b0-r;
    C := EllipticCurve([K|0,d*(10+r),0,10*r*d^2,0]);
    assert #TwoTorsionSubgroup(C) eq 4;
    upper := Valuation(#TwoSelmerGroup(C),2)-2;
    assert upper eq (b0 eq 1 select 3 else 4);
end for;
