// Proposition 27. Run the prime-power and factor-base checks in type13_prym_class first.
SetColumns(0);
SetSeed(1);
Q := Rationals();
R<T> := PolynomialRing(Q);
K<th> := NumberField(T^3-20*T^2/3+13*T/12-95/48);
P<x> := PolynomialRing(K);
f := x^5*(908800*th^2/27-18444800*th/81+7356800/81)
    + x^4*(-108800*th^2/9+2840320*th/27-743680/27)
    + x^3*(2800*th^2-133760*th/9+21740/9)
    + x^2*(-980*th^2/9+6520*th/9-50/9)
    + x*(-25*th^2/24-253*th/24-247/96);
g := ExactQuotient(f,x);
g := g/LeadingCoefficient(g);
assert IsSquarefree(f) and IsIrreducible(g);
Rel := ext<K|g>;
Abs := AbsoluteField(Rel);
L := NumberField(T^12-4*T^11-68*T^10-18698*T^9+336998*T^8-500324*T^7
    +7039932*T^6-194015520*T^5+1041300636*T^4-3497537936*T^3
    +20028904184*T^2-109223576168*T+269574001992);
assert IsIsomorphic(Abs,L);
nreal,ncomplex := Signature(L);
assert [nreal,ncomplex] eq [2,5];
OL := MaximalOrder(L);
assert Discriminant(OL) eq -172707117717588564754287819295488000000;
assert MinkowskiBound(MaximalOrder(K)) eq 419;
function CertifiedBound(order, reference)
    field := NumberField(order);
    if AbsoluteDegree(field) eq 12 and Discriminant(order) eq Discriminant(MaximalOrder(reference)) then
        assert IsIsomorphic(field,reference);
        return 1000;
    end if;
    return MinkowskiBound(order);
end function;
search := map<PowerStructure(RngOrd)->Integers() | order:->Min(1000,MinkowskiBound(order))>;
proof := map<PowerStructure(RngOrd)->Integers() | order:->CertifiedBound(order,L)>;
SetClassGroupBoundMaps(search,proof);
assert Invariants(ClassGroup(OL)) eq [2,2];
assert IsTrivial(ClassGroup(MaximalOrder(K)));
J := Jacobian(HyperellipticCurve(f));
sel := TwoSelmerGroup(J : Al := "TwoSelmerGroupOld", Fields := {K,L});
tors := TwoTorsionSubgroup(J);
assert Invariants(sel) eq [2,2,2,2,2,2] and Invariants(tors) eq [2];
print "Class group (2,2), 2-Selmer dimension 6, rational 2-torsion dimension 1";
