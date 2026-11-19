// Run the explicit-formula and factor-base certificates for Type 18 first.
SetColumns(0);
SetSeed(3);
Q := Rationals();
A<x> := PolynomialRing(Q);
M := NumberField(x^12-3*x^11-5625*x^10+44481*x^9+9843099*x^8-99843252*x^7
    -5998855201*x^6+67870752036*x^5+854290494639*x^4-9750057898095*x^3
    -2947764849192*x^2+213431600705238*x-420745260940011);
function ProofBound(O)
    if AbsoluteDegree(NumberField(O)) eq 12 and Discriminant(O) eq 5^6*11^9*1231^5*2417^5 then
        assert IsIsomorphic(NumberField(O), M);
        return 1000;
    end if;
    return MinkowskiBound(O);
end function;
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd)->Integers() | O :-> Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd)->Integers() | O :-> ProofBound(O)>);
K<h> := NumberField(x^4-x^3-991*x^2+13345*x-19455);
r := (-h^3-23*h^2+585*h+5002)/73;
R<u> := PolynomialRing(K);
f := u^4-165*u^3+8778*u^2-172810*u+1057221;
q := ExactQuotient(f,u-r);
k := 25-r;
E := EllipticCurve([K|0,k*Coefficient(q,2),0,k^2*Coefficient(q,1),k^3*Coefficient(q,0)]);
O := MaximalOrder(M);
assert Signature(M) eq 12 and Discriminant(O) eq 5^6*11^9*1231^5*2417^5;
assert AbelianInvariants(ClassGroup(O)) eq [2];
assert Evaluate(f,r) eq 0 and IsIrreducible(q);
assert IsIsomorphic(AbsoluteField(ext<K|q>), M);
assert Invariants(TwoSelmerGroup(E)) eq [2,2,2];
print "The class group has order two and the 2-Selmer group has dimension three.";
