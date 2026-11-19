// Run the explicit-formula and PARI factor-base checks before these descents.
// Directions and Branches may be restricted to run individual descents.
if not assigned Directions then Directions := [false,true]; end if;
if not assigned Branches then Branches := [1,9,25]; end if;
SetColumns(0); SetSeed(3);
Q := Rationals(); A<x> := PolynomialRing(Q);
M := NumberField(x^12-6*x^11-66*x^10-1415*x^9-3987*x^8+128772*x^7+1543187*x^6
    -30854580*x^5+671675061*x^4-3249300991*x^3+46170279879*x^2-461672478855*x+1325918174175);
function ProofBound(O)
    if AbsoluteDegree(NumberField(O)) eq 12 and Discriminant(O) eq 5528109^5 then
        assert IsIsomorphic(NumberField(O),M);
        return 1000;
    end if;
    return MinkowskiBound(O);
end function;
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd)->Integers()|O:->ProofBound(O)>);
assert Discriminant(MaximalOrder(M)) eq 5528109^5 and Signature(M) eq 0;
assert AbelianInvariants(ClassGroup(MaximalOrder(M))) eq [3];
K<h> := NumberField(x^4-x^3+28*x^2-54*x+231);
r := (-2*h^3+h^2-33*h+129)/15;
R<u> := PolynomialRing(K);
f := u^4-18*u^3+114*u^2-97*u+405;
q := ExactQuotient(f,u-r);
assert IsIrreducible(q) and IsIsomorphic(AbsoluteField(ext<K|q>),M);
for reciprocal in Directions do
    g := reciprocal select ReciprocalPolynomial(q)/Coefficient(q,0) else q;
    assert IsMonic(g);
    for B0 in Branches do
        t := reciprocal select 1/B0 else B0;
        k := Evaluate(g,t);
        E := EllipticCurve([K|0,k*Coefficient(g,2),0,k^2*Coefficient(g,1),k^3*Coefficient(g,0)]);
        S := TwoSelmerGroup(E);
        dim := reciprocal select 3 else (B0 eq 9 select 4 else 5);
        assert Invariants(S) eq [2:i in [1..dim]];
        print "B0", B0, "reciprocal", reciprocal, "2-Selmer dimension", dim;
    end for;
end for;
