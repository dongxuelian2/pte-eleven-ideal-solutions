// Unconditional 2-Selmer rank bounds for the three Type 11 genus-three covers.
SetColumns(0);
SetSeed(1);
Q := Rationals();
R<r> := PolynomialRing(Q);
search := map<PowerStructure(RngOrd) -> Integers() | O :-> Min(1000,MinkowskiBound(O))>;
proof := map<PowerStructure(RngOrd) -> Integers() | O :-> MinkowskiBound(O)>;
SetClassGroupBoundMaps(search,proof);

d := r^2-11*r-1;
a := r^2-261*r-1;
b := 7*r^2-2*r+18;
c := 9*r^2-69*r-49;
for f in [r*d*a*b,r*d*b*c,r*a*b*c] do
    assert Degree(f) eq 7 and IsSquarefree(f);
    assert Sort([Degree(term[1]) : term in Factorization(f)]) eq [1,2,2,2];
    H := HyperellipticCurve(f);
    S := TwoSelmerGroup(Jacobian(H) : Al:="TwoSelmerGroupOld");
    assert AbelianInvariants(S) eq [2,2,2,2,2];
    // Four irreducible factors give three independent rational 2-torsion classes.
    print "Type 11 odd-degree genus-three cover, rank at most",Ngens(S)-3;
end for;
print "PASS: all three rank bounds";
