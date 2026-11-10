// Software-assisted upper bounds; run in a fresh Magma session.
// The independent lower bounds are check_shift_curve_rank_obstructions.py.
// No GRH assumption: every nontrivial descent field is quadratic, and
// its unconditional Minkowski bound is checked below against 1000.
SetColumns(0);
SetSeed(1);
SetClassGroupBounds(1000);
Q := Rationals();
P<x> := PolynomialRing(Q);
cases := [
    <4,3*x*(240*x^2+40*x+167)*(256*x^2+97*x+466),32,4,[-231,-155]>,
    <5,x*(80*x^2-840*x+589)*(4*x-49)*(4*x-9),64,8,[505]>,
    <15,2*x*(2*x-23)*(4*x-31)*(160*x^2-20*x+49),128,8,[-215]>,
    <16,6*x*(x-6)*(4*x-9)*(160*x^2-1335*x+2304),32,8,[34185]>
];
for data in cases do
    typ,f,selmer_order,torsion_order,expected_discs := Explode(data);
    discs := [Integers() | ];
    for factor in Factorization(f) do
        g,e := Explode(factor);
        assert e eq 1 and Degree(g) in {1,2};
        if Degree(g) eq 2 then
            K := NumberField(g/LeadingCoefficient(g));
            assert MinkowskiBound(K) lt 1000;
            Append(~discs,Discriminant(MaximalOrder(K)));
        end if;
    end for;
    Sort(~discs);
    assert discs eq expected_discs;
    J := Jacobian(HyperellipticCurve(f));
    S := TwoSelmerGroup(J : Al := "TwoSelmerGroupOld");
    T := TwoTorsionSubgroup(J);
    assert #S eq selmer_order and #T eq torsion_order;
    print "TYPE",typ,"FIELD_DISCRIMINANTS",discs,
          "SELMER_ORDER",#S,"TWO_TORSION_ORDER",#T;
end for;
print "PASS: exact two-Selmer upper bounds with unconditional class-group bounds.";
