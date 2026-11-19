// Logarithm spans and the rational tangent for both reciprocal covers.
SetColumns(0);
Q := Rationals(); A<x> := PolynomialRing(Q);
K<h> := NumberField(x^4-x^3+28*x^2-54*x+231);
r := (-2*h^3+h^2-33*h+129)/15;
R<u> := PolynomialRing(K);
f := u^4-18*u^3+114*u^2-97*u+405;
q := ExactQuotient(f,u-r);
g := ReciprocalPolynomial(q)/Coefficient(q,0);

for B0 in [9,25] do
    t0 := Q!1/B0; k := Evaluate(g,t0);
    if B0 eq 9 then
        p := 131;
        ts := [K|t0,-59/64,
            (1079145271477+193174581916*h+69465846023*h^2-2741758426*h^3)/3034235677725];
    else
        p := 2;
        ts := [K|t0,(40557875124-6111404418*h+447986201*h^2-437817952*h^3)/24006000375,
            (-664143+293451*h-67207*h^2+2264*h^3)/1012500];
    end if;
    E := EllipticCurve([K|0,k*Coefficient(g,2),0,k^2*Coefficient(g,1),k^3*Coefficient(g,0)]);
    pts := [E![k*t,Sqrt(k^3*Evaluate(g,t)),1] : t in ts];
    places := Decomposition(MaximalOrder(K),p);
    degrees := Sort([Valuation(Norm(v[1]),p):v in places]);
    assert degrees eq (p eq 2 select [4] else [1,1,2]);
    completions := [* *]; axis := []; m := p eq 2 select 2^26*3*5*7*17 else p^10;
    for v in places do
        assert v[2] eq 1;
        L,em := Completion(K,v[1]:Precision:=192);
        O := Integers(L); F,red := ResidueClassField(O);
        if Degree(L) gt 1 then
            assert Degree(MinimalPolynomial(red(O!L.1))) eq Degree(L);
        end if;
        assert Valuation(em(k)) eq 0 and &and[Valuation(em(c)) ge 0:c in Coefficients(g)];
        EL := EllipticCurve([L|em(a):a in aInvariants(E)]);
        if p ne 2 then
            EF := EllipticCurve([F|red(O!em(a)):a in aInvariants(E)]);
            m *:= #EF;
        end if;
        Append(~completions,<em,EL>);
        coeff := Degree(L) eq 1 select [em(1/k)] else Eltseq(em(1/k));
        axis cat:= [Integers()!c mod p^5:c in coeff];
    end for;
    rows := [];
    for P in pts do
        coordinates := [* *];
        for v in completions do
            em,EL := Explode(v);
            H := m*(EL![em(P[1]),em(P[2]),1]); assert H ne EL!0;
            z := -H[1]/H[2]; assert Valuation(z) ge 8;
            coeff := Degree(Parent(z)) eq 1 select [z] else Eltseq(z);
            for a in coeff do Append(~coordinates,a); end for;
        end for;
        depth := Minimum([Valuation(a):a in coordinates]);
        assert &and[AbsolutePrecision(a) gt depth+5:a in coordinates];
        Append(~rows,[Integers()!(a/p^depth) mod p^5:a in coordinates]);
    end for;
    minors := [Determinant(Matrix(Integers(),[[v[i]:i in [1..4]|i ne j]:v in rows]))
               : j in [1..4]];
    common := Minimum([Valuation(a,p):a in minors|a ne 0]);
    assert common eq (p eq 2 select 1 else 0);
    assert Valuation(Determinant(Matrix(Integers(),rows cat [axis])),p) eq common;
end for;
print "Both rank-three logarithm spans are transverse to their rational parameter disks.";
