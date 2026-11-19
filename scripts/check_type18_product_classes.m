// The class of the discriminant parameter and the five marked vertex classes.
SetColumns(0); SetSeed(2);
Q := Rationals(); A<x> := PolynomialRing(Q);
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>);
H := x^4-165*x^3+8778*x^2-172810*x+1057221;
seeds := [1,9,25,49,81]; U0 := 945^2;
assert x*H-U0 eq &*[x-a:a in seeds];
S<T> := PolynomialRing(A);
d := Discriminant(T^5-165*T^4+8778*T^3-172810*T^2+1057221*T-x);
assert Degree(d) eq 4 and IsSquarefree(x*d);
primes := [13,17,19,29,37,73,89,139,193];
for p in primes do
    F := GF(p); Rp<X> := PolynomialRing(F);
    assert p mod 10 ne 1 and #{F!s:s in seeds} eq 5;
    possible := {};
    for u in {a^2:a in F} do
        roots := Roots(X*(Rp!H)-u);
        if &+[Integers()|v[2]:v in roots] eq 5 and &and[IsSquare(v[1]):v in roots] then
            Include(~possible,u);
        end if;
    end for;
    assert possible eq {F!U0};
end for;

function Compatible(s,mu,f,z,p)
    R<U> := PolynomialRing(GF(p)); fp := R!f;
    assert IsSquarefree(U*fp) and Evaluate(fp,z) ne 0 and z mod p ne 0;
    hh := Polynomial(Q,Eltseq(s@@mu));
    hh /:= Coefficient(hh,0);
    if &or[Valuation(c,p) lt 0:c in Coefficients(hh)] then return true; end if;
    h := R!hh; if Degree(GCD(h,fp)) gt 0 then return true; end if;
    return &and[IsSquare(Resultant(v[1],h*(z-U)/z)):v in Factorization(fp)];
end function;

for which in [1,2] do
    f := which eq 1 select d else H;
    parameters := which eq 1 select [U0] else seeds;
    ps := which eq 1 select [13,19] else [13,17,29,37,73,89,139,193];
    counts := which eq 1 select [2,1] else [30,22,14,13,10,8,6,5];
    C := HyperellipticCurve(x*f);
    ss,mu := TwoCoverDescent(C); et := Domain(mu);
    assert #ss eq (which eq 1 select 5 else 38);
    for i in [1..#ps] do
        ss := {s:s in ss|&or[Compatible(s,mu,f,z,ps[i]):z in parameters]};
        assert #ss eq counts[i];
    end for;
    want := {mu(z-et.1):z in parameters};
    assert ss eq want and #want eq #parameters;
    if which eq 2 then
        for z in seeds do
            assert {s:s in ss|Compatible(s,mu,f,z,89)} eq {mu(z-et.1)};
        end for;
    end if;
end for;
print "The discriminant class is fixed and each of the five vertex classes occurs once.";
