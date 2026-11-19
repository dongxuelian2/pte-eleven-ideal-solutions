// Enforce a common sign at all three places over 131.
SetColumns(0);
Q := Rationals(); A<x> := PolynomialRing(Q);
pol := x^4-x^3+28*x^2-54*x+231; K<h> := NumberField(pol);
R<u> := PolynomialRing(K); r := (-2*h^3+h^2-33*h+129)/15;
f := x^4-18*x^3+114*x^2-97*x+405;
q := ExactQuotient(R!f,u-r); g := ReciprocalPolynomial(q)/Coefficient(q,0);
t := Q!1/9; k := Evaluate(g,t);
E := EllipticCurve([K|0,k*Coefficient(g,2),0,k^2*Coefficient(g,1),k^3*Coefficient(g,0)]);
P0 := E![k*t,k^2,1];
t3 := (1079145271477+193174581916*h+69465846023*h^2-2741758426*h^3)/3034235677725;
G := [P0] cat [E![k*z,Sqrt(k^3*Evaluate(g,z)),1]:z in [K|-59/64,t3]];

for p in [17,131,197,269] do
    F := GF(p); Rp<X> := PolynomialRing(F); factors := Factorization(Rp!f);
    assert IsSquarefree(Rp!f) and p mod 20 ne 1;
    seeds := {};
    for B in {a^2:a in F} do
        P := (X-B)*((X-B-4)^2-16*B)*((X-B-16)^2-64*B);
        roots := Roots(ExactQuotient(P-Evaluate(P,1),X-1));
        if &+[Integers()|v[2]:v in roots] eq 4 and &and[IsSquare(v[1]):v in roots]
           and &and[IsSquare(Evaluate(v[1],B)*Evaluate(v[1],9)):v in factors] then
            Include(~seeds,B);
        end if;
    end for;
    assert seeds eq {F!9};
end for;

function Saturate(points,ell)
    assert #DivisionPoints(E!0,ell) eq 1;
    coeffs := ell eq 2 select [0,1] else [-1,0,1];
    while true do
        changed := false;
        for a,b,c in coeffs do
            v := [a,b,c]; ix := [i:i in [1..3]|v[i] ne 0];
            if IsEmpty(ix) or v[ix[#ix]] ne 1 then continue; end if;
            halves := DivisionPoints(&+[v[i]*points[i]:i in ix],ell);
            if not IsEmpty(halves) then
                points[ix[#ix]] := halves[1]; changed := true; break a;
            end if;
        end for;
        if not changed then return points; end if;
    end while;
end function;

signs := {<a,b,c>:a,b,c in [-1,1]};
for ell in [2,3] do
    primePowers := ell eq 2 select 4 else 2; m := ell^primePowers;
    primes := ell eq 2 select [17,197,131] else [269,131];
    basis := Saturate(G,ell);
    candidates := [<a,b,c>:a,b,c in [0..m-1]];
    at131 := [* *];
    for p in primes do
        F := GF(p); Rp<z> := PolynomialRing(F);
        fac := Factorization(Rp!pol);
        if p eq 131 then assert [Degree(v[1]):v in fac] eq [1,1,2]; end if;
        for v in fac do
            d := v[1];
            if Degree(d) eq 1 then L := F; root := -Coefficient(d,0);
            else L<root> := ext<F|d>; end if;
            image := func<a | Evaluate(Polynomial([L!c:c in Eltseq(K!a)]),root)>;
            EL := EllipticCurve([L|image(a):a in aInvariants(E)]);
            N := #EL; multiple := N div GCD(N,m);
            pts := [multiple*(EL![image(P[1]),image(P[2]),1]):P in basis cat [P0]];
            multiples := [[a*pts[j]:a in [0..m-1]]:j in [1..3]];
            allowed := AssociativeArray();
            for c in candidates do
                pt := &+[multiples[j][1+c[j]]:j in [1..3]];
                allowed[c] := {s:s in [-1,1]|pt eq s*pts[4]};
            end for;
            candidates := [c:c in candidates|not IsEmpty(allowed[c])];
            if p eq 131 then Append(~at131,allowed); end if;
        end for;
    end for;
    assert #candidates eq (ell eq 2 select 4 else 36);
    signs := {s:s in signs | &or[&and[s[j] in at131[j][c]:j in [1..3]]:c in candidates]};
end for;
assert signs eq {<1,1,1>,<-1,-1,-1>};
print "The simultaneous reductions force a common sign above 131.";
