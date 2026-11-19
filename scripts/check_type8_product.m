// Type 8 product covers and the three possible two-cover classes.
SetColumns(0); SetSeed(4);
Q := Rationals(); A<x> := PolynomialRing(Q);
SetClassGroupBoundMaps(map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>,
                       map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>);
g := 5*(256*x^4-1152*x^3+1824*x^2-388*x+405); f := x*g;
C1 := HyperellipticCurve(g); E := EllipticCurve(C1,Points(C1,1/4)[1]);
M := MinimalModel(E); assert aInvariants(M) eq [Q|0,0,0,-105150,7613125];
lo,hi := RankBounds(M); assert lo eq 5 and hi eq 5;
C := HyperellipticCurve(f); J := Jacobian(C);
lo,hi := RankBounds(J); assert hi eq 4 and Invariants(TorsionSubgroup(J)) eq [2];
H,eta := TwoCoverDescent(C); L<t> := Domain(eta); G := Codomain(eta);
bs := [Q|1/4,3/4,9/4,4,25/4]; assert &and[IsSquare(Evaluate(f,b)):b in bs];
assert Invariants(sub<G|[eta(b-t):b in bs]>) eq [2,2,2,2,2];
target := {eta(i/4-t):i in [1,9,25]}; assert #H eq 7 and #target eq 3;

function SplitAt(b)
    F := Parent(b); R<u> := PolynomialRing(F);
    w := Sqrt(b); B := &*[u-(w+j)^2:j in [-2..2]];
    K := (B-Evaluate(B,F!(1/4))) div (u-1/4);
    roots := Roots(K);
    return &+[Integers()|p[2]:p in roots] eq 4 and &and[IsSquare(p[1]):p in roots];
end function;

function Survives(h,eta,f,p)
    F := GF(p); R<v> := PolynomialRing(F); fp := R!f;
    assert p mod 20 ne 1 and IsSquarefree(fp);
    hp := Polynomial(Q,Eltseq(h@@eta));
    hp *:= LCM([Denominator(c):c in Coefficients(hp)]);
    hp /:= GCD([Integers()!c:c in Coefficients(hp)]);
    ep := R!hp; if Degree(GCD(ep,fp)) gt 0 then return true; end if;
    for b in F do
        if not IsSquare(b) or not SplitAt(b) then continue; end if;
        odd := {}; ok := true;
        for pair in Factorization(fp) do
            g := pair[1]/LeadingCoefficient(pair[1]);
            if Evaluate(g,b) eq 0 then continue; end if;
            sq := IsSquare(Evaluate(g,b)*Resultant(g,ep));
            if IsOdd(Degree(g)) then Include(~odd,sq);
            elif not sq then ok := false; end if;
        end for;
        if ok and #odd le 1 then return true; end if;
    end for;
    return false;
end function;

for pair in [[13,6],[17,5],[19,4],[43,3]] do
    p,n := Explode(pair); H := {h:h in H|Survives(h,eta,f,p)}; assert #H eq n;
end for;
assert H eq target;
print "Type 8 product ranks 5 and 4, and three global two-cover classes verified.";
