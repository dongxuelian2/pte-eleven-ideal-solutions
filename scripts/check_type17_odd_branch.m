// Single-root cubic ranks for Types 1 and 17; the two Type 17 parity covers.
SetColumns(0); SetSeed(1);
Q := Rationals(); R<x> := PolynomialRing(Q);
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd) -> Integers() | O :-> Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd) -> Integers() | O :-> MinkowskiBound(O)>);
PP<U,b,h> := ProjectiveSpace(Q,2);
cubics := [
    64*U^3-(1152*b+1328*h)*U^2+(7232*b^2+14080*b*h+7564*h^2)*U
    -16384*b^3-38672*b^2*h-33144*b*h^2-11025*h^3,
    64*U^3+(48*b-1328*h)*U^2+(32*b^2-920*b*h+7564*h^2)*U
    +16*b^3-472*b^2*h+4131*b*h^2-11025*h^3
];
basepoints := [[Q | 25/4,1,1],[Q | 1/4,16,1]];
models := [[Q | 1,0,0,-51865,3207225],[Q | 0,0,0,-3207,50506]];
for i in [1..2] do
    C := Curve(PP,cubics[i]);
    assert IsNonsingular(C);
    E := EllipticCurve(C,C!basepoints[i]);
    assert aInvariants(MinimalModel(E)) eq models[i];
    lo,hi := RankBounds(E); assert lo eq 4 and hi eq 4;
end for;

f := 11025*x^3-4131*x^2+472*x-16;
E0 := EllipticCurve([Q | 0,-4131,0,11025*472,-16*11025^2]);
G := E0![11025/16,165375/64,1];
gens,rank_proved,full := Generators(E0);
assert rank_proved and full and #gens eq 1 and gens[1] in {G,-G};
assert #TorsionSubgroup(E0) eq 1;

K<w> := NumberField(x^3+x^2-240*x+1120);
r := (2*w^2+5*w+140)/3675;
assert IsIrreducible(f) and Evaluate(f,r) eq 0;
v := Q!1/16;
for d in [K!1, v-r] do
    g := EllipticCurve([K | 0,d*(5-r),0,-5*r*d^2,0]);
    sel := TwoSelmerGroup(g);
    assert #TwoTorsionSubgroup(g) eq 4;
    if d eq 1 then assert #sel eq 32; continue; end if;
    assert #sel eq 8;
    E := MinimalModel(g);
    yes,iso := IsIsomorphic(g,E); assert yes;
    P := iso(g![d*v,9*d^2/16,1]);
    assert Order(P) eq 0;
    tors,inc := TorsionSubgroup(E); assert AbelianInvariants(tors) eq [2,2];

    // For each translate of P by 2-torsion, a split place over 17
    // witnesses a nonsquare difference from a branch abscissa.
    witnessed := {};
    for root in [6,11] do
        F := GF(17); ww := F!root;
        assert ww^3+ww^2-240*ww+1120 eq 0;
        rr := (2*ww^2+5*ww+140)/3675; dd := F!v-rr;
        Ebar := EllipticCurve([F | 0,dd*(5-rr),0,-5*rr*dd^2,0]);
        roots := [F | 0,-5*dd,rr*dd]; assert #Seqset(roots) eq 3;
        translates := [Ebar!0] cat [Ebar![e,0,1] : e in roots];
        U := Ebar![dd*v,9*dd^2/16,1];
        for i in [1..4] do
            A := U+translates[i];
            if A ne Ebar!0 and exists{e : e in roots | not IsSquare(A[1]-e)} then
                Include(~witnessed,i);
            end if;
        end for;
    end for;
    assert witnessed eq {1,2,3,4};
    group := AbelianGroup([2,2,0]);
    mp := map<group->E | a :-> Eltseq(a)[1]*inc(tors.1)
                             +Eltseq(a)[2]*inc(tors.2)+Eltseq(a)[3]*P>;
    P1 := ProjectiveSpace(Q,1);
    oldmap := map<g->P1 | [g.1,d*g.3]>;
    phi := Extend(Expand(Inverse(iso)*oldmap));
    V,idx := Chabauty(mp,phi : IndexBound:=2);
    assert #V eq 5 and IsPowerOf(idx,2);
    assert {phi(mp(a)): a in V} eq {P1 | [1,0],[0,1],[-5,1],[1/16,1]};
    print "Odd coset: five points; index obstruction",idx;
end for;
print "Type 17: every valid lift would map to an even multiple on E0.";

// If the product point is 2Q, the two squared linear factors at 2Q
// give a genus-three necessary condition on the abscissa of Q.
a := Q!11025; b2 := Q!-4131; c := Q!472; e := Q!-16;
F := a^2*x^4-2*a*c*x^2-8*a*e*x+c^2-4*b2*e;
H := 4*a*f;
FF := FieldOfFractions(R); u := FF!x;
assert (FF!F)/H eq Evaluate(Derivative(f),u)^2/(4*a*Evaluate(f,u))-b2/a-2*u;
C := HyperellipticCurve(F*(F+5*H));
assert Genus(C) eq 3 and IsSquarefree(F*(F+5*H));
assert RankBound(Jacobian(C)) eq 2;
print "Even coset: the genus-three halving curve has Jacobian rank at most two.";
