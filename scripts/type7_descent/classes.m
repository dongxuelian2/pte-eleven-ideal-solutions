// Global squareclasses and the residue disks used in the Type 7 descent.
SetColumns(0);
SetSeed(3);
Q := Rationals();
SetClassGroupBoundMaps(map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>,
                       map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>);
Z := Integers();
A<x> := PolynomialRing(Q);
Base := EllipticCurve([Q|0,-395,0,38475,-820125]);
anchors := [Q|0,2,1/3,3,2/5]; // The first entry denotes the origin.
twists := [15,2,3,3,10];

function BaseCoset(p,t0)
    E := ChangeRing(Base,GF(p));
    if t0 eq 0 then R := E!0; else R := Points(E,135*t0)[1]; end if;
    return {R + 216*P : P in Points(E)};
end function;

function AllowsTwist(P,d,p)
    F := GF(p);
    if P[3] eq 0 then return IsSquare(F!(15*d)); end if;
    t := P[1]/135;
    if t ne 0 then return IsSquare(t/d); end if;
    R<u> := PolynomialRing(F);
    return exists{v : v in Roots(27*u^2+168*u+288) | IsSquare(v[1]/d)};
end function;

for i in [1..5] do
    ds := {1,2,3,5,6,10,15,30};
    for p in [11,17,19,41] do
        coset := BaseCoset(p,anchors[i]);
        ds := {d : d in ds | exists{P : P in coset | AllowsTwist(P,d,p)}};
    end for;
    assert ds eq {twists[i]};
end for;

function DiskParameters(p,d,t0)
    F := GF(p); R<v> := PolynomialRing(F); good := {F|};
    for P in BaseCoset(p,t0) do
        assert P[3] ne 0;
        t := P[1]/135;
        if not IsSquare(t/d) then continue; end if;
        aa := 45*t^3-145*t^2+135*t-27;
        bb := -8*(45*t^2-80*t+21); cc := 144*(5*t-2);
        if aa eq 0 then Include(~good,t); continue; end if;
        for ur in Roots(aa*v^2+bb*v+cc) do
            u := ur[1]; a := 1+u; b := t*u; c := (25+(10-6*t)*u)/4;
            if not &and[IsSquare(h):h in [a,b,c]] then continue; end if;
            ss := (16*a+20*c+81)/12;
            pp := -(128*a^2-640*a*c+48*a+80*c^2-1380*c-1701)/144;
            roots := Roots(v^2-ss*v+pp);
            if &+[Z|e[2]:e in roots] eq 2 and &and[IsSquare(e[1]):e in roots] then
                Include(~good,t);
            end if;
        end for;
    end for;
    return good;
end function;

for row in [[2,17],[3,17],[3,83],[4,61],[5,389]] do
    i,p := Explode(row);
    assert DiskParameters(p,twists[i],anchors[i]) eq {GF(p)!anchors[i]};
end for;

function LocalFakeClass(f,h,z,p)
    F := GF(p); R<v> := PolynomialRing(F); fp := R!f;
    assert IsSquarefree(fp) and Evaluate(fp,z) ne 0;
    hp := Polynomial(Q,Eltseq(h));
    hp *:= LCM([Denominator(c):c in Coefficients(hp)]);
    // Remove a common rational power of p before testing the component norms.
    hp /:= GCD([Z!c:c in Coefficients(hp)]);
    ep := R!hp;
    if Degree(GCD(fp,ep)) gt 0 then return true; end if;
    odd := {};
    for pair in Factorization(fp) do
        g := pair[1]/LeadingCoefficient(pair[1]);
        sq := IsSquare(Evaluate(g,z)*Resultant(g,ep));
        if IsOdd(Degree(g)) then Include(~odd,sq);
        elif not sq then return false; end if;
    end for;
    return #odd le 1;
end function;

for row in [[2,1,1],[3,1,3]] do
    d,v,w := Explode(row); z0 := Q!(v/w);
    f := 5*(27*d^3*x^6-79*d^2*x^4+57*d*x^2-9);
    H,eta := TwoCoverDescent(HyperellipticCurve(f)); L<e> := Domain(eta);
    assert #H eq (d eq 2 select 4 else 8);
    keep := {h:h in H | exists{z:z in [z0,-z0] | LocalFakeClass(f,h@@eta,GF(17)!z,17)}};
    assert keep eq {eta(z0-e),eta(-z0-e)} and #keep eq 2;
    if d eq 3 then
        assert not LocalFakeClass(f,L!(z0-e),GF(83)!(-z0),83);
    end if;
end for;
print "Type 7 squareclasses, disk parameters and two-cover classes verified.";
