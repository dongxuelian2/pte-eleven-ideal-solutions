// Reduction of the rank-four Type 1 product group modulo twice the group.
SetColumns(0);
SetSeed(1);
Q:=Rationals();
R<x>:=PolynomialRing(Q);
a:=16384;
f:=a*x^3+38672*x^2+33144*x+11025;
E:=EllipticCurve([Q|0,38672,0,a*33144,a^2*11025]);
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd)->Integers() | O:->Min(1000,MinkowskiBound(O))>,
    map<PowerStructure(RngOrd)->Integers() | O:->MinkowskiBound(O)>);
basis,proved,full:=Generators(E);
assert proved and full and #basis eq 4 and #TorsionSubgroup(E) eq 1;
reps:=[&+[E|basis[j]:j in [1..4]|IsOdd(mask div 2^(j-1))]:mask in [0..15]];
anchors:=[];
for b in [0,1] do
    ok,y:=IsSquare(Evaluate(f,b));
    assert ok;
    P:=E![a*b,a*y,1];
    indices:=[i:i in [1..16]|#DivisionPoints(P-reps[i],2) gt 0];
    assert #indices eq 1;
    Append(~anchors,indices[1]);
end for;
assert #Seqset(anchors) eq 2;
function SquareRootsOnly(g)
    roots:=Roots(g);
    return &+[Integers()|pair[2]:pair in roots] eq 3 and
           &and[IsSquare(pair[1]):pair in roots];
end function;
surviving:={1..16};
sizes:=[];
for p in [7,67,73,107,257] do
    assert Valuation(Discriminant(E),p) eq 0;
    F:=GF(p);
    S<z>:=PolynomialRing(F);
    EF:=ChangeRing(E,F);
    G,points:=AbelianGroup(EF);
    H,reduce:=quo<G|[2*G.i:i in [1..Ngens(G)]]>;
    parity:=AssociativeArray();
    for g in G do parity[points(g)]:=reduce(g);end for;
    good:={H|};
    if IsSquare(F!5) and SquareRootsOnly(64*z^3-1152*z^2+7232*z-16384) then
        Include(~good,H!0);
    end if;
    for b in F do
        if not IsSquare(b) or not IsSquare(5*b+4) then continue;end if;
        g:=64*z^3-(1152*b+1328)*z^2+(7232*b^2+14080*b+7564)*z
           -(16384*b^3+38672*b^2+33144*b+11025);
        if not SquareRootsOnly(g) then continue;end if;
        ok,y:=IsSquare(Evaluate(ChangeRing(f,F),b));
        assert ok;
        Include(~good,parity[EF![a*b,a*y,1]]);
    end for;
    surviving:={i:i in surviving|parity[EF!reps[i]] in good};
    Append(~sizes,#surviving);
end for;
assert sizes eq [8,6,5,4,2] and surviving eq Seqset(anchors);
print "Every Type 1 lift lies in the class of the point over b=0 or b=1.";

// The abscissa character at each split prime of the cubic field.
p:=73;
F:=GF(p);
S<z>:=PolynomialRing(F);
rs:=Roots(ChangeRing(f,F));
assert #rs eq 3 and not IsSquare(F!5);
for b0 in [0,1] do
    candidates:={};
    for b in F do
        if not IsSquare(b) or not IsSquare(5*b+4) then continue;end if;
        if not &and[IsSquare((b-u[1])*(b0-u[1])):u in rs] then continue;end if;
        g:=64*z^3-(1152*b+1328)*z^2+(7232*b^2+14080*b+7564)*z
           -(16384*b^3+38672*b^2+33144*b+11025);
        if SquareRootsOnly(g) then Include(~candidates,b);end if;
    end for;
    assert candidates eq {F!b0};
end for;
print "Each remaining class has one center residue at 73.";
