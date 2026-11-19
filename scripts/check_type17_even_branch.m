// The norm bound for the degree-twelve field is proved in section 5.60.
// Run the explicit-formula and factor-base checks before this script.
SetColumns(0);
SetSeed(1);
Q := Rationals();
R<x> := PolynomialRing(Q);
a := 11025;
g := x^3 - 4131*x^2 + 472*a*x - 16*a^2;
A := x^4 - 2*472*a*x^2 + 128*a^2*x + a^2*(472^2-4*4131*16);
B := A + 20*a*g;

function FullBound(O)
    return MinkowskiBound(O);
end function;
search := map<PowerStructure(RngOrd) -> Integers() | O :-> Min(1000,MinkowskiBound(O))>;
SetClassGroupBoundMaps(search,
    map<PowerStructure(RngOrd) -> Integers() | O :-> FullBound(O)>);

// In odd degree the norm selects a representative modulo scalars and squares.
for quartic in [A,B] do
    assert IsIrreducible(g) and IsIrreducible(quartic);
    C := HyperellipticCurve(g*quartic);
    covers,mu := TwoCoverDescent(C);
    alg := Domain(mu);
    cubicField<r> := NumberField(g);
    quarticField<t> := NumberField(quartic);
    assert not IsSquare(a/16-r);
    toCubic := hom<alg -> cubicField | r>;
    toQuartic := hom<alg -> quarticField | t>;
    assert mu(alg!1) in covers;
    assert #covers eq (quartic eq A select 2 else 1);
    for element in covers do
        v := element@@mu;
        v /:= Norm(toCubic(v))*Norm(toQuartic(v));
        assert IsSquare(toCubic(v)) eq (element eq mu(alg!1));
        assert IsSquare(toQuartic(v)) eq (element eq mu(alg!1));
    end for;
end for;

K<w> := NumberField(x^4-194*x^2-96*x+8302);
theta := 3*(2*w^3+4*w^2-240*w-532);
P<U> := PolynomialRing(K);
h := ExactQuotient(Evaluate(A,U),U-theta);
assert IsMonic(h) and IsIrreducible(h);
L := NumberField(x^12-472*x^10+75282*x^8-4450612*x^6
                 +64570537*x^4-277288596*x^2+234027000);
M := AbsoluteField(ext<K|h>);
assert IsIsomorphic(M,L);
disc := 58160039859394517078535559761100800000;
OL := MaximalOrder(L);
assert Discriminant(OL) eq disc and IsTotallyReal(L);
function CertifiedBound(O, reference)
    if AbsoluteDegree(NumberField(O)) eq 12 and Discriminant(O) eq disc then
        assert IsIsomorphic(NumberField(O),reference);
        return 1000;
    end if;
    return MinkowskiBound(O);
end function;
SetClassGroupBoundMaps(search,
    map<PowerStructure(RngOrd) -> Integers() | O :-> CertifiedBound(O,L)>);
assert AbelianInvariants(ClassGroup(OL)) eq [2];
E := EllipticCurve([K|0,Coefficient(h,2),0,Coefficient(h,1),Coefficient(h,0)]);
S,mu := TwoSelmerGroup(E);
alg := Domain(mu);
assert #S eq 4 and Evaluate(h,alg.1) eq 0;
places := Decomposition(MaximalOrder(K),2);
assert #places eq 1 and places[1][2] eq 4;
p := places[1][1];
loc := LocalTwoSelmerMap(alg,p);
D,emb := Completion(K,p : Precision:=300);
hd := Polynomial([D|emb(Coefficient(h,i)): i in [0..3]]);
E2 := EllipticCurve([D|0,Coefficient(hd,2),0,Coefficient(hd,1),Coefficient(hd,0)]);
roots := Roots(hd);
assert #roots eq 1 and Valuation(emb(w)) eq 1 and Valuation(D.1) eq 1;

// A center in K approximating the unique local two-torsion abscissa.
root := roots[1][1];
center := K!0;
for i in [Valuation(root)..32] do
    if Valuation(root-emb(center)) eq i then center +:= w^i; end if;
end for;
xs := [center+1,center+w^2,w^4*(1+w^2+w^3),
       w^6*(1+w+w^2+w^3),center+w^14*(1+w^2)];
V := VectorSpace(GF(2),Ngens(Codomain(loc)));
rows := Matrix(GF(2),[Eltseq(loc(alg!z-alg.1)): z in xs]);
assert Rank(rows) eq 5;
globalCoordinates := Matrix(GF(2),[Eltseq(Solution(rows,V!Eltseq(loc(s@@mu)))):
                                    s in OrderedGenerators(S)]);

// Basis of 2^11 times the maximal ideal: 2^12, 2^11*pi, 2^11*pi^2, 2^11*pi^3.
orders := [12,11,11,11];
logRows := [];
for z in xs do
    flag, y := IsSquare(Evaluate(hd,emb(z)));
    assert flag;
    T := (2^12*105)*(E2![emb(z),y,1]);
    assert T ne E2!0;
    formal := -T[1]/T[2];
    assert Valuation(formal) ge 44 and AbsolutePrecision(formal) ge 60;
    coeff := Eltseq(formal);
    assert #coeff eq 4;
    assert &and[Valuation(coeff[i]) ge orders[i] and
                AbsolutePrecision(coeff[i]) gt orders[i]: i in [1..4]];
    Append(~logRows,[Integers()!(coeff[i]/2^orders[i]) mod 2: i in [1..4]]);
end for;
logs := Matrix(GF(2),logRows);
assert Rank(logs) eq 4;
plane := globalCoordinates*logs;
axis := Matrix(GF(2),1,4,[1,0,0,0]);
assert Rank(plane) eq 2 and Rank(VerticalJoin(plane,axis)) eq 3;
print "Type 17: Selmer plane and rational tangent span dimension three.";
