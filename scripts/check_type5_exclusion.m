// Genuine type-5 configurations: scalar-square lemma in PROOF.md 5.42.
// This does not classify arbitrary points in an even-degree fake class.
SetColumns(0);
SetSeed(1);
cg_search:=map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>;
cg_proof:=map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>;
SetClassGroupBoundMaps(cg_search,cg_proof);
Q := Rationals();
P<x> := PolynomialRing(Q);
qs := [80*x^2-840*x+589, 80*x^2-40*x+4869,
       80*x^2-1064*x+5125];
discs := [];
for q in qs do
    K0 := NumberField(q/LeadingCoefficient(q));
    assert MinkowskiBound(K0) lt 100000;
    Append(~discs, Discriminant(MaximalOrder(K0)));
end for;
assert discs eq [505,-95,-31];
f := &*qs;
assert LeadingCoefficient(f) eq 512000 and not IsSquare(Q!512000);
C := HyperellipticCurve(f);
H, mp := TwoCoverDescent(C);
A<theta> := Domain(mp);
delta49 := A!(49/4)-theta;
delta81 := A!(81/4)-theta;
assert #H eq 2 and {mp(delta49),mp(delta81)} eq H;
assert mp(delta49) ne mp(delta81);

K<r> := NumberField(qs[1]/80);
R<X> := PolynomialRing(K);
g := (X-r)*Evaluate(qs[3]/80,X);
gamma := Evaluate(g,49/4);
B<thetaB> := quo<R|g>;
j := hom<A->B|thetaB>;
assert Norm(j(delta49)) eq gamma;
E := EllipticCurve([K|0,gamma*Coefficient(g,2),0,
                        gamma^2*Coefficient(g,1),gamma^3*Coefficient(g,0)]);
pt := E![gamma*49/4,gamma^2,1];
assert Order(pt) eq 0;
lo,hi := RankBounds(E);
assert lo eq 1 and hi eq 1;
free := Saturation([pt],2:TorsionFree:=true);
assert #free eq 1 and Order(free[1]) eq 0;
T,tm := TorsionSubgroup(E);
assert #T eq 2;
G := AbelianGroup([2,0]);
m := map<G->E|z:->Eltseq(z)[1]*tm(T.1)+Eltseq(z)[2]*free[1]>;
P1 := ProjectiveSpace(Q,1);
phi := map<E->P1|[E.1,gamma*E.3]>;
V,index_bound := Chabauty(m,phi:IndexBound:=2);
assert #V eq 3 and IsPowerOf(index_bound,2);
phi := Extend(phi);
images := {phi(m(v)):v in V};
assert images eq {P1|[49/4,1],[1,0]};
// At c=49/4, the two B paired-center squares are 4 and 0.
c := Q!49/4;
U := (12*c-83)/16;
VV := (4*c-49)*(4*c-9)/256;
assert U eq 4 and VV eq 0;

// Second normalized fake class, retaining the original c-square condition.
K<r> := NumberField(qs[2]/80);
R<X> := PolynomialRing(K);
c0 := Q!81/4;
gamma := c0-r;
g := (X-r)*Evaluate(qs[3]/80,X);
B<thetaB> := quo<R|g>;
j := hom<A->B|thetaB>;
assert Norm(j(delta81)) eq Evaluate(g,c0);
h := X*(X-r)*Evaluate(qs[3],X);
assert gamma*Evaluate(h,c0) eq (gamma*9/2*128)^2;
b := K!(-1064); d := K!5125;
L := -gamma*r*d;
E := EllipticCurve([K|0,gamma*(d-r*b),0,
                        L*gamma*(b-80*r),L^2*80*gamma]);
pt := E![L/c0,L*gamma*9/2*128/c0^2,1];
assert Order(pt) eq 0;
lo,hi := RankBounds(E);
assert lo eq 1 and hi eq 1;
free := Saturation([pt],2:TorsionFree:=true);
assert #free eq 1 and Order(free[1]) eq 0;
T,tm := TorsionSubgroup(E);
assert #T eq 2;
G := AbelianGroup([2,0]);
m := map<G->E|z:->Eltseq(z)[1]*tm(T.1)+Eltseq(z)[2]*free[1]>;
phi := map<E->P1|[L*E.3,E.1]>;
phi := Extend(phi);
V,index_bound := Chabauty(m,phi:IndexBound:=2);
assert #V eq 3 and IsPowerOf(index_bound,2);
assert {phi(m(v)):v in V} eq {P1|[81/4,1],[0,1]};
// c=0 gives negative S. At c=81/4, the A and B seeds coincide.
assert (20*(Q!0)-85)/16 lt 0;
c := Q!81/4;
assert (20*c-85)/16 eq 20 and (80*c^2-840*c+589)/256 eq 64;
assert (12*c-83)/16 eq 10 and (4*c-49)*(4*c-9)/256 eq 9;
print "PASS: both normalized fake classes are excluded for genuine type-5 inputs.";
print "All rational heights are covered; the original size-eleven PTE problem remains open.";
