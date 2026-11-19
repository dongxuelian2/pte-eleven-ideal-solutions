// Global groups and elliptic Chabauty for Theorem 30, Type 3.
SetColumns(0);
SetSeed(1);
Q := Rationals();
R<t> := PolynomialRing(Q);
SetClassGroupBoundMaps(
    map<PowerStructure(RngOrd) -> Integers() | O :-> Min(1000, MinkowskiBound(O))>,
    map<PowerStructure(RngOrd) -> Integers() | O :-> MinkowskiBound(O)>);

function RationalCover(scale)
    return EllipticCurve([Q | 0, -57*scale, 0, 1656*scale^2, 0]);
end function;

E6 := RationalCover(6);
lo, hi := RankBounds(E6);
assert lo eq 0 and hi eq 0 and #TorsionSubgroup(E6) eq 2;
assert Order(E6![0,0,1]) eq 2;

E10 := RationalCover(10);
lo, hi := RankBounds(E10);
assert lo eq 1 and hi eq 1 and #TorsionSubgroup(E10) eq 2;
P0 := E10![160,4000,1];
assert Order(P0) eq 0 and Order(E10![0,0,1]) eq 2;
assert IsSquare(Q!(1656*10^2)/46) and IsSquare(P0[1]/10);
assert not IsSquare(Q!10) and not IsSquare(Q!46) and not IsSquare(Q!(10*46));

K<w> := NumberField(t^2+15);
r := (57+15*w)/2;
assert r^2-57*r+1656 eq 0;
q := t^2-131/8*t+3381;
scale := 64*(16-r)*3/5;
E := EllipticCurve([K | 0, scale*Evaluate(Derivative(q),r), 0,
                        scale^2*Evaluate(q,r), 0]);
assert aInvariants(E) eq [K | 0,45300-18900*w,0,488160000*w-4186080000,0];
X0 := scale*(16-r);
P := E![X0,360*X0,1];
assert X0 eq 7200*w-26400 and Order(P) eq 0;
sel := TwoSelmerGroup(E);
tors, incl := TorsionSubgroup(E);
assert #sel eq 4 and #tors eq 2 and #TwoTorsionSubgroup(E) eq 2;
T := incl(tors.1);
assert T eq E![0,0,1];

// At (19,w-17), P and P+T both have nonsquare x, hence neither is twice a point.
F := GF(19);
ww := F!17;
assert ww^2+15 eq 0;
a2 := 45300-18900*ww;
a4 := 488160000*ww-4186080000;
xp := 7200*ww-26400;
assert a4*(a2^2-4*a4) ne 0 and IsSquare(a4);
assert not IsSquare(xp) and not IsSquare(a4/xp);

G := AbelianGroup([2,0]);
mp := map<G -> E | v :-> Eltseq(v)[1]*T + Eltseq(v)[2]*P>;
P1 := ProjectiveSpace(Q,1);
phi := Extend(map<E -> P1 | [E.1+scale*r*E.3, scale*E.3]>);
V, bad_index := Chabauty(mp,phi : IndexBound:=2);
assert #V eq 3 and IsPowerOf(bad_index,2);
assert {phi(mp(v)) : v in V} eq {P1 | [1,0], [16,1]};
print "Type 3: rank 0 in the first class; the second has t = infinity or 16.";
print "Chabauty returned",#V,"points; index obstruction",bad_index;
