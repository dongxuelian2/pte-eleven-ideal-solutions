// Complete rational-point computation for the elliptic quotient in PROOF 5.45.
SetColumns(0);SetSeed(1);Q:=Rationals();P<x>:=PolynomialRing(Q);
cg_search:=map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>;
cg_proof:=map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>;
SetClassGroupBoundMaps(cg_search,cg_proof);
f:=-240*x^4+760*x^2+81;C:=HyperellipticCurve(f);
E,phi:=EllipticCurve(C,C![0,9,1]);
assert aInvariants(E) eq [Q|0,-1520/81,0,655360/6561,0];
lo,hi:=RankBounds(E);assert lo eq 0 and hi eq 0;
Emin:=MinimalModel(E);assert aInvariants(Emin) eq [Q|1,0,0,-28,272];
assert Discriminant(Emin) eq -31104000;
assert #ChangeRing(Emin,GF(7)) eq 10 and #ChangeRing(Emin,GF(11)) eq 10;
T,tm:=TorsionSubgroup(E);assert AbelianInvariants(T) eq [10];
inv:=Inverse(phi);seen:={};
for t in T do
 pt:=inv(tm(t));assert pt[3] ne 0;
 v:=pt[1]/pt[3];y:=pt[2]/pt[3]^2;
 assert y^2 eq Evaluate(f,v);
 Include(~seen,<v,y>);
 a:=(y+20*v^2+29)/8;
 if v eq 0 then assert not IsSquare(a);
 else assert v^2 in {Q!1/4,Q!9/4};end if;
end for;
expected:={<Q!0,Q!9>,<Q!0,Q!-9>};
for v in [Q!1/2,Q!-1/2] do for y in [Q!16,Q!-16] do Include(~expected,<v,y>);end for;end for;
for v in [Q!3/2,Q!-3/2] do for y in [Q!24,Q!-24] do Include(~expected,<v,y>);end for;end for;
assert seen eq expected and #seen eq 10;
print "TYPE2_ALL_RATIONAL_QUOTIENT_POINTS",Sort(Setseq(seen));
print "TYPE2_ALL_HEIGHT_EXCLUSION_OK";
