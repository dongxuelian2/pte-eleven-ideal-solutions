SetColumns(0);SetSeed(1);Q:=Rationals();P<x>:=PolynomialRing(Q);
cg_search:=map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>;
cg_proof:=map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>;SetClassGroupBoundMaps(cg_search,cg_proof);
K<r>:=NumberField(x^2-49/2*x-283/80);
ss := (20*r-245)/64;
assert ss^2 eq 15 and Discriminant(MaximalOrder(K)) eq 60;
assert MinkowskiBound(K) eq 3 and ClassNumber(K) eq 2;
assert Norm(3+ss) eq -6 and Norm(4+ss) eq 1;
assert (3+ss)^2 eq 6*(4+ss);
assert (49/4-r)/(6*(1/4-r)) eq (2/(3+ss))^2;
rows:=[<Q!1/4,3*(-80*x^2+808*x+2875),96>,<Q!49/4,3*(2*x-23)*(4*x-31),9>];
for row in rows do
 c0,q,s0:=Explode(row);gamma:=c0-r;a,b,d:=Explode(Reverse(Eltseq(q)));
 L:=gamma*a;E:=EllipticCurve([K|0,gamma*(b-r*a),0,L*gamma*(d-r*b),-L^2*gamma*r*d]);
 pt:=E![L*c0,L*gamma*s0,1];print "C0",c0,"POINT_ORDER",Order(pt);assert Order(pt) eq 0;
 lo,hi:=RankBounds(E);assert lo eq 1 and hi eq 1;free:=Saturation([pt],2:TorsionFree:=true);assert #free eq 1;
 T,tm:=TorsionSubgroup(E);print "TORSION",T;ng:=Ngens(T);G:=AbelianGroup(AbelianInvariants(T) cat [0]);
 m:=map<G->E|z:->&+[Eltseq(z)[i]*tm(T.i):i in [1..ng]]+Eltseq(z)[ng+1]*free[1]>;
 P1:=ProjectiveSpace(Q,1);phi:=map<E->P1|[E.1,L*E.3]>;phi:=Extend(phi);
 V,bd:=Chabauty(m,phi:IndexBound:=2);assert IsPowerOf(bd,2);
images := {phi(m(v)):v in V};
expected := c0 eq 1/4 select {P1|[1/4,1],[1,0]}
                         else {P1|[31/4,1],[23/2,1],[49/4,1],[1,0]};
assert images eq expected;
assert #V eq (c0 eq 1/4 select 3 else 5);
print "CERTIFIED_C_VALUES",images,"INDEX_BOUND",bd;
end for;

assert not IsSquare(Q!31/4) and not IsSquare(Q!23/2);
print "PASS: unit-normalized cases exhaust type15 and leave only nonsquares or collisions.";
