SetColumns(0);SetSeed(1);Q:=Rationals();P<x>:=PolynomialRing(Q);
cg_search:=map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>;
cg_proof:=map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>;SetClassGroupBoundMaps(cg_search,cg_proof);
qs:=[3*(240*x^2+40*x+167),-80*x^2+220*x+31,3*(-720*x^2+2152*x+275)];f:=&*qs;C:=HyperellipticCurve(f);H,mp:=TwoCoverDescent(C);A<th>:=Domain(mp);assert #H eq 3;surv:=[];
for h in H do
 delta:=h@@mp;good:=true;
 for q in qs do K<r>:=NumberField(q/LeadingCoefficient(q));proj:=hom<A->K|r>;nn:=Norm(proj(delta))*LeadingCoefficient(q); if not IsSquare(nn) then good:=false;end if;end for;
 if good then Append(~surv,h);end if;
end for;
assert #surv eq 1 and Seqset(surv) eq {mp(A!(1/4)-th)};

SetSeed(1);
K<ss>:=NumberField(x^2-x-1);r:=(14+27*ss)/20;gamma:=1/4-r;c0:=Q!1/4;
a:=720;b:=120;d:=501;L:=gamma*a;
E0:=EllipticCurve([K|0,gamma*(b-r*a),0,L*gamma*(d-r*b),-L^2*gamma*r*d]);E:=MinimalModel(E0);iso:=Isomorphism(E0,E);
assert RootNumber(E) eq -1;lo,hi:=RankBounds(E);assert lo eq 1 and hi eq 1;
pt:=iso(E0![L*c0,L*gamma*24,1]);assert Order(pt) eq 0;
free:=Saturation([pt],2:TorsionFree:=true);assert #free eq 1;
T,tm:=TorsionSubgroup(E);assert #T eq 2;ng:=Ngens(T);G:=AbelianGroup(AbelianInvariants(T) cat [0]);
m:=map<G->E|z:->&+[Eltseq(z)[i]*tm(T.i):i in [1..ng]]+Eltseq(z)[ng+1]*free[1]>;
P1:=ProjectiveSpace(Q,1);phi0:=map<E0->P1|[E0.1,L*E0.3]>;phi:=Expand(Inverse(iso)*phi0);phi:=Extend(phi);
V,bd:=Chabauty(m,phi:IndexBound:=2);assert #V eq 3 and IsPowerOf(bd,2);assert {phi(m(v)):v in V} eq {P1|[1/4,1],[1,0]};

print "PASS: individual norms select one fake class; normalized elliptic images are 1/4 or infinity.";
