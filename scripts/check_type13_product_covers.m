// Exact rank upper bounds accompanying the Type13 two-coset reduction.
SetColumns(0);SetSeed(1);Q:=Rationals();P<c>:=PolynomialRing(Q);
function UnconditionalBound(O)
 b:=MinkowskiBound(O);print "CLASS_PROOF",Degree(NumberField(O)),Discriminant(O),b;return b;
end function;
search:=map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>;
proof:=map<PowerStructure(RngOrd)->Integers()|O:->UnconditionalBound(O)>;
SetClassGroupBoundMaps(search,proof);
N:=-240*c^3+1600*c^2-260*c+475;f:=6*c*(25-4*c)*N;
assert Degree(f) eq 5 and GCD(f,Derivative(f)) eq 1;
H:=HyperellipticCurve(f);J:=Jacobian(H);lo,hi:=RankBounds(J);assert hi eq 2;
T,tm:=TwoTorsionSubgroup(J);assert #T eq 4;
DP:=J![H![1/4,135,1],H![1,0,0]];DQ:=J![H![19/4,2565/2,1],H![1,0,0]];
DT0:=J![H![0,0,1],H![1,0,0]];DTA:=J![H![25/4,0,1],H![1,0,0]];
assert Order(DT0) eq 2 and Order(DTA) eq 2;
print "PRODUCT_JACOBIAN_RANK_UPPER",hi;
// The companion Python checker verifies the Kummer independence and cosets.
K<th>:=NumberField(c^3-Q!20/3*c^2+Q!13/12*c-Q!95/48);R<x>:=PolynomialRing(K);
assert Evaluate(N,th) eq 0 and th ne 0 and th ne K!25/4;
alpha:=K!25/4;delta:=alpha-th;d:=-6*delta;
E:=EllipticCurve([K|0,-d*(alpha+th),0,d^2*alpha*th,0]);
Ep,phi:=IsogenyFromKernel(E,x-d*th);psi:=DualIsogeny(phi);
S,sm:=TwoSelmerGroup(Ep);T2,t2m:=TwoTorsionSubgroup(Ep);
assert AbelianInvariants(S) eq [2,2] and #T2 eq 2;
assert Valuation(#S,2)-Valuation(#T2,2) eq 1;
Ae:=d*(2*th-alpha);Be:=d^2*th*(th-alpha);
Es:=EllipticCurve([K|0,-2*Ae,0,Ae^2-4*Be,0]);
assert IsIsomorphic(Es,Ep);
assert Ae^2-4*Be eq d^2*alpha^2;
print "ALPHA_ISOGENOUS_RANK_UPPER",1;
print "PASS: exact rank upper bounds for the product Jacobian and alpha genus-one cover.";
