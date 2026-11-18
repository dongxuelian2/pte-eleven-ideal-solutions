// Type 13 alpha. Certify the two local Kummer classes used in Section 5.50.
SetColumns(0);
SetSeed(1);
Q := Rationals();
R<u> := PolynomialRing(Q);
search := map<PowerStructure(RngOrd) -> Integers() | O :-> Min(1000, MinkowskiBound(O))>;
proof := map<PowerStructure(RngOrd) -> Integers() | O :-> MinkowskiBound(O)>;
SetClassGroupBoundMaps(search, proof);

K<uu> := NumberField(u^3+u^2+315*u+1983);
O := MaximalOrder(K);
theta := (uu^2-3*uu+342)/60;
alpha := K!25/4;
delta := alpha-theta;
d := -6*delta;
a := 2*d*(alpha-2*theta);
b := d^2*alpha^2;
a2 := 60^2*a;
a4 := 60^4*b;
E := EllipticCurve([K|0,a2,0,a4,0]);
assert forall{z : z in aInvariants(E) | z in O};
places := Decomposition(O,2);
assert #places eq 1 and places[1][2] eq 3;
pr := places[1][1];
pi := (uu+1)^2/2;
assert pi^3+314*pi^2+26632*pi-347778 eq 0;
assert Valuation(pi,pr) eq 1;
assert Valuation(pi^10*a,pr) eq 4 and Valuation(pi^20*b,pr) eq 2;
assert Norm(theta*(theta-alpha)) eq 35^2*19/2^10;

S, to_S := TwoSelmerGroup(E);
assert AbelianInvariants(S) eq [2,2];
A := Domain(to_S);
tau := A.1;
pol := Parent(Modulus(A)).1;
assert Modulus(A) eq pol^3+a2*pol^2+a4*pol;
loc_map := LocalTwoSelmerMap(A,pr);
images := [loc_map(g@@to_S) : g in Generators(S)];
global_image := sub<Codomain(loc_map) | images>;
assert #global_image eq #S;
print "PASS: global 2-Selmer localization is injective";

// Exact approximations in Q(pi), with every error bounded by Hensel's lemma.
z := 25677+35360*pi+43700*pi^2;
assert Valuation(z-1,pr) ge 6;
z_residual := Valuation(z^2-(1-6/delta),pr);
assert z_residual gt 2*Valuation(2*z,pr);

// The reference point has -delta*X/Y=1 on the unscaled model.
r := -150*delta*(1+z);
x := r^2;
y := -60*delta*x;
assert Valuation(r,pr) eq 0 and Valuation(y,pr) eq 0;
point_accuracy := z_residual-Valuation(2*z,pr)-3;
assert point_accuracy gt 2*Valuation(y,pr)+6;
assert loc_map(A!x-tau) eq Codomain(loc_map)!0;

product := -y/r;
radicand := 2*x+a2+2*product;
s := 30760+25828*pi+61728*pi^2;
assert Valuation(s,pr) ge 0;
residual := Min(point_accuracy,Valuation(s^2-radicand,pr));
assert residual gt 2*Valuation(2*s,pr);
half_accuracy := Min(point_accuracy,residual-Valuation(2*s,pr));
correction := r*s+product;
x_half := x+correction;
y_half := -y-(r+s)*correction;
assert Valuation(x_half,pr) ge 0 and Valuation(y_half,pr) ge 0;
assert half_accuracy gt 2*Valuation(y_half,pr)+6;
assert loc_map(A!x_half-tau) notin global_image;
print "PASS: reference locally divisible by two, half outside global 2-Selmer image";
print "Certified reference and half precision", point_accuracy, half_accuracy;
print "PASS: Type 13 alpha local certificate";
