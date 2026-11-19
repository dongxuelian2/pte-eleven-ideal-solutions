// Complete all-height exclusion of symmetric shift-chain type 14.
// Type (m,aa,bb)=(3,(2),(1,1,1,2)); see PROOF section 5.3.
// Requires Magma. No GRH or guessed Mordell-Weil index is used.
SetColumns(0);
SetSeed(1);
Q := Rationals();
P<b> := PolynomialRing(Q);

// Every class-group proof uses the actual unconditional Minkowski bound.
// The callback prints all orders for which the descent requests a proof.
function UnconditionalClassBound(O)
    bound := MinkowskiBound(O);
    print "CLASS_PROOF", Degree(NumberField(O)), Discriminant(O), bound;
    return bound;
end function;
cg_search := map< PowerStructure(RngOrd) -> Integers() |
    O :-> Min(1000,MinkowskiBound(O)) >;
cg_proof := map< PowerStructure(RngOrd) -> Integers() |
    O :-> UnconditionalClassBound(O) >;
SetClassGroupBoundMaps(cg_search,cg_proof);

h := 900*b^5-7100*b^4-21675*b^3+71658*b^2+3321*b-93312;
f := -3*h;
assert Degree(f) eq 5 and GCD(f,Derivative(f)) eq 1;
C := HyperellipticCurve(f);
J := Jacobian(C);
lower,upper := RankBounds(J);
assert upper eq 1;
print "CERTIFIED_RANK_UPPER", upper;

pinf := C![1,0,0];
p9 := C![9,3240,1];
D := J![p9,pinf];
assert Order(D) eq 0; // Magma convention: zero means infinite order.

// The modern Chabauty intrinsic includes saturation of <D> internally.
// At rank one this returns ALL rational points; D need not be a generator.
// Official specification: magma/handbook/text/1620, Chabauty(P : ptC).
points := Chabauty(D : ptC := pinf);
expected := { pinf, p9, C![9,-3240,1] };
assert points eq expected;
assert { Q | p[1]/p[3] : p in points | p[3] ne 0 } eq { Q | 9 };
print "COMPLETE_RATIONAL_POINTS", points;

// Exact singleton polynomial and discriminant identity, before specialization.
S<X> := PolynomialRing(P);
e1 := (16*b+129)/12;
e2 := -(128*b^2-1440*b-4779)/144;
e3 := (256*b^3-3072*b^2+5616*b+18225)/576;
K := X^3-e1*X^2+e2*X-e3;
assert Discriminant(K) eq -4*b*h/243;
K9 := Polynomial([ Evaluate(Coefficient(K,i),9) : i in [0..3] ]);
assert K9 eq (4*P.1-1)*(4*P.1-9)*(4*P.1-81)/64;
assert Evaluate((5*b+3)/3,9) eq 16;

// Signs of centers only permute the reflected chains. At b=9 they give
// u=+/-4, v=+/-3 and singleton absolute centers 1/2,3/2,9/2.
A := { Q | (2*j-5)/2 : j in [0..5] }
     join { Q | sg*4+d : sg in [-1,1], d in [-1/2,1/2] };
B := { Q | sg*3+d : sg in [-1,1], d in [-1/2,1/2] }
     join { Q | sg*r : sg in [-1,1], r in [1/2,3/2,9/2] };
assert #A eq 10 and A eq B;
assert A eq { Q | (2*j-9)/2 : j in [0..9] };
print "ONLY_FINITE_BRANCH_IS_COLLIDING", A;
print "PASS: type 14 excluded at every rational height.";
