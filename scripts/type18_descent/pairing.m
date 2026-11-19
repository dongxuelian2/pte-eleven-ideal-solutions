// Explicit everywhere locally soluble quartics on the elliptic curve marked by 25.
SetColumns(0);
Q := Rationals(); A<x> := PolynomialRing(Q);
K<h> := NumberField(x^4-x^3-991*x^2+13345*x-19455);
r := (-h^3-23*h^2+585*h+5002)/73;
R<u> := PolynomialRing(K);
q := ExactQuotient(u^4-165*u^3+8778*u^2-172810*u+1057221,u-r);
k := 25-r;
E := EllipticCurve([K|0,k*Coefficient(q,2),0,k^2*Coefficient(q,1),k^3*Coefficient(q,0)]);

f1 := 1/73*(326191125074*h^3 + 6862578283530*h^2 - 228375184609562*h + 359306674735384)*u^4
    + 1/219*(-68146848515935*h^3 - 1146954650973021*h^2 + 42423477924819832*h - 67064446139806428)*u^3
    + 1/73*(532419837520287*h^3 + 9948662909456332*h^2 - 349661419343905572*h + 551521283298925706)*u^2
    + 1/219*(-17476367657300689*h^3 - 311397084968974479*h^2 + 11197832892950970757*h - 17680324780047370650)*u
    + 1/219*(70099333572212947*h^3 + 1276741669345440462*h^2 - 45426339226399038625*h + 71690263015899809025);
f2 := 1/73*(21676496499*h^3 - 813389539308*h^2 + 8232971883903*h - 11541580840841)*u^4
    + 1/219*(2668272288424*h^3 - 100136656976148*h^2 + 1013669592704180*h - 1421062099339992)*u^3
    + 1/73*(13684809418338*h^3 - 513673043233358*h^2 + 5200629586198914*h - 7290949365165148)*u^2
    + 1/219*(280668222680764*h^3 - 10541234662342296*h^2 + 106777148298199184*h - 149707943110710504)*u
    + 1/219*(719625048806446*h^3 - 27039042667262238*h^2 + 273964615213522625*h - 384132484185360732);

for f in [f1,f2] do
    assert Degree(f) eq 4 and IsSquarefree(f);
    a,b,c,d,e := Explode([Coefficient(f,i):i in [4,3,2,1,0]]);
    I := 12*a*e-3*b*d+c^2;
    J := 72*a*c*e+9*b*c*d-27*a*d^2-27*b^2*e-2*c^3;
    assert IsIsomorphic(E,EllipticCurve([K|0,0,0,-27*I,-27*J]));
    assert HasPointsEverywhereLocally(f,2);
end for;
// These prime divisors speed the auxiliary conic factorization.
StoreFactor([865363,529735207,9479541337913,12512171853893,73954328282694579111817]);
SetSeed(3,20191985);
assert CasselsTatePairing(HyperellipticCurve(f1),HyperellipticCurve(f2)) eq 1;
print "The two locally soluble quartics have Cassels-Tate pairing one.";
