// Optional Magma cross-check of the exact proof in PROOF.md Section 5.38.
// This computes isogenies and curve isomorphisms, not ranks or Selmer groups.
SetColumns(0);
Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=3*x*(240*x^2+40*x+167)*(256*x^2+97*x+466);
g:=-6477*(128*x^2-233)*(240*x^2-167)*(13040*x^2+138176*x+2441);
C:=HyperellipticCurve(f); D:=HyperellipticCurve(g);
for pair in [<C,D>,<D,C>] do
    E,F:=Explode(pair);
    neighbors:=RichelotIsogenousSurfaces(E);
    assert #neighbors eq 1;
    assert Type(neighbors[1]) eq CrvHyp;
    assert IsIsomorphic(neighbors[1],F);
    assert #DoubleRichelotIsogenies(E) eq 0;
end for;
Js,Es,Ks:=TwoPowerIsogenies(Jacobian(C));
assert #Js eq 1 and #Es eq 0 and #Ks eq 0;
assert IsIsomorphic(Curve(Js[1]),D);
assert not IsIsomorphic(C,D);
print "PASS: two distinct Jacobian vertices, one edge, no double-Richelot edges or split surfaces.";
