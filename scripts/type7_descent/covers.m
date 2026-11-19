// Rank and point certificates for the cubic-field covers.
SetColumns(0);
SetSeed(2);
Q := Rationals(); A<x> := PolynomialRing(Q);
SetClassGroupBoundMaps(map<PowerStructure(RngOrd)->Integers()|O:->Min(1000,MinkowskiBound(O))>,
                       map<PowerStructure(RngOrd)->Integers()|O:->MinkowskiBound(O)>);
K<h> := NumberField(x^3-8*x-2); e := (-5*h^2+4*h+53)/27;
assert 27*e^3-79*e^2+57*e-9 eq 0;
g1 := e-79/27; g0 := e^2-79/27*e+57/27;
for row in [[15,0,1,0],[2,2,1,3],[3,1,3,3],[3,3,1,2],[10,2,5,2]] do
    d,r,s,expected := Explode(row); t0 := Q!(r/s);
    delta := r eq 0 select K!135 else t0-e; k := 135*d*delta;
    E := EllipticCurve([K|0,k*g1,0,k^2*g0,0]);
    lo,hi := RankBounds(E); assert lo eq expected and hi eq expected;
    if r eq 0 then assert #TorsionSubgroup(E) eq 2; end if;
end for;
L<r> := NumberField(x^3-5/9*x^2-1/3*x+1/9);
pminus := func<v|9*v^3+5*v^2-3*v-1>;
z0 := Q!1/3; k := (z0-r)*pminus(z0); m := k*pminus(r);
E := EllipticCurve([L|0,k*(27*r^2+10*r-3),0,k*m*(27*r+5),9*k*m^2]);
lo,hi := RankBounds(E); assert lo eq 1 and hi eq 1;
P := E![m/(z0-r),m*k/(z0-r)^2,1]; assert Order(P) eq 0;
print "Type 7 cubic-field ranks verified.";
