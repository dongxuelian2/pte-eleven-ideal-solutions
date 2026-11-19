// Proposition 41. All coefficients in the sieve are integral away from 2 and 3.
SetColumns(0);
SetSeed(2);
function Cubic(a,c,z)
    return 256*a^3-1280*a^2*c-64*a^2*z+2560*a*c^2-2560*a*c*z-2976*a*z^2
           -960*c^3+3920*c^2*z+780*c*z^2+3159*z^3;
end function;

Q := Rationals();
P<a,c,z> := ProjectiveSpace(Q,2);
C := Curve(P,Cubic(a,c,z));
assert IsNonsingular(C);
E,phi := EllipticCurve(C,C![4,1,4]);
M := MinimalModel(E);
assert aInvariants(M) eq [Q|0,1,0,-13533,-323937];
lo,hi := RankBounds(M);
assert lo eq 2 and hi eq 2 and #TorsionSubgroup(M) eq 1;
P0 := phi(C![-14,1,4]); Q0 := phi(C![36,9,4]);
for l in [2,3], i,j in [0..l-1] do
    if i+j gt 0 then assert #DivisionPoints(i*P0+j*Q0,l) eq 0; end if;
end for;

function Sieve(classes,m,primes)
    for p in primes do
        F := GF(p);
        S<x> := PolynomialRing(F);
        P<a,c,z> := ProjectiveSpace(F,2);
        C := Curve(P,Cubic(a,c,z));
        assert IsNonsingular(C);
        E,phi := EllipticCurve(C,C![4,1,4]);
        G,points := AbelianGroup(E);
        W,rho := quo<G|[m*g:g in Generators(G)]>;
        tab := AssociativeArray();
        for g in G do tab[points(g)] := rho(g); end for;
        A := tab[phi(C![-14,1,4])]; B := tab[phi(C![36,9,4])];
        good := {W|};
        for pt in Points(C) do
            if pt[3] ne 0 then
                a := pt[1]/pt[3]; c := pt[2]/pt[3]; b := (10*a-4*c+15)/6;
                if not &and[IsSquare(t):t in [a,b,c]] then continue; end if;
                s := (16*a+20*c+81)/12;
                t := -(128*a^2-640*a*c+48*a+80*c^2-1380*c-1701)/144;
                roots := Roots(x^2-s*x+t);
                if &+[Integers()|r[2]:r in roots] ne 2 or
                   not &and[IsSquare(r[1]):r in roots] then continue; end if;
            end if;
            Include(~good,tab[phi(pt)]);
        end for;
        classes := {k:k in classes|(k mod m)*A+(k div m)*B in good};
    end for;
    return classes;
end function;

S9 := Sieve({0..80},9,[7,11,13,17,19,23,29,31,37,41,61,73,139,149,181,193,337]);
assert S9 eq {0,9,17,38,60,65,66};
S27 := Sieve({i+27*j:i,j in [0..26]|(i mod 9)+9*(j mod 9) in S9},27,
             [41,181,353,359,389,491,503,563,577,673,761,821,859,937]);
assert S27 eq {0,27,53,677,678};
m := 216;
S216 := Sieve({i+m*j:i,j in [0..m-1]|(i mod 27)+27*(j mod 27) in S27},m,
              [7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,101,103,107,109,113]);
assert S216 eq {0,216,431,46226,46227};
for row in [[0,0,1,1], [0,1,9,9], [-1,1,4,49], [2,-2,4,1], [3,-2,1,25]] do
    i,j,a,c4 := Explode(row);
    assert (i mod m)+m*(j mod m) in S216;
    assert i*P0+j*Q0 eq phi(C![a,c4/4,1]);
end for;
print "Type 7 has rank two and five classes modulo 216.";
