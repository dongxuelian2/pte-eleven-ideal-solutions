// Type 6 on the smooth cubic of section 5.32.
SetColumns(0); SetSeed(1);
function Cubic(a,c,z)
    return a^3-5*a^2*c-9*a^2*z+10*a*c^2+15*a*c*z+24*a*z^2
           -10*c^3+10*c^2*z-20*c*z^2-16*z^3;
end function;
Q:=Rationals(); P<a,c,z>:=ProjectiveSpace(Q,2);
C:=Curve(P,Cubic(a,c,z));assert IsNonsingular(C);
E,phi:=EllipticCurve(C,C![1,0,1]);M:=MinimalModel(E);
assert aInvariants(M) eq [Q|1,1,1,-2913,80031];
lo,hi:=RankBounds(M);assert lo eq 2 and hi eq 2;
assert AbelianInvariants(TorsionSubgroup(M)) eq [2];
den:=3*a^2-8*a*c-7*a*z+6*c^2+2*c*z+4*z^2;
num:=16*a^3-48*a^2*c-29*a^2*z+48*a*c^2-8*a*c*z+9*a*z^2
     -16*c^3+38*c^2*z-14*c*z^2+4*z^3;
assert num-(5*a-c+5*z)*den eq Cubic(a,c,z);

quart:=a*c^3-5*a^2*c^2+10*a^3*c-10*a^4+(3*c^2-25*a*c+50*a^2)*z^2-40*z^4;
assert Cubic(4*z^2+a*c,a^2,z^2) eq a^2*quart;
assert IsNonsingular(Curve(P,quart));

P0:=phi(C![4,0,1]);Q0:=phi(C![4,1,1]);
T0:=phi(C![1,1,1])-P0;
assert T0 ne E!0 and 2*T0 eq E!0;
assert phi(C![9,1,1]) eq P0-Q0+T0;
assert phi(C![16,4,1]) eq -P0+2*Q0+T0;
for i,j,t in [0,1] do
    if i+j+t gt 0 then assert #DivisionPoints(i*P0+j*Q0+t*T0,2) eq 0;end if;
end for;
for i,j in [0..2] do
    if i+j gt 0 then assert #DivisionPoints(i*P0+j*Q0,3) eq 0;end if;
end for;

function Sieve(surv,m,primes)
    for p in primes do
        F:=GF(p); ambient<A,D,Z>:=ProjectiveSpace(F,2);
        C:=Curve(ambient,Cubic(A,D,Z));assert IsNonsingular(C);
        E,phi:=EllipticCurve(C,C![1,0,1]);G,points:=AbelianGroup(E);
        W,rho:=quo<G|[m*G.i:i in [1..Ngens(G)]]>;
        tab:=AssociativeArray();for u in G do tab[points(u)]:=rho(u);end for;
        g:=tab[phi(C![4,0,1])];h:=tab[phi(C![4,1,1])];
        tau:=tab[phi(C![1,1,1])]-g;
        good:={W|};S<x>:=PolynomialRing(F);
        for pt in Points(C) do
            if pt[3] ne 0 then
                a:=pt[1]/pt[3];c:=pt[2]/pt[3];b:=(5*a-c+5)/4;
                if not &and[IsSquare(u):u in [a,b,c]] then continue;end if;
                if a eq 4 and c eq 0 and not IsSquare(F!(40/3)) then continue;end if;
                sum:=(-3*a+15*c+22)/4;
                prod:=(8*a^2-40*a*c-67*a+80*c^2+95*c+149)/16;
                rs:=Roots(x^2-sum*x+prod);
                if &+[Integers()|u[2]:u in rs] ne 2 or
                   not &and[IsSquare(u[1]):u in rs] then continue;end if;
            end if;
            Include(~good,tab[phi(pt)]);
        end for;
        surv:={k:k in surv|(k mod m)*g+((k div m) mod m)*h+(k div m^2)*tau in good};
    end for;
    return surv;
end function;
s9:=Sieve({0..80},9,[7,11,13,17,41,61,73,307]);
assert s9 eq {1,26,47,73};
s27:=Sieve({i+27*j:i,j in [0..26]|(i mod 9)+9*(j mod 9) in s9},
           27,[197,271,397,647]);
assert s27 eq {1,80,703};
m:=216;
s216:=Sieve({i+m*j+m^2*t:i,j in [0..m-1],t in [0,1]|
             (i mod 27)+27*(j mod 27) in s27},m,[7,13,23,31,41,43,67,73]);
assert s216 eq {m^2+1,m^2+(m-1)+2*m,m^2+1+(m-1)*m};
print "Three of 93312 global classes remain. Their representatives are (1,1), (9,1) and (16,4).";
