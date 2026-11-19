// Local Kummer classes and primitive halves on the sextic cover.
SetColumns(0);
Q := Rationals(); A<x> := PolynomialRing(Q);
K<h> := NumberField(x^6-96*x^4+144*x^2-54);
n := 4/27*h^5-127/9*h^3+32/3*h;
r := -5/81*h^4+53/9*h^2-41/9;
assert n^2 eq 2/3 and 9*r^3-5*r^2-3*r+1 eq 0;
ss := [K|-73/3*h^5-187/9*h^4+6950/3*h^3+1978*h^2-1664*h-1405,
    -8/9*h^5+10/9*h^4+254/3*h^3-317/3*h^2-66*h+79,
    8/9*h^5+17/9*h^4-254/3*h^3-180*h^2+64*h+140,
    h^5+5/9*h^4-286/3*h^3-53*h^2+80*h+47];
kap := 314/3*h^5+70*h^4-29924/3*h^3-6677*h^2+8049*h+5583;
U := 1222/15+4438/45*h-3581/30*h^2-13061/90*h^3+169/135*h^4+247/162*h^5;
V := 109023/194+146688/485*h-152419/194*h^2-466127/970*h^3+16011/1940*h^4+9801/1940*h^5;
dec := Decomposition(MaximalOrder(K),2);
assert #dec eq 1 and dec[1][2] eq 6;
L,em := Completion(K,dec[1][1]:Precision:=600);
n := em(n); r := em(r); s := em(ss[2]); kap := em(kap); U := em(U); V := em(V);
Pminus := func<v|9*v^3+5*v^2-3*v-1>;
k := (n-r)*Pminus(n); m := k*Pminus(r); j := -m/(2*r);
a0 := k*(27*r^2+10*r-3);
a := a0+3*j; b := 3*j^2+2*j*a0+k*m*(27*r+5); D := a^2-4*b;
E := EllipticCurve([L|0,a,0,b,0]); T := E![0,0,1];
P0 := E![m/(n-r)-j,m*k/(n-r)^2,1];

function Sq(t)
    assert AbsolutePrecision(t) gt Valuation(t)+12;
    return IsSquare(t);
end function;

assert not Sq(b) and not Sq(D);
assert not Sq((-n-r)*Pminus(-n)/k);
assert {i:i in [0..15] | Sq(&*[L|em(ss[j]):j in [1..4] | IsOdd(i div 2^(j-1))])} eq {0,2};
assert Sq(s);
for sign in [-1,1] do
    assert not Sq(kap*(2*s*U+a*V+sign*2*Sqrt(s)));
end for;

function Half(P)
    X := P[1]; Y := P[2];
    if Sq(X) then
        rr := Sqrt(X);
        for root in [rr,-rr] do
            v := -Y/root; c := 2*X+a+2*v;
            if Sq(c) then
                t := Sqrt(c); delta := root*t+v;
                H := E![X+delta,-Y-(root+t)*delta,1];
                Doubled := 2*H;
                for coord in [1,2] do
                    err := Doubled[coord]-P[coord];
                    assert AbsolutePrecision(err) gt 40 and Valuation(err) gt 40;
                end for;
                return true,H;
            end if;
        end for;
    end if;
    return false,P;
end function;

for z in [5,9] do
    w := Sqrt(k*(n*z-r)*Pminus(n*z));
    if Valuation(w-k) lt Valuation(w+k) then w := -w; end if;
    P := E![m/(n*z-r)-j,m*w/(n*z-r)^2,1]-P0;
    assert Valuation(-P[1]/P[2]) eq (z eq 5 select 9 else 15);
    for i in [1..(z eq 5 select 1 else 2)] do
        yes,H := Half(P); if not yes then yes,H := Half(P+T); end if;
        assert yes; P := H;
    end for;
    if not Sq(P[1]) then P +:= T; end if;
    assert Sq(P[1]) and not Half(P) and not Half(P+T);
    u := Sqrt(P[1]/s); v := P[2]/(s*u);
    l := (2*s*U+a*V)*u^2+(a*U+2*b/s*V)-2*v;
    assert not Sq(l*kap) and not Sq(l*kap*D);
end for;
R<X> := PolynomialRing(L);
g := (n+4*n*X-r)*Pminus(n+4*n*X)/k;
assert [Valuation(Coefficient(g,i)):i in [1..4]] eq [17,32,47,60];
B := [L|1];
for i in [1..11] do
    Append(~B,-&+[L|(2*i-j)*Coefficient(g,j)*B[i-j+1]:j in [1..Min(i,4)]]/(2*i));
end for;
assert [Valuation(2*n*B[i]/i):i in [1..12]] eq [9,14,25,24,35,40,51,44,55,60,71,70];
print "Dyadic injectivity test, primitive halves and analytic coefficients verified.";
