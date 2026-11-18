// Type 13, Proposition 26. Run with Magma V2.30 or later.
// The determinant and both finite-field sieves use unconditional arithmetic.
SetColumns(0);
Q := Rationals();
R<T> := PolynomialRing(Q);
K<th> := NumberField(T^3-20*T^2/3+13*T/12-95/48);
O := MaximalOrder(K);
P<x> := PolynomialRing(K);

function N(c)
    return -240*c^3+1600*c^2-260*c+475;
end function;

function Quad(c,t)
    return t^2-(10*c+50)*t/3+5*c^2-30*c+19;
end function;

function Ht(c,t)
    return 2*Quad(c,t)*(2*t-(10*c+50)/3)-4*N(c)/9;
end function;

function Model(x,theta)
    return x^5*(908800*theta^2/27-18444800*theta/81+7356800/81)
        + x^4*(-108800*theta^2/9+2840320*theta/27-743680/27)
        + x^3*(2800*theta^2-133760*theta/9+21740/9)
        + x^2*(-980*theta^2/9+6520*theta/9-50/9)
        + x*(-25*theta^2/24-253*theta/24-247/96);
end function;

// The tangent ruling on Bruin's double cover, on a monic Mumford chart.
function Chart(c,t,r,theta,X)
    delta := 1/4-theta;
    q := Quad(c,t);
    s := q/r;
    qt := 2*t-(10*c+50)/3;
    qc := 10*c-30-10*t/3;
    dt := -(2*q*qc-4*t*(-720*c^2+3200*c-260)/9)/Ht(c,t);
    dr := (t+(c-theta)*dt)/(2*delta*r);
    ds := (qc+qt*dt-s*dr)/r;
    a := dt/delta-dr^2;
    b := dt^2-10*dt/3+5-dr*ds;
    k := -320*delta/3-ds^2;
    u := X^2+(2*delta*b/k)*X+delta^2*a/k;
    v := delta^2*((dr+X*ds/delta)*((c-theta)/(2*delta)-X*((10*c+50)/3-2*t)/delta)
        +(r+X*s/delta)*(-1/(2*delta)-X*(2*dt-10/3)/delta));
    return u, v mod u;
end function;

function Branch(c,base,tangentRoot)
    S := Parent(c);
    t := S!(tangentRoot^2/4);
    for i in [1..4] do
        t -:= (Quad(c,t)^2-4*t*N(c)/9)/Ht(c,t);
    end for;
    a := S!(tangentRoot/2);
    z := S!1;
    for i in [1..4] do
        a := (a+t/a)/2;
        z := (z+(c-base)/((1/4-base)*z))/2;
    end for;
    return t, a*z;
end function;

f := Model(x,th);
J := Jacobian(HyperellipticCurve(f));
roots := [1,3,7,9];
us := []; vs := []; ds := [];
for k in roots do
    u,v := Chart(K!1/4,K!(k^2)/4,K!k/2,th,x);
    assert (f-v^2) mod u eq 0;
    Append(~us,u); Append(~vs,v); Append(~ds,J![u,v]);
end for;
assert ds[1] eq ds[2]+ds[3]-ds[4];

// Recover the two pullbacks exactly in characteristic zero. The forms
// r dc/H_t and s dc/H_t form the anti-invariant basis by adjunction.
S<h> := PowerSeriesRing(K,8);
PS<W> := PolynomialRing(S);
controls := [];
for k in roots do
    c := S!(1/4)+h;
    t,r := Branch(c,th,K!k);
    u,v := Chart(c,t,r,th,W);
    a := Coefficient(u,1); b := Coefficient(u,0);
    A := Coefficient(v,0); B := Coefficient(v,1);
    den := 2*(A^2-a*A*B+b*B^2);
    l1 := (Derivative(a)*A-Derivative(b)*B)/den;
    l0 := (Derivative(b)*(A-a*B)+Derivative(a)*b*B)/den;
    s := Quad(c,t)/r;
    assert Valuation(-l1+2*s/((1/4-th)*Ht(c,t))) ge 3;
    assert Valuation(-l0+a*l1-2*r/Ht(c,t)) ge 3;
    Append(~controls,[Coefficient(r/Ht(c,t),0),Coefficient(s/Ht(c,t),0)]);
end for;
assert Rank(Matrix(K,controls)) eq 2;

// Evaluate the formal logarithm modulo 23 after translating to a split chart.
function Small(poly)
    return IsZero(poly) or Min([Valuation(z): z in Coefficients(poly)]) ge 15;
end function;

function ReduceDivisor(u,v,ff)
    while Degree(u) gt 2 do
        q,rem := Quotrem(ff-v^2,u);
        assert Small(rem);
        u := q/LeadingCoefficient(q);
        v := (-v) mod u;
    end while;
    assert Small((ff-v^2) mod u);
    return <u,v>;
end function;

function Plus(A,B,ff)
    u,v := Explode(A);
    w,z := Explode(B);
    gcd,inv,_ := XGCD(u mod w,w);
    assert Degree(gcd) eq 0;
    inv /:= Coefficient(gcd,0);
    return ReduceDivisor(u*w,v+u*(((z-v)*inv) mod w),ff);
end function;

function Twice(A,ff)
    u,v := Explode(A);
    q,rem := Quotrem(ff-v^2,u);
    assert Small(rem);
    gcd,inv,_ := XGCD(2*v,u);
    assert Degree(gcd) eq 0;
    inv /:= Coefficient(gcd,0);
    return ReduceDivisor(u^2,v+u*((q*inv) mod u),ff);
end function;

// Each entry lists two abscissas and their two ordinate residues at 23.
charts := [
    [[1,6,11,6],[1,6,11,6],[1,6,11,6],[1,8,11,4]],
    [[1,2,1,8],[1,2,1,8],[1,2,1,8],[4,6,5,3]],
    [[3,4,2,8],[2,3,3,2],[2,3,3,2],[3,4,2,8]]
];
expected := [
    [4,11,0,15,7,4], [22,17,2,8,3,11],
    [20,4,21,20,12,22], [8,13,0,10,15,17]
];
for pair in Factorization(ideal<O|23>) do
    L,em := Completion(K,pair[1]:Precision:=200);
    j := Index([11,13,21],(Integers()!em(th)) mod 23);
    assert j ne 0;
    PL<z> := PolynomialRing(L);
    ff := Polynomial(L,[em(cf):cf in Eltseq(f)]);
    order := [640,720,408][j];
    fb := Polynomial(GF(23),[(Integers()!cf) mod 23:cf in Eltseq(ff)]);
    assert #Jacobian(HyperellipticCurve(fb)) eq order;
    for k in [1..4] do
        chart := charts[j][k];
        xs := [L|chart[1],chart[2]];
        ys := [];
        for i in [1,2] do
            ok,y := IsSquare(Evaluate(ff,xs[i]));
            assert ok;
            if Valuation(y-chart[i+2]) eq 0 then y := -y; end if;
            Append(~ys,y);
        end for;
        B := <(z-xs[1])*(z-xs[2]),ys[1]+(z-xs[1])*(ys[2]-ys[1])/(xs[2]-xs[1])>;
        sign := k eq 4 select 1 else -1;
        D := <Polynomial(L,[em(cf):cf in Eltseq(us[k])]),
              Polynomial(L,[sign*em(cf):cf in Eltseq(vs[k])])>;
        A := Plus(B,D,ff);
        bits := Reverse(Intseq(order,2));
        for bit in bits[2..#bits] do
            A := Plus(Twice(A,ff),<B[1],-B[2]>,ff);
            if bit eq 1 then A := Plus(A,D,ff); end if;
        end for;
        assert Degree(A[1]) eq 2;
        assert Min([Valuation(cf):cf in Coefficients(A[1]-B[1])]) ge 1;
        assert Min([Valuation(cf):cf in Coefficients(A[2]-B[2])]) ge 1;
        vec := [L!0,L!0];
        for i in [1,2] do
            dx := -Evaluate(A[1],xs[i])/Evaluate(Derivative(A[1]),xs[i]);
            assert Valuation(dx) ge 1;
            for t in [0,1] do vec[t+1] +:= dx*xs[i]^t/(2*ys[i]*23); end for;
        end for;
        assert Min([AbsolutePrecision(v):v in vec]) gt 0;
        assert [(Integers()!v) mod 23:v in vec] eq expected[k][2*j-1..2*j];
    end for;
end for;

F := GF(23);
D := -Matrix(F,expected[1..3]);
assert Rank(D) eq 3;
SF<h> := PowerSeriesRing(F,8);
jet := [ ZeroMatrix(F,3,6) : n in [1..3] ];
for j -> theta in [F|11,13,21] do
    c := SF!(F!1/4)+h;
    for i -> k in roots[1..3] do
        t,r := Branch(c,theta,F!k);
        s := Quad(c,t)/r;
        forms := [-2*s/(1/4-theta),2*r];
        for l in [1,2], n in [1..3] do
            jet[n][i,2*j-2+l] := [640,720,408][j]*Coefficient(forms[l]/Ht(c,t),n-1)/n;
        end for;
    end for;
end for;
U<u> := PolynomialRing(F);
StackRows := func<A,B|VerticalJoin(ChangeRing(A,U),B)>;
jetU := [ChangeRing(M,U): M in jet];
poly := Determinant(StackRows(D,jetU[1]+u*jetU[2]+u^2*jetU[3]));
assert poly mod u^3 eq 7*u^2;

relations := [];
for j in [1..3] do
    sign := DiagonalMatrix(F,[ ((l+1) div 2) eq j select -1 else 1 : l in [1..6] ]);
    mixed := ChangeRing(D*sign,U)+u*ChangeRing(jet[1]*sign,U);
    det := Determinant(StackRows(D,mixed));
    assert det eq [11,21,10][j]*u+[20,18,3][j]*u^2;
    root := [19,18,12][j];
    M := VerticalJoin(D,(D+root*jet[1])*sign);
    assert Rank(M) eq 5;
    rel := Eltseq(Basis(Nullspace(M))[1]);
    Append(~relations,[Integers()!(a/rel[6]):a in rel]);
end for;
assert relations eq [[8,7,8,15,4,1],[12,5,22,10,3,1],[19,1,21,13,8,1]];

L := GF(1181);
PL<W> := PolynomialRing(L);
assert #Roots(144*W^3-240*W^2-80*W+240) eq 1;
assert IsSquarefree(144*W^3-240*W^2-80*W+240);
thetaRoots := [L|99,198];
keySpace := CartesianPower(L,4);
emptyHits := [ {@ keySpace | @} : j in [1..3] ];
survivors := emptyHits;
sizes := [];
for e -> theta in thetaRoots do
    assert 48*theta^3-320*theta^2+52*theta-95 eq 0;
    ff := Model(W,theta);
    assert Degree(ff) eq 5 and IsSquarefree(ff);
    JJ := Jacobian(HyperellipticCurve(ff));
    ord := #JJ;
    assert ord eq [1383036,1338324][e] and Valuation(ord,23) eq 1;
    proj := func<P|(ord div 23)*P>;
    fixed := [];
    for k in roots[1..3] do
        U0,V0 := Chart(L!1/4,L!k^2/4,L!k/2,theta,W);
        Append(~fixed,proj(JJ![U0,V0]));
    end for;
    assert Order(fixed[1]) eq 23;
    fixedLogs := [Index([i*fixed[1]:i in [0..22]],a)-1:a in fixed];
    assert fixedLogs eq [[1,11,21],[1,16,17]][e];
    localHits := emptyHits;
    centers := 0;
    for c in L do
        if not IsSquare(c) or not IsSquare((25-4*c)/6) then continue; end if;
        ok,z := IsSquare((c-theta)/(1/4-theta));
        if not ok then continue; end if;
        G := W^3-(5*c+25)*W^2/3+(-80*c^2+2080*c+1816)*W/144-N(c)/144;
        entries := Roots(G);
        if &+[Integers()|a[2]:a in entries] ne 3 then continue; end if;
        if not forall{a:a in entries|IsSquare(a[1])} then continue; end if;
        centers +:= 1;
        vals := &cat[[Sqrt(a[1]):i in [1..a[2]]]:a in entries];
        logs := AssociativeArray();
        for eps in CartesianPower({-1,1},3) do
            v := &+[eps[i]*vals[i]:i in [1..3]];
            U0,V0 := Chart(c,v^2,z*v,theta,W);
            value := proj(JJ![U0,V0]);
            logs[v] := Index([i*fixed[1]:i in [0..22]],value)-1;
            assert logs[v] ge 0;
        end for;
        for g in Sym(3), eps in CartesianPower({-1,1},3) do
            a,b,d := Explode([eps[i]*vals[i^g]:i in [1..3]]);
            varying := [logs[-a-b+d],logs[a-b+d],logs[-a+b+d]];
            key := <c,a,b,d>;
            for j -> rel in relations do
                left := &+[rel[i]*fixedLogs[i]:i in [1..3]];
                right := &+[rel[i+3]*varying[i]:i in [1..3]];
                if (left+right) mod 23 eq 0 or (left-right) mod 23 eq 0 then
                    Include(~localHits[j],key);
                end if;
            end for;
        end for;
    end for;
    assert centers eq [7,4][e];
    Append(~sizes,[#v:v in localHits]);
    survivors := e eq 1 select localHits else
        [survivors[j] meet localHits[j]: j in [1..3]];
end for;
assert sizes eq [[8,20,16],[18,16,22]];
assert forall{s:s in survivors | IsEmpty(s)};
print "AUXILIARY SIGNED TRIPLES", sizes;
print "PASS: every valid beta lift forces six independent Prym divisor classes";
