// Run after class_bounds.gp. The place order used for the logarithms is r = 146, 51.
SetColumns(0);
SetSeed(2);
Q := Rationals();
A<t> := PolynomialRing(Q);
K<r> := NumberField(t^2+551);
e := (155+r)/64;
B<v,w,z> := PolynomialRing(K,3);
R<x> := PolynomialRing(K);

function Cubic(a,c,z)
    return a^3-5*a^2*c-9*a^2*z+10*a*c^2+15*a*c*z+24*a*z^2
           -10*c^3+10*c^2*z-20*c*z^2-16*z^3;
end function;

function SymmetricMatrix(q)
    vv := [v,w,z];
    return Matrix(K,3,3,[MonomialCoefficient(q,vv[i]*vv[j])/
               (i eq j select 1 else 2): i,j in [1..3]]);
end function;

aa := e^3-5*e^2+10*e-10;
h := (6*e^2-5*e-20)/(2*aa);
q1 := 3*z^2+v*w-e*v^2;
q2 := v^2-h*z^2;
quart := ExactQuotient(Cubic(4*z^2+v*w,v^2,z^2),v^2);
q3 := ExactQuotient(q2^2-quart/aa,q1);
matrices := [SymmetricMatrix(q): q in [q1,q2,q3]];
pencil := ChangeRing(matrices[1],R)+2*x*ChangeRing(matrices[2],R)
          +x^2*ChangeRing(matrices[3],R);
g := -Determinant(pencil);
assert Degree(g) eq 6 and IsSquarefree(g);
assert 32*e^2-155*e+192 eq 0;
assert (-6*t^2+5*t+20)^2-4*(t^3-5*t^2+10*t-10)*(9*t-10)
    eq 5*t*(32*t^2-155*t+192);

K4 := NumberField(t^4-t^3+63*t^2+313*t+1544);
K8 := NumberField(t^8-3*t^7+4*t^6-328*t^5+3937*t^4-26339*t^3
                  +91290*t^2-140562*t+86400);
factors := Factorization(g);
assert [Degree(u[1]):u in factors] eq [2,4];
assert IsIsomorphic(AbsoluteField(ext<K|factors[1][1]>),K4);
assert IsIsomorphic(AbsoluteField(ext<K|factors[2][1]>),K8);
function GeneratorBound(O,K4,K8)
    L := NumberField(O);
    if AbsoluteDegree(L) eq 4 and Discriminant(O) eq 7590025 then
        assert IsIsomorphic(L,K4);
        return 5;
    elif AbsoluteDegree(L) eq 8 and Discriminant(O) eq 943857328138240000 then
        assert IsIsomorphic(L,K8);
        return 83;
    end if;
    return MinkowskiBound(O);
end function;
search := map<PowerStructure(RngOrd)->Integers()|O:->Min(83,MinkowskiBound(O))>;
proof := map<PowerStructure(RngOrd)->Integers()|O:->GeneratorBound(O,K4,K8)>;
SetClassGroupBoundMaps(search,proof);
for delta in [-e,8-e] do
    J := Jacobian(HyperellipticCurve(delta*g));
    S := TwoSelmerGroup(J:Al:="TwoSelmerGroupOld",Fields:={K4,K8});
    assert Ngens(S)-Ngens(TwoTorsionSubgroup(J)) eq 2;
end for;

function Slope(v,w)
    return -(w^3-10*v*w^2+30*v^2*w-40*v^3-25*w+100*v)
            /(3*v*w^2-10*v^2*w+10*v^3+6*w-25*v);
end function;

// The tangent meets the quadric pencil at the roots of U.
function Chart(v,w,root,delta)
    dw := Slope(v,w);
    s := delta*(v^2-h)/root;
    dr := delta*(w+v*dw-2*e*v)/(2*root);
    ds := (2*delta*v-dr*s)/root;
    position := Vector(R,[v,w,1]);
    tangent := Vector(R,[1,dw,0]);
    value := root+x*s;
    slope := dr+x*ds;
    support := InnerProduct(tangent*(delta*pencil),tangent)-slope^2;
    U := support/LeadingCoefficient(support);
    V := ((slope*position-value*tangent)*pencil)[2] mod U;
    return U,V;
end function;

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

function LogRow(D,ff,order,p)
    L := BaseRing(Parent(ff));
    X := Parent(ff).1;
    F := GF(p);
    fb := Polynomial(F,[F!(Integers()!cf):cf in Eltseq(ff)]);
    // Translate to a split support disjoint from intermediate divisors.
    xs0 := [u:u in [0..p-1]|Evaluate(fb,F!u) ne 0 and IsSquare(Evaluate(fb,F!u))];
    for k in [1..#xs0-1] do
        try
            xs := [L|xs0[k],xs0[k+1]];
            ys := [Sqrt(Evaluate(ff,u)):u in xs];
            B := <(X-xs[1])*(X-xs[2]),ys[1]+(X-xs[1])*(ys[2]-ys[1])/(xs[2]-xs[1])>;
            Dn := Plus(B,D,ff);
            bits := Reverse(Intseq(order,2));
            for bit in bits[2..#bits] do
                Dn := Plus(Twice(Dn,ff),<B[1],-B[2]>,ff);
                if bit eq 1 then Dn := Plus(Dn,D,ff); end if;
            end for;
            assert Min([Valuation(cf):cf in Coefficients(Dn[1]-B[1])
                      cat Coefficients(Dn[2]-B[2])]) ge 1;
            dx := [-Evaluate(Dn[1],u)/Evaluate(Derivative(Dn[1]),u):u in xs];
            ret := [&+[dx[i]*xs[i]^j/(2*ys[i]*p):i in [1,2]]: j in [0,1]];
            assert Min([AbsolutePrecision(l):l in ret]) gt 0;
            return [F!(Integers()!l):l in ret];
        catch err
            if k eq #xs0-1 then error err; end if;
        end try;
    end for;
    error "no split chart";
end function;

p := 197;
F := GF(p);
E := EllipticCurve([F|0,-775,0,153600,0]);
assert AbelianInvariants(AbelianGroup(E)) eq [2,108];
expectedLogs := [[123,77,4,92],[184,142,27,148],[75,54,2,128]];
expectedTans := [[66,116,10,167],[69,56,174,18],[148,119,66,137]];
for pair in Factorization(ideal<Integers(K)|p>) do
    L,em := Completion(K,pair[1]:Precision:=200);
    j := Index([146,51],(Integers()!em(r)) mod p);
    assert j ne 0;
    red := func<poly|Polynomial(F,[F!(Integers()!em(cf)):cf in Eltseq(poly)])>;
    for i in [1..3] do
        d := i eq 2 select 8-e else -e;
        v0 := [K|1,1,2][i];
        w0 := [K|-3,5,6][i];
        root := [e,8-e,(7+5*r)/32][i];
        assert root^2 eq d*(3+v0*w0-e*v0^2);
        U,V := Chart(v0,w0,root,d);
        assert (d*g-V^2) mod U eq 0;
        ff := Polynomial(L,[em(cf):cf in Eltseq(d*g)]);
        assert Degree(red(d*g)) eq 6 and IsSquarefree(red(d*g));
        order := #Jacobian(HyperellipticCurve(red(d*g)));
        D := <Polynomial(L,[em(cf):cf in Eltseq(U)]),Polynomial(L,[em(cf):cf in Eltseq(V)])>;
        block := [2*j-1..2*j];
        assert LogRow(D,ff,order,p) eq [F!expectedLogs[i][k]:k in block];
        dw := Slope(v0,w0);
        dr := d*(w0+v0*dw-2*e*v0)/(2*root);
        sheet := d*(v0^2-h)/root;
        ds := (2*d*v0-dr*sheet)/root;
        direction := Vector(K,[1,dw,0]);
        assert Valuation(em(d*InnerProduct(direction*matrices[3],direction)-ds^2)) eq 0;
        assert Valuation(em(2*root*(3*v0*w0^2-10*v0^2*w0+10*v0^3+6*w0-25*v0))) eq 0;
        assert &and[Valuation(em(c)) ge 0:c in &cat[Eltseq(M):M in matrices]];
        U1,V1 := Chart(v0+p,w0+p*dw,root+p*dr,d);
        u := red(U);
        assert IsSquarefree(u) and GCD(red(V),u) eq 1;
        _,inv,_ := XGCD(2*red(V),u);
        du := red((U1-U)/p);
        tangent := [-order*Coefficient(((Parent(u).1)^k*du*inv) mod u,1):k in [0,1]];
        assert tangent eq [F!expectedTans[i][k]:k in block];
    end for;
end for;
for i in [1,3] do
    rows := [expectedLogs[1],expectedLogs[3],expectedTans[i]];
    assert Determinant(Matrix(F,[a[1..3]:a in rows])) eq (i eq 1 select 7 else 100);
end for;
print "Both rank bounds and both transverse collision disks verified.";
