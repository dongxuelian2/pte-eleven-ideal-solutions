// Run after beta_selmer.m. Confirm the canonical forms by the tangent chart.
S<s> := PowerSeriesRing(K,10);
v := 1+s; w := S!5; delta := 8-e;
F := func<v,w|v*w^3-5*v^2*w^2+10*v^3*w-10*v^4+3*w^2-25*v*w+50*v^2-40>;
Fw := func<v,w|3*v*w^2-10*v^2*w+10*v^3+6*w-25*v>;
for i in [1..4] do w -:= F(v,w)/Fw(v,w); end for;
rw := S!delta;
for i in [1..4] do rw := (rw+delta*(3+v*w-e*v^2)/rw)/2; end for;
sw := delta*(v^2-h)/rw;
RX<X> := PolynomialRing(S);
pos := Vector(RX,[v,w,1]); dir := Vector(RX,[1,Derivative(w),0]);
M := ChangeRing(matrices[1],RX)+2*X*ChangeRing(matrices[2],RX)+X^2*ChangeRing(matrices[3],RX);
value := rw+X*sw; slope := Derivative(rw)+X*Derivative(sw);
u := InnerProduct(dir*(delta*M),dir)-slope^2; u /:= LeadingCoefficient(u);
y := ((slope*pos-value*dir)*M)[2] mod u;
u0 := Coefficient(u,0); u1 := Coefficient(u,1);
a := Coefficient(y,0); b := Coefficient(y,1); c := a-b*u1;
den := 2*(a*c+b^2*u0);
da := Derivative(u0); db := Derivative(u1);
constant := (da*c+db*b*u0)/den; linear := (-da*b+db*c+db*b*u1)/den;
forms := [constant-u1*linear,linear];  // Reciprocal pencil coordinates.
gamma := (58773-2161*r)/128000;
assert Determinant(Matrix(K,[[Coefficient(u,n):u in [rw,sw]]:n in [0,1]])) ne 0;
for j in [1,2], n in [0..3] do
    assert Coefficient(forms[j]*Fw(v,w),n) eq Coefficient([gamma*rw,-gamma*sw][j],n);
end for;

function SmallLog(D,B,xs,ys,ff,n)
    A := Plus(B,D,ff);
    bits := Reverse(Intseq(n,2));
    for bit in bits[2..#bits] do
        A := Plus(Twice(A,ff),<B[1],-B[2]>,ff);
        if bit eq 1 then A := Plus(A,D,ff); end if;
    end for;
    L := BaseRing(Parent(ff));
    ans := [L!0,L!0];
    disc := Discriminant(A[1]);
    root := Sqrt(disc);
    roots := [(-Coefficient(A[1],1)+sign*root)/2:sign in [-1,1]];
    for j in [1,2] do
        close := [z:z in roots|Valuation(z-xs[j]) gt Valuation(xs[1]-xs[2])+8];
        assert #close eq 1;
        xx := close[1]; dx := xx-xs[j];
        assert Valuation(Evaluate(A[2],xx)/ys[j]-1) gt 2;
        S<z> := PowerSeriesRing(L,9);
        w := Evaluate(ff,xs[j]+dx*z)/ys[j]^2-1;
        assert &and[Valuation(Coefficient(w,k)) ge 6*k:k in [1..6]];
        rec := &+[(-1)^k*Binomial(2*k,k)/4^k*w^k:k in [0..7]];
        for k in [0,1] do
            form := rec*(xs[j]+dx*z)^k;
            ans[k+1] +:= dx/(2*n*ys[j])*&+[Coefficient(form,l)/(l+1):l in [0..7]];
        end for;
        assert Valuation(dx/(2*n*ys[j]))+Min([0,Valuation(xs[j]),Valuation(dx)])+24 ge 8;
    end for;
    assert &and[AbsolutePrecision(z) ge 8:z in ans];
    return ans;
end function;

revU := 1+Coefficient(U0,1)*x+Coefficient(U0,0)*x^2; revU /:= LeadingCoefficient(revU);
revV := (Coefficient(V0,0)*x^3+Coefficient(V0,1)*x^2) mod revU;
assert (f-revV^2) mod revU eq 0;
for pair in Factorization(ideal<O|2>) do
    L,em := Completion(K,pair[1]:Precision:=8500); P<X> := PolynomialRing(L);
    j := (Integers()!em(r) mod 256) eq 155 select 1 else 2;
    conv := func<a|Polynomial(L,[em(cf):cf in Eltseq(a)])>;
    ff := conv(f); rt := OrderedRoots(ff)[1]; pts := [* *];
    for i in (j eq 1 select [352,4704] else [9]) do
        y := Sqrt(Evaluate(ff,L!i));
        Append(~pts,<(X-i)*(X-rt),y*(X-rt)/(i-rt)>);
    end for;
    if j eq 2 then
        uq := X^2+5/2; rem := ff mod uq;
        a := Coefficient(rem,0); b := Coefficient(rem,1); z := Sqrt(a^2+5*b^2/2);
        ok,v0 := IsSquare((a+z)/2); if not ok then v0 := Sqrt((a-z)/2); end if;
        Append(~pts,<uq,v0+X*b/(2*v0)>);
    end if;
    Append(~pts,<conv(revU),conv(revV)>);
    loc := [L|]; scale := j eq 1 select 32 else 1;
    for i in [1..3000 by 2] do
        a := L!(i*scale);
        if IsSquare(Evaluate(ff,a)) then Append(~loc,a); end if;
        if #loc gt 17 then break; end if;
    end for;
    logs := [];
    for D in pts do
        count := #logs;
        for b in [4..#loc-1] do
            try
                xs := loc[b..b+1]; ys := [Sqrt(Evaluate(ff,a)):a in xs];
                B := <(X-xs[1])*(X-xs[2]),ys[1]+(X-xs[1])*(ys[2]-ys[1])/(xs[2]-xs[1])>;
                Append(~logs,SmallLog(D,B,xs,ys,ff,2^18*9*5*7));
                break;
            catch err
                if b eq #loc-1 then error err; end if;
            end try;
        end for;
        assert #logs eq count+1;
    end for;
    N := Matrix(L,logs[1..2])^-1;
    assert Min([Valuation(a):a in Eltseq(N)]) ge -1;
    dc := Vector(L,logs[3])*N;
    tc := em(gamma)/40*Vector(L,[em(delta),-em(1-h)])*N;
    assert [(Integers()!(z/2)) mod 2:z in Eltseq(dc)] eq (j eq 1 select [1,1] else [1,0]);
    assert [(Integers()!(8*z)) mod 2:z in Eltseq(tc)] eq (j eq 1 select [1,0] else [0,0]);
    terms := [[Valuation(em(gamma*delta)*N[1,i]):i in [1,2]],
              [Valuation(em(gamma*(1-h))*N[2,i]):i in [1,2]],
              [Valuation(em(gamma)*N[2,i]):i in [1,2]]];
    assert terms eq (j eq 1 select [[-5,-2],[-5,-2],[-5,-2]] else [[1,2],[4,3],[0,-1]]);
end for;
assert Rank(Matrix(GF(2),[[1,1,1,0],[1,1,0,1],[1,0,0,0]])) eq 3;
assert Rank(Matrix(GF(2),[[1,1,1,0],[1,1,0,1],[0,0,1,0]])) eq 3;
print "The tangent and the opposite-sign displacement both miss the global logarithm space.";
