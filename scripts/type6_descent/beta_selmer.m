// Run after prym.m in the same Magma process. X here is the reciprocal
// of the pencil coordinate, so the leading coefficient is supported over 2,3,5.
f := &+[Coefficient((8-e)*g,6-i)*x^i:i in [0..6]];
monic := f/LeadingCoefficient(f);
Al<theta> := quo<R|monic>;
_,_ := AbsoluteAlgebra(Al:Fields:={K4,K8});
O := Integers(K);
bad := [2,3,5,7,31];
primes := [t[1]:t in &cat[Factorization(ideal<O|p>):p in bad]];
Sbad := Seqset(primes);
cl,clmap := ClassGroup(O);
assert #Subfields(K8,2) eq 1 and #Subfields(K8,4) eq 0;
assert #cl eq 26 and Order(primes[1] @@ clmap) eq 26;
assert Seqset(PrimeDivisors(Numerator(Abs(Norm(Discriminant(f))))) cat
              PrimeDivisors(Denominator(Norm(Discriminant(f))))) subset Seqset(bad);
for cf in Coefficients(monic) do
    assert Seqset(PrimeDivisors(Denominator(cf))) subset Seqset(bad);
end for;
G,mu := pSelmerGroup(Al,2,Sbad);
GK,nu := pSelmerGroup(2,Sbad);
assert Ngens(G) eq 40 and Ngens(GK) eq 11;
basis := [G.i @@ mu:i in [1..Ngens(G)]];
data := [Eltseq(nu(Resultant(monic,R!Eltseq(a)))):a in basis]
        cat [Eltseq(nu(LeadingCoefficient(f)))];
F2 := GF(2);
blocks := [
  <2,155,[Rationals()|224,352,1632,2656,4704,18016]>,
  <2,101,[Rationals()|9,45,101]>,
  <3,2,[Rationals()|5,7/3,8]>, <3,1,[Rationals()|2,8/3,27]>,
  <5,2,[Rationals()|9]>, <5,3,[Rationals()|4]>,
  <7,4,[Rationals()|4,3/7]>, <7,3,[Rationals()|1,2/7,1/7]>,
  <31,10,[Rationals()|1,3,7]>, <31,21,[Rationals()|7/31,29/31]>
];
fac := Factorization(monic);
u2 := fac[1][1]; q4 := fac[2][1];
torpol := CRT([(-LeadingCoefficient(f)*q4) mod u2,u2 mod q4],[u2,q4]);
U0,V0 := Chart(K!1,K!5,8-e,8-e);
collision := 1+Coefficient(U0,1)*theta+Coefficient(U0,0)*theta^2;
torclass := Evaluate(torpol,theta);
for I in primes do
    p := Minimum(I);
    L,em := Completion(K,I:Precision:=200);
    key := Integers()!em(r) mod (p eq 2 select 256 else p);
    xs := [b[3]:b in blocks|b[1] eq p and b[2] eq key][1];
    ff := Polynomial(L,[em(a):a in Eltseq(f)]);
    lm := LocalTwoSelmerMap(Al,I);
    W := VectorSpace(F2,Ngens(Codomain(lm)));
    loc := func<a|W!Eltseq(lm(Al!a))>;
    scalars := p eq 2 select [K|-1,2,5] else [K|p,Integers()!PrimitiveElement(GF(p))];
    origin := loc(xs[1]-theta);
    assert &and[IsSquare(Evaluate(ff,L!a)):a in xs];
    H := sub<W|[loc(a):a in scalars] cat [loc(a-theta)-origin:a in xs]
               cat [loc(torclass),loc(collision)]>;
    if p eq 2 and key eq 101 then
        // X = sqrt(-10)/2. Test whether f(X) is a square in this quadratic field.
        X := Parent(ff).1;
        rem := ff mod (X^2+5/2);
        a := Coefficient(rem,0); b := Coefficient(rem,1);
        z := Sqrt(a^2+5*b^2/2);
        assert b ne 0 and (IsSquare((a+z)/2) or IsSquare((a-z)/2));
        H := sub<W|H,loc(theta^2+5/2)>;
    end if;
    assert Dimension(H) eq #Factorization(ff)+(p eq 2 select 3 else 0);
    _,proj := quo<W|H>;
    block := [Eltseq(proj(loc(a))):a in basis] cat [Eltseq(proj(origin))];
    data := [data[i] cat block[i]:i in [1..#data]];
end for;
N := Nullspace(Matrix(F2,data));
assert Dimension(N) eq 14;
odd := [a:a in Basis(N)|a[#data] eq 1][1];
oddRep := &*[Al|basis[i]:i in [1..#basis]|odd[i] eq 1];

// Identify the odd coset modulo local torsion at the two dyadic places.
function OrderedRoots(ff)
    roots := [-Coefficient(q[1],0):q in Factorization(ff)|Degree(q[1]) eq 1];
    Sort(~roots,func<a,b|(Integers()!(a*32) mod 2^32)-(Integers()!(b*32) mod 2^32)>);
    return roots;
end function;
for pair in Factorization(ideal<O|2>) do
    I := pair[1]; L,em := Completion(K,I:Precision:=560);
    ff := Polynomial(L,[em(cf):cf in Eltseq(f)]);
    roots := OrderedRoots(ff); lm := LocalTwoSelmerMap(Al,I);
    W := VectorSpace(F2,Ngens(Codomain(lm)));
    loc := func<a|W!Eltseq(lm(Al!a))>;
    VV,proj := quo<W|[loc(a):a in [-1,2,5]]>;
    vclass := func<a|proj(loc(a))>;
    // At the chosen root the squareclass is that of the derivative. For the
    // other roots the shifted polynomial bounds each relative displacement.
    near := [];
    for rt in roots do
        fp := Evaluate(Derivative(ff),rt);
        a := rt+2^160/fp;
        q := 2^Max(0,-Valuation(a));
        xx := Rationals()!((Integers()!(q*a)) mod 2^300)/q;
        assert IsSquare(Evaluate(ff,L!xx)) and Valuation((xx-rt)*fp/2^160-1) ge 3;
        assert &and[Valuation((xx-rt)^j*Coefficient(Evaluate(ff,Parent(ff).1+rt),j+1)/fp) ge 3*j:
                    j in [1..5]];
        Append(~near,xx);
    end for;
    w0 := vclass(near[1]-theta);
    rows := [vclass(a-theta)-w0:a in near[2..#near]];
    j := (Integers()!em(r) mod 256) eq 155 select 1 else 2;
    if j eq 1 then
        rows cat:= [vclass(a-theta)-w0:a in [352,4704]];
    else
        rows cat:= [vclass(9-theta)-w0,vclass(theta^2+5/2)];
    end if;
    M := Matrix(rows);
    assert Rank(M) eq #Factorization(ff);
    assert Eltseq(Solution(M,vclass(oddRep)-w0))[#rows-1..#rows] eq
           (j eq 1 select [F2|1,1] else [F2|0,1]);
    assert Eltseq(Solution(M,vclass(collision))) eq
           (j eq 1 select [F2|0,1,1,0,0] else [F2|0,0,0]);
    cc := Solution(M,vclass(torclass));
    if j eq 1 then assert cc eq Solution(M,vclass(collision));
    else assert Eltseq(cc) eq [F2|1,0,0]; end if;
end for;
print "The degree-one coset gives the half-class (1,1,0,1) modulo local torsion.";
