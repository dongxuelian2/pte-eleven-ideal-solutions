// The 4b=1 two-cover class, with X equal to four times a singleton square.
Q := Rationals();
A<B,X> := PolynomialRing(Q,2);
P := (X-B)*((X-B-4)^2-16*B)*((X-B-16)^2-64*B);
F := ExactQuotient(P-Evaluate(P,[B,1]),X-1);
assert Evaluate(F,[1,X]) eq (X-25)*(X-9)^2*(X-1);
D<U,V> := PolynomialRing(Q,2);
shifted := Evaluate(F,[1+U,9+V]);
assert &+[MonomialCoefficient(shifted,m)*m : m in Monomials(shifted) | TotalDegree(m) le 2]
       eq -64*(2*V^2-15*U^2);
k := GF(23); R<t> := PolynomialRing(k);
sq := {a^2:a in k}; seeds := {};
for s in sq do
    p := (t-s)*((t-s-4)^2-16*s)*((t-s-16)^2-64*s);
    g := ExactQuotient(p-Evaluate(p,1),t-1);
    factors := Roots(g);
    if &+[Integers()|a[2]:a in factors] eq 4 and &and[a[1] in sq:a in factors] then
        Include(~seeds,s);
    end if;
end for;
assert seeds eq {k|1,2,9};
assert {s:s in seeds | (s+4)*5 in sq} eq {k!1};
assert Evaluate(t^4-18*t^3+114*t^2-97*t+405,-4) eq 0;
assert Evaluate(4*t^3-54*t^2+228*t-97,-4) ne 0;
assert k!30 notin sq;
print "The branch represented by b=1/4 contains only the collision.";
