SetClassGroupBounds("GRH");
R<x>:=PolynomialRing(Rationals()); K<a>:=NumberField(x^5-660*x^3+70983*x-102060); O:=MaximalOrder(K);
time Cl:=ClassGroup(O); assert #Cl eq 1;
th:=(-165*a^4+70983*a^2-255150*a-750870)/2;
S:=[];
for p in [3,5,7,773,3087467] do
 D:=Decomposition(O,p);
 for item in D do
  I:=item[1]; e:=item[2]; f:=Valuation(Norm(I),p);
  keep:=(p eq 3 and e eq 1) or (p in [5,7] and f eq 4) or (p eq 773 and f eq 1) or (p eq 3087467 and f eq 1 and Valuation(th-306176,I) gt 0);
  if keep then Append(~S,I); end if;
 end for;
end for;
print "IDEAL_COUNT",#S; assert #S eq 9;
time U,mu,B:=SUnitGroup(S:Raw:=true); B:=Eltseq(B); print "RAW_SIZE",#B,Ngens(U);
for i in [1..Ngens(U)] do
 e:=[z mod 2:z in Eltseq(mu(U.i))]; u:=&*[K!B[j]^e[j]:j in [1..#B]];
 print "GEN",i,Eltseq(u);
 print "VAL",[Valuation(u,I) mod 2:I in S];
end for;
