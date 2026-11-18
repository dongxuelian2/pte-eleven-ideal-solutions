\\ A certified factor-base generating bound for the sextic descent field.
default(parisizemax,1000000000);
p = x^6+228*x^4+12996*x^2+7077920;
B = bnfinit(p,1);
if(B.disc != -5540347111683392000, error("wrong absolute discriminant"));
if(vecmax(apply(a->idealnorm(B,a),B[5])) > 100, error("oversized factor-base ideal"));
if(bnfcertify(B,1) != 1, error("factor-base coverage failed"));
print("Sextic field: the primes of norm at most 100 generate the class group.");
