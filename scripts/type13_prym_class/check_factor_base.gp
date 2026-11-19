\\ PARI/GP 2.19.1. This calls the library's unconditional relation check.
\\ Run with gp -fq scripts/type13_prym_class/check_factor_base.gp.
default(parisizemax,1000000000);
pol = x^12-4*x^11-68*x^10-18698*x^9+336998*x^8-500324*x^7+7039932*x^6-194015520*x^5+1041300636*x^4-3497537936*x^3+20028904184*x^2-109223576168*x+269574001992;
b = bnfinit(pol,1);
if(b.disc != -172707117717588564754287819295488000000, error("wrong field discriminant"));
if(vecmax(apply(a->idealnorm(b.nf,a),b[5])) > 1000, error("factor base exceeds 1000"));
install("bnftestprimes", "vGG");
bnftestprimes(b,48000001);
print("PASS: every prime ideal of norm at most 48000001 reduces to the factor base");
