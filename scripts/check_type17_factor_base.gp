\\ PARI/GP 2.19.1. The explicit formula in section 5.60 gives the bound 2^20.
default(parisizemax,1000000000);
pol = x^12-472*x^10+75282*x^8-4450612*x^6+64570537*x^4-277288596*x^2+234027000;
b = bnfinit(pol,1);
if(b.disc != 58160039859394517078535559761100800000 || b.sign != [12,0], error("unexpected field"));
if(vecmax(apply(a->idealnorm(b,a),b[5])) > 1000, error("factor base exceeds 1000"));
install("bnftestprimes", "vGG");
bnftestprimes(b,2^20);
print("All prime ideals of norm at most 2^20 reduce to the factor base.");
