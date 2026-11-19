\\ After the explicit-formula check for Type 18.
default(parisizemax,1000000000);
pol=x^12-3*x^11-5625*x^10+44481*x^9+9843099*x^8-99843252*x^7-5998855201*x^6+67870752036*x^5+854290494639*x^4-9750057898095*x^3-2947764849192*x^2+213431600705238*x-420745260940011;
b=bnfinit(pol,1);
if(b.disc!=5^6*11^9*1231^5*2417^5 || b.sign!=[12,0],error("wrong field"));
if(vecmax(apply(P->idealnorm(b.nf,P),b[5]))>1000,error("factor base exceeds bound"));
install("bnftestprimes","vGG");
bnftestprimes(b,85000000);
print("All ideals of prime norm power up to 85000000 reduce to the factor base.");
