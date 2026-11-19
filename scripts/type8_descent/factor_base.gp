\\ Section 5.73 gives the norm bound 2^17 using the explicit formula.
f=x^12-6*x^11-66*x^10-1415*x^9-3987*x^8+128772*x^7+1543187*x^6-30854580*x^5+671675061*x^4-3249300991*x^3+46170279879*x^2-461672478855*x+1325918174175;
b=bnfinit(f,1);
if(b.disc!=5528109^5 || b.sign!=[0,6],error("wrong field"));
if(vecmax(apply(P->idealnorm(b,P),b[5]))>1000,error("factor base exceeds 1000"));
install("bnftestprimes","vGG");
bnftestprimes(b,2^17);
print("Every prime ideal of norm at most 2^17 reduces to the factor base.");
