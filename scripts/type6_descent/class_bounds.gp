\\ Certify the ideal-class generation bounds used in prym.m.
default(parisizemax,1800000000);
checkfield(P, disc, bound) = {
  my(b=bnfinit(P,1));
  if(b.disc != disc || vecmax(apply(z->idealnorm(b,z),b[5])) > bound,
     error("unexpected discriminant or factor base"));
  if(bnfcertify(b,1) != 1, error("class-group generation failed"));
  print("Certified degree ",poldegree(P)," with generator norms <= ",bound);
};
checkfield(x^4-x^3+63*x^2+313*x+1544,7590025,5);
checkfield(x^8-3*x^7+4*x^6-328*x^5+3937*x^4-26339*x^3+91290*x^2-140562*x+86400,943857328138240000,83);
m=ellinit([1,1,1,-2913,80031]); red=elllocalred(m,2);
if(red[1]!=1 || red[2]!=18 || red[4]!=14,error("unexpected reduction at two"));
