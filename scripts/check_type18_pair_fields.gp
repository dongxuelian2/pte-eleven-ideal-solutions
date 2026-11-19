\\ The two quadratic steps from the cubic partition field to the ordered-pair field.
default(parisizemax,1200000000);

unitbasis(bn)={
  my(nf=bn.nf,U=concat([-1],bn.fu),A=matrix(#U,0)*Mod(1,2));
  if(nf.r1!=poldegree(nf.pol) || #U!=nf.r1,error("unit count"));
  for(j=1,#U,if(denominator(nfalgtobasis(nf,U[j]))!=1 || abs(nfeltnorm(nf,U[j]))!=1,error("unit")));
  forprime(p=3,997,
    my(pp=idealprimedec(nf,p));
    for(j=1,#pp,
      my(md=nfmodprinit(nf,pp[j]));
      A=matconcat([A,vector(#U,i,1-issquare(nfmodpr(nf,U[i],md)))~*Mod(1,2)]);
      if(matrank(A)==#U,return(U))
    )
  );
  error("unit independence");
};

normdata(bn,D,U)={
  my(nf=bn.nf);
  if(vecmin(nfeltsign(nf,D))!=1 || nfeltispower(nf,D,2),error("quadratic extension"));
  my(relative=rnfdisc(nf,x^2-D)[1],places=idealfactor(nf,relative)[,1]);
  my(characters=matrix(#U,#places,i,j,(1-nfhilbert(nf,U[i],D,places[j]))/2)*Mod(1,2));
  [#places,matrank(characters)];
};

pairfield(f,original,root,expected)={
  my(a=-polcoef(f,3),b=polcoef(f,2),c=-polcoef(f,1));
  if(!polisirreducible(f) || !polisirreducible(original) || subst(original,y,Mod(root,f))!=0,error("quartic model"));
  my(sumf=subst(factor(polresultant(f,subst(f,y,x-y),y))[2,1],x,y));
  if(poldegree(sumf)!=6 || !polisirreducible(sumf),error("pair sums"));
  my(partf=subst(factor(polresultant(sumf,x-y*(a-y),y))[1,1],x,y));
  if(poldegree(partf)!=3 || !polisirreducible(partf),error("partitions"));
  my(cr=polredbest(partf,1),base=bnfinit(cr[1],1),D0=a^2-4*cr[2]);
  if(!bnfcertify(base) || base.no!=1,error("cubic class group"));
  my(lr=polredbest(sumf,1),L=bnfinit(lr[1],1),s=lr[2],v=s*(a-s));
  if(subst(partf,y,v)!=0,error("pair-complement field"));
  my(prod=(c-s*(b-a*s+s^2))/(a-2*s),D1=s^2-4*prod);
  if(subst(f,y,x)%(x^2-s*x+prod)!=0,error("root pair"));
  my(first=normdata(base,D0,unitbasis(base)),second=normdata(L,D1,unitbasis(L)));
  if([first,second]!=expected,error("ramification and norm indices"));
  print("Quadratic steps ",first," and ",second,"; ordered-pair class group 2-rank ",second[1]-1-second[2]);
};

P=y^5-165*y^4+8778*y^3-172810*y^2+1057221*y;
disc=polresultant(P-x,deriv(P),y);
pairfield(y^4-2*y^3-1137*y^2+19706*y-89410,subst(disc,x,y),(447676416*y^3+4775215104*y^2-450777710592*y+3196597907745)/18625,[[4,3],[7,3]]);
pairfield(y^4-y^3-991*y^2+13345*y-19455,P/y,(-y^3-23*y^2+585*y+5002)/73,[[3,2],[5,3]]);
