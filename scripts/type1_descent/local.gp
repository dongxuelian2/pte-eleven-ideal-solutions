localclass(x,p,pi,bid)={
  my(v=idealval(nf,x,p));
  concat(ideallog(nf,x/pi^v,bid),[v]~)*Mod(1,2);
};
appendindependent(M,v)={
  my(N=matconcat([M,v]));
  if(matrank(N)>matsize(M)[2],N,M);
};
localimage(es,p,pi,bid)={
  my(scale=min(idealval(nf,es[2],p),idealval(nf,es[3],p))\2);
  my(e=es/pi^(2*scale),M=matrix(0,0),x,v,c);
  my(dim=2+(p[1]==2)*p[3]*p[4]);
  for(i=1,3,
    c=vector(2,j,if(i!=j,e[i]-e[j],prod(k=1,3,if(k==i,1,e[i]-e[k]))));
    M=appendindependent(M,concat(localclass(c[1],p,pi,bid),localclass(c[2],p,pi,bid))));
  my(maxv=vecmax([idealval(nf,e[2],p),idealval(nf,e[3],p),idealval(nf,e[2]-e[3],p)])+2*p[3]+3);
  for(attempt=1,10000,
    if(matsize(M)[2]==dim,return(M));
    v=if(attempt<60,(attempt-1)\8,random(maxv+2)-1);
    x=e[1+random(3)]+pi^v*nfbasistoalg(nf,vector(#nf.zk,k,random(60))~);
    if(prod(j=1,3,x-e[j]) && nfislocalpower(nf,p,prod(j=1,3,x-e[j]),2),
      M=appendindependent(M,concat(localclass(x,p,pi,bid),localclass(x-e[2],p,pi,bid)))));
  error("incomplete local basis",p);
};
