\\ Cassels' tangent-plane formula for split two-torsion.
\\ All symbol arguments below are exact field elements.
isolationIndex=1;
localroot(u,p,pi,md,precision)={
  my(v=idealval(nf,u,p),unit=u/pi^v,e=p[3],z,fac);
  if(v%2,error("odd square valuation"));
  if(p[1]!=2,
    z=nfbasistoalg(nf,nfmodprlift(nf,sqrt(nfmodpr(nf,unit,md)),md)),
    my(f=p[4],g=ffgen(nfmodpr(nf,1,md)));
    my(reps=vector(2^f,k,nfbasistoalg(nf,nfmodprlift(nf,sum(i=0,f-1,bittest(k-1,i)*g^i),md))));
    for(mask=0,2^(f*(2*e+1))-1,
      z=sum(k=0,2*e,reps[1+((mask>>(f*k))%2^f)]*pi^k);
      if(idealval(nf,z^2-unit,p)>2*e,break)));
  if(idealval(nf,z^2-unit,p)<=if(p[1]==2,2*e,0),error("initial Hensel root"));
  fac=Mat([p,precision+2*e]);
  while(idealval(nf,z^2-unit,p)<precision,
    z=nfbasistoalg(nf,idealchinese(nf,fac,[(z+unit/z)/2])));
  z*pi^(v/2);
};
localx(es,ds,p,pi)={
  my(vs=apply(z->idealval(nf,z,p),[es[2],es[3],es[2]-es[3]]),v,x);
  for(attempt=1,50000,
    v=vecmin(vs)-3+random(vecmax(vs)-vecmin(vs)+12+2*p[3]);
    x=es[1+random(3)]+pi^v*nfbasistoalg(nf,vector(#nf.zk,k,random(60))~);
    if(!setsearch(Set(es),x) && prod(k=1,3,nfislocalpower(nf,p,(x-es[k])/ds[k],2)),return(x)));
  error("missing local lift",p);
};
interval(z)={
  my(t=2^256,k=floor(real(nfeltembed(nf,z)[isolationIndex])*t));
  if(sg(z-k/t)<0 || sg((k+1)/t-z)<0,error("real interval"));
  [k/t,(k+1)/t];
};
sqrtinterval(z)={
  my(I=interval(z),t=2^106);
  [sqrtint(floor(I[1]*t^2))/t,(1+sqrtint(ceil(I[2]*t^2)))/t];
};
intervalprod(I,J)={
  my(v=[I[1]*J[1],I[1]*J[2],I[2]*J[1],I[2]*J[2]]);
  [vecmin(v),vecmax(v)];
};
realrow(es,ds,co)={
  my(low=select(j->sg(ds[j])>0,[1,2,3]),high=select(j->sg(ds[j])<0,[1,2,3]));
  my(a=if(#low,vecsort(apply(j->es[j],low),(u,v)->sg(u-v))[#low],-oo));
  my(b=if(#high,vecsort(apply(j->es[j],high),(u,v)->sg(u-v))[1],oo));
  my(x=if(a==-oo,b-1,if(b==oo,a+1,(a+b)/2)));
  if(!prod(i=1,3,sg((x-es[i])/ds[i])>0),error("real lift"));
  my(z=vector(3,i,sqrtinterval((x-es[i])/ds[i])));
  my(L=vector(3,j,intervalprod(interval(co[j][1]),z[j%3+1])+intervalprod(interval(co[j][2]),z[(j+1)%3+1])+interval(co[j][3])));
  if(vecmin(apply(v->v[1]*v[2],L))<=0,error("real tangent interval"));
  vector(3,j,sign(L[j][1]));
};
pairmatrix(es,reps,points)={
  my(ds=apply(row->apply(z->Mod(z,nf.pol),row),reps),n=#ds,M=matrix(n,n));
  my(places=Set(S),fa,co,P,jj,kk,p,pi,md,x,z,L,err,prec=64,ss);
  for(i=1,n,
    if(#nfroots(nf,'x^2-prod(j=1,3,ds[i][j]))!=2,error("Kummer norm"));
    for(j=1,3,
      fa=idealfactor(nf,ds[i][j]);
      for(k=1,matsize(fa)[1],if(!setsearch(Set(S),fa[k,1]) && fa[k,2]%2,error("ramification")));
      places=setunion(places,Set(Vec(fa[,1])))));
  for(i=1,n,
    co=vector(3,j,jj=j%3+1;kk=jj%3+1;P=apply(z->Mod(z,nf.pol),points[i][j]);
      my(a=[ds[i][jj],-ds[i][kk],es[jj]-es[kk]]);
      if(sum(k=1,3,a[k]*P[k]^2),error("tangent conic"));
      vector(3,k,a[k]*P[k]));
    my(T=places);
    for(j=1,3,T=setunion(T,Set(primesof(idealadd(nf,idealadd(nf,co[j][1],co[j][2]),co[j][3])))));
    ss=realrow(es,ds[i],co);
    my(signs=vector(n,k,prod(j=1,3,if(ss[j]<0,sg(ds[k][j]),1))));
    for(h=1,#T,
      p=T[h];pi=nfbasistoalg(nf,p[2]);md=nfmodprinit(nf,p);
      x=localx(es,ds[i],p,pi);
      z=vector(3,j,localroot((x-es[j])/ds[i][j],p,pi,md,prec));
      L=vector(3,j,co[j][1]*z[j%3+1]+co[j][2]*z[(j+1)%3+1]+co[j][3]);
      for(j=1,3,
        err=min(idealval(nf,co[j][1]*z[j%3+1],p),idealval(nf,co[j][2]*z[(j+1)%3+1],p))+prec-(p[1]==2)*p[3];
        if(err<=idealval(nf,L[j],p)+(p[1]==2)*2*p[3],error("tangent precision")));
      for(k=1,n,signs[k]*=prod(j=1,3,nfhilbert(nf,L[j],ds[k][j],p))));
    for(k=1,n,M[i,k]=(1-signs[k])/2));
  if(M!=M~ || vecmax(vector(n,i,abs(M[i,i])))>0,error("alternation"));
  M*Mod(1,2);
};
