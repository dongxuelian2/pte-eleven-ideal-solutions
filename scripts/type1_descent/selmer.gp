\\ Run from the repository root with PARI/GP.
default(realprecision,240);
setrand(1);
bn=bnfinit(y^3-y^2-330*y-149578,1);
if(bnfcertify(bn)!=1,error("cubic class group and units"));
nf=bn.nf;
w=Mod(y,nf.pol);
r=(-w^2+20*w-5426)/7168;
if(16384*r^3+38672*r^2+33144*r+11025,error("wrong root"));
S=concat(apply(l->idealprimedec(nf,l),[2,3,5,7,11,12101,1122137]));
if(bnfsunit(bn,S)[5][1]!=1,error("S-class group"));
units=apply(z->nfbasistoalg(nf,nffactorback(nf,z)),bnfunits(bn,S)[1]);
ng=#units;
if(ng!=#S+2,error("S-unit rank"));

primesof(z)=if(z,Vec(idealfactor(nf,z)[,1]),[]);
sg(z)=if(z==0,0,nfeltsign(nf,z)[1]);
signbit(z)=(1-sg(z))/2;
signrow=Mat(apply(signbit,units))*Mod(1,2);
rootsE(b)={my(d=5*(b-r));[0,d*r,-4*d/5];};
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
    x=e[1+random(3)]+pi^v*(1+random(60)+random(20)*w+random(20)*w^2);
    if(prod(j=1,3,x-e[j]) && nfislocalpower(nf,p,prod(j=1,3,x-e[j]),2),
      M=appendindependent(M,concat(localclass(x,p,pi,bid),localclass(x-e[2],p,pi,bid)))));
  error("incomplete local basis",p);
};
selmerdim(es)={
  my(H=matid(2*ng),p,pi,bid,globalmap,image,middle,ordered);
  if(setminus(Set(concat([primesof(2),primesof(es[2]),primesof(es[3]),primesof(es[2]-es[3])])),Set(S))!=[],error("bad reduction outside S"));
  for(i=1,#S,
    p=S[i]; pi=nfbasistoalg(nf,p[2]);
    bid=idealstar(nf,idealpow(nf,p,if(p[1]==2,2*p[3]+1,1)),1,2);
    globalmap=Mat(apply(z->localclass(z,p,pi,bid),units));
    image=localimage(es,p,pi,bid);
    H=H*matker(matker(image~)~*matconcat([globalmap,0;0,globalmap])*H));
  ordered=vecsort(es,(a,b)->sg(a-b));
  middle=(ordered[1]+ordered[2])/2;
  image=Mat(vector(2,j,signbit(middle-es[j]))~)*Mod(1,2);
  H=H*matker(matker(image~)~*matconcat([signrow,0;0,signrow])*H);
  matsize(H)[2];
};
