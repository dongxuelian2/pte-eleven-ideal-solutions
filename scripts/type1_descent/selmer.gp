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
read("scripts/type1_descent/local.gp");
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
