\\ Run from the repository root with gp -s 300M -fq.
default(parisizemax,1500000000);
default(realprecision,320);
setrand(7);
read("scripts/type1_descent/local.gp");
read("scripts/type1_descent/pairing.gp");

f(B)=16384*B^3+38672*B^2+33144*B+11025;
singleton(B,Z)=64*Z^3-(1152*B+1328)*Z^2+(7232*B^2+14080*B+7564)*Z-f(B);
h=-800*y^6+2625*y^5-12050*y^4+529*y^3+54448*y^2+53744*y+14400;
ratio=poldisc(singleton(y,x))/h;
if(ratio!=2^26 || !polisirreducible(h),error("discriminant"));
small=nfinit(y^3-y^2-330*y-149578);
w=Mod(y,small.pol);rc=(-w^2+20*w-5426)/7168;
p=select(q->q[3]==2,idealprimedec(small,2))[1];
if(idealval(small,rc,p)!=-4 || nfislocalpower(small,p,1-1/rc,2),error("odd center"));
pr=polredbest(h,1);bn=bnfinit(pr[1],1);nf=bn.nf;r=pr[2];
if(subst(h,y,r) || nf.disc!=14755078257381100 || bn.no!=1 || bnfcertify(bn,1)!=1,error("sextic field certificate"));
sg(z)=if(z==0,0,nfeltsign(nf,z)[isolationIndex]);
primesof(z)=if(z,Vec(idealfactor(nf,z)[,1]),[]);
above(R)=concat(apply(l->idealprimedec(nf,l),R));
localdata(p)=idealstar(nf,idealpow(nf,p,if(p[1]==2,2*p[3]+1,1)),1,2);
squareclasses(P)={
  my(fu=bnfunits(bn,P)[1]);
  my(U=apply(v->prod(j=1,matsize(v)[1],if(v[j,2]%2,nfbasistoalg(nf,v[j,1]),1)),fu));
  my(n=#U);
  if(n!=#P+nf.sign[1]+nf.sign[2],error("S-unit count"));
  my(M=matrix(#P,n,i,j,idealval(nf,U[j],P[i])%2)*Mod(1,2));
  forprime(l=3,1200,
    if(matrank(M)==n,return(U));
    foreach(idealprimedec(nf,l),q,
      my(pi=nfbasistoalg(nf,q[2]),bid=localdata(q));
      M=matconcat([M;Mat(apply(z->localclass(z,q,pi,bid),U))])));
  error("S-unit independence");
};
centersOK(B,p)={
  my(z=Mod(B,p),fac);
  if(!issquare(z) || !issquare(5*z+4) || !issquare(f(z)),return(0));
  fac=factor(singleton(z,'x));
  prod(j=1,matsize(fac)[1],poldegree(fac[j,1])==1 && issquare(-polcoef(fac[j,1],0)));
};
primeResidue(z,p)={
  my(md=nfmodprinit(nf,p),t=nfmodpr(nf,z,md));
  if(t^p[1]!=t,return(-1));
  my(a=lift(nfbasistoalg(nf,nfmodprlift(nf,t,md))));
  if(poldegree(a)>0,error("residue representative"));
  lift(Mod(polcoef(a,0),p[1]));
};
charsupport=Set(Vec(factor(abs(poldisc(h)*pollead(h)*polcoef(h,0)))[,1]));
if(charsupport!=[2,3,5,23,41,2111,11801,144461],error("character support"));
character()={
  my(U=squareclasses(above(charsupport)),n=#U,norms=apply(z->nfeltnorm(nf,z),U));
  my(ns=concat(charsupport,[-1]),N=matrix(#ns,n,i,j,if(ns[i]==-1,norms[j]<0,valuation(norms[j],ns[i])%2)));
  my(H=matker(N*Mod(1,2)),pi,bid,V,pp,rv,res,pool,base,span);
  foreach(above([2,73]),p,
    if((p[1]==2 && 6*p[3]-idealval(nf,r,p)<=2*p[3]) || (p[1]==73 && idealval(nf,r,p)),error("anchored square"));
    pi=nfbasistoalg(nf,p[2]);bid=localdata(p);
    V=Mat(apply(z->localclass(z,p,pi,bid),U));H=H*matker(V*H));
  for(j=2,#charsupport,my(l=charsupport[j]);pp=idealprimedec(nf,l);
    V=matrix(#pp,n,i,k,idealval(nf,U[k],pp[i])%2)*Mod(1,2);
    rv=vector(#pp,i,idealval(nf,r,pp[i]));
    res=vector(#pp,i,if(rv[i]>=0,primeResidue(r,pp[i]),-1));
    base=vector(#pp,i,max(rv[i],0)%2)~*Mod(1,2);
    pool=List([base,vector(#pp,i,rv[i]%2)~*Mod(1,2)]);
    for(i=1,#pp,if(res[i]>=0 && centersOK(res[i],l),listput(pool,vector(#pp,k,k==i)~*Mod(1,2))));
    span=Mat(Vec(pool));H=H*matker(matker(span~)~*V*H));
  if(matsize(H)[2]!=0,error("nontrivial center character"));
};
character();
print("Discriminant character (b-r)/(-r) is trivial.");
read("scripts/type1_descent/sextic_tangents.gp");
read("scripts/type1_descent/sextic_four.gp");
