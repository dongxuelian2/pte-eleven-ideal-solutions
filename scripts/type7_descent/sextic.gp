\\ The two isogeny Selmer spaces and a lift of the dyadic kernel class.
default(realprecision,180);
setrand(5);
F=y^6-96*y^4+144*y^2-54;
bn=bnfinit(F,1); nf=bn.nf;
if(bnfcertify(bn)!=1 || bn.no!=1,error("sextic class group and units"));
n=Mod(4/27*y^5-127/9*y^3+32/3*y,F);
r=Mod(-5/81*y^4+53/9*y^2-41/9,F);
pminus(z)=9*z^3+5*z^2-3*z-1;
k=(n-r)*pminus(n); m=k*pminus(r); j0=-m/(2*r);
a0=k*(27*r^2+10*r-3); a=a0+3*j0; b=3*j0^2+2*j0*a0+k*m*(27*r+5);
disc2=a^2-4*b; coeff=[[a,b],[-2*a,disc2]];
read("scripts/type1_descent/local.gp");
support(z)=if(z,Vec(idealfactor(nf,z)[,1]),[]);
sg(z,i)=if(z==0,0,nfeltsign(nf,z)[i]);
positiveonly(ab,i)=sg(ab[2],i)>0 && (sg(ab[1],i)<0 || sg(ab[1]^2-4*ab[2],i)<0);
bad=Set(concat([support(2),support(b),support(disc2)]));
su=apply(t->nfbasistoalg(nf,nffactorback(nf,t)),bnfunits(bn,bad)[1]);
ss=apply(v->nfbasistoalg(nf,Col(v)),[[52,-50,-42,13,-16,13],[5,3,0,-2,-2,0],[1,4,2,1,0,0],[4,1,2,1,0,1]]);
expect=Mod(-8/9*y^5+10/9*y^4+254/3*y^3-317/3*y^2-66*y+79,F);
if(ss[2]!=expect,error("number-field integral basis"));
bidat(p)=idealstar(nf,idealpow(nf,p,if(p[1]==2,2*p[3]+1,1)),1,2);
uniform(p)=nfbasistoalg(nf,p[2]);
localimages(p,pi,bid)={
  my(ab=vector(2,j,my(h=min(idealval(nf,coeff[j][1],p)\2,idealval(nf,coeff[j][2],p)\4));[coeff[j][1]/pi^(2*h),coeff[j][2]/pi^(4*h)]));
  my(target=2+(p[1]==2)*p[3]*p[4],H=vector(2,j,Mat(localclass(ab[j][2],p,pi,bid))));
  for(attempt=1,160000,
    if(matrank(H[1])+matrank(H[2])==target,return(H));
    for(j=1,2,
      my(v=random(max(0,idealval(nf,ab[j][2],p))+20)-10);
      my(X=pi^v*nfbasistoalg(nf,vector(6,i,random(60))~),f=X*(X^2+ab[j][1]*X+ab[j][2]));
      if(f && nfislocalpower(nf,p,f,2),
        my(M=matconcat([H[j],localclass(X,p,pi,bid)]));if(matrank(M)>matrank(H[j]),H[j]=M))));
  error("incomplete pair of local images",p);
};
bounds()={
  my(H=vector(2,i,matid(#su)*Mod(1,2)));
  for(i=1,#bad,
    my(p=bad[i],pi=uniform(p),bid=bidat(p),L=localimages(p,pi,bid),mp=Mat(apply(t->localclass(t,p,pi,bid),su)));
    for(j=1,2,if(matrank(L[j])<matsize(mp)[1],H[j]=H[j]*matker(matker(L[j]~)~*mp*H[j]))));
  for(i=1,nf.r1,for(j=1,2,if(positiveonly(coeff[j],i),
    my(row=Mat(apply(t->(1-sg(t,i))/2,su))*Mod(1,2));H[j]=H[j]*matker(row*H[j]))));
  H;
};
H=bounds();
if(apply(matrank,H)!=[4,1],error("isogeny descent dimensions"));
prodmask(V,i)=prod(j=1,#V,if(bittest(i,j-1),V[j],1));
insquare(z)=#nfroots(nf,x^2-z)>0;
basis1=vector(4,i,prod(j=1,#su,su[j]^lift(H[1][j,i])));
masks=vector(4,j,select(i->insquare(ss[j]*prodmask(basis1,i)),[0..15]));
if(vecsort(apply(v->#v,masks))!=[1,1,1,1],error("class representatives"));
if(matrank(matrix(4,4,i,j,bittest(masks[j][1],i-1))*Mod(1,2))!=4,error("dependent representatives"));
if(insquare(disc2),error("dual torsion class"));

\\ Affine lift verification for s=ss[2].  W=1 at the chosen conic point.
s=ss[2];kap=Mod(314/3*y^5+70*y^4-29924/3*y^3-6677*y^2+8049*y+5583,F);
U=Mod(1222/15+4438/45*y-3581/30*y^2-13061/90*y^3+169/135*y^4+247/162*y^5,F);
V=Mod(109023/194+146688/485*y-152419/194*y^2-466127/970*y^3+16011/1940*y^4+9801/1940*y^5,F);
if(s*U^2+a*U*V+b/s*V^2!=1,error("conic"));
lin=[2*s*U+a*V,a*U+2*b/s*V,-2];
affinebad=Set(concat([bad,support(a),support(s),support(kap),support(idealadd(nf,idealadd(nf,lin[1],lin[2]),lin[3]))]));
quart(u)=s*u^4+a*u^2+b/s;
localroot(t,p,pi,md)={
  my(v=idealval(nf,t,p),un=t/pi^v,e=idealval(nf,2,p),z=0);
  if(v%2,error("odd square valuation"));
  if(p[1]!=2,z=nfbasistoalg(nf,nfmodprlift(nf,sqrt(nfmodpr(nf,un,md)),md)),
    if(p[4]!=1,error("dyadic residue degree"));
    for(mask=0,2^(2*e+1)-1,z=sum(i=0,2*e,bittest(mask,i)*pi^i);if(idealval(nf,z^2-un,p)>2*e,break)));
  if(idealval(nf,z^2-un,p)<=2*e,error("square root seed"));
  while(idealval(nf,z^2-un,p)<72,z=(z+un/z)/2;z=nfbasistoalg(nf,idealchinese(nf,Mat([p,72+2*e]),[z])));
  z*pi^(v/2);
};
tangentvalue(p,pi)={
  my(md=nfmodprinit(nf,p),h=(idealval(nf,b,p)-2*idealval(nf,s,p))\4,e=idealval(nf,2,p));
  for(i=1,60000,
    my(u=pi^(h-5+random(12+2*p[3]))*nfbasistoalg(nf,vector(6,k,random(60))~),t=quart(u));
    if(!t || !nfislocalpower(nf,p,t,2),next);
    my(z=localroot(t,p,pi,md),L=lin[1]*u^2+lin[2]+lin[3]*z);
    if(!L,next);
    my(err=idealval(nf,z^2-t,p)-idealval(nf,z,p)-e);
    if(err>idealval(nf,z,p)+e && idealval(nf,lin[3],p)+err>idealval(nf,L,p)+2*e,return(L)));
  error("local tangent value",p);
};
checkfinite()={
  for(i=1,#affinebad,
    my(p=affinebad[i],pi=uniform(p),bid=bidat(p),L=localimages(p,pi,bid)[2],t=tangentvalue(p,pi));
    if(matrank(matconcat([L,localclass(kap*t,p,pi,bid)]))!=matrank(L),error("affine condition",p)));
};
checkreal()={
  for(i=1,nf.r1,if(positiveonly(coeff[2],i),
    my(ok=0);
    for(j=0,161,
      my(u=if(j==0,0,if(j%2,2^((j+1)\2),1/2^(j\2))),v=quart(u),h=lin[1]*u^2+lin[2]);
      if(j==161 && sg(-a/(2*s),i)>0,
        u=round(2^80*sqrt(nfeltembed(nf,-a/(2*s))[i]))/2^80;v=quart(u);h=lin[1]*u^2+lin[2]);
      if(sg(v,i)<=0 || sg(h^2-lin[3]^2*v,i)==0,next);
      my(signL=if(sg(h^2-lin[3]^2*v,i)>0,sg(h,i),sg(lin[3],i)));
      if(signL!=sg(kap,i),error("real tangent class"));ok=1;break);
    if(!ok,error("real reference point"))));
};
checkfinite();checkreal();
print("Certified sextic isogeny bounds 4 and 1 and the affine lift of the local kernel.");
