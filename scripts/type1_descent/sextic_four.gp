es=[0,-5*r^2,4*r];
S=above([2,3,5,41,2111,11801,144461]);
if(setminus(Set(concat([primesof(2),primesof(es[2]),primesof(es[3]),primesof(es[2]-es[3])])),Set(S))!=[],error("bad reduction"));
rows=apply(row->apply(t->Mod(t,nf.pol),row),sexticReps);
ds=rows[6];n=5;
for(i=1,6,if(!#nfroots(nf,'x^2-prod(k=1,3,rows[i][k])),error("Kummer norm")));
mask(n)=vector(5,k,bittest(n,k-1))~*Mod(1,2);
rep(n)=vector(2,j,prod(k=1,5,if(bittest(n,k-1),rows[k][j],1)));
same(a,b)=prod(j=1,2,#nfroots(nf,'x^2-a[j]/b[j])>0);
for(i=1,31,if(same(rep(i),[1,1]),error("dependent basis")));
if(!same(ds,rep(13)),error("primitive class representative"));
{for(i=1,5,for(j=1,3,my(F=idealfactor(nf,rows[i][j]));
  for(k=1,matsize(F)[1],if(!setsearch(Set(S),F[k,1]) && F[k,2]%2,error("ramified basis")))));}
tors=[27,7,28];
{for(i=1,3,my(v=vector(2,j,if(i==j,prod(k=1,3,if(i==k,1,es[i]-es[k])),es[i]-es[j])));
  if(!same(rep(tors[i]),v),error("torsion masks")));}
realimage()={
  my(E=vecsort(es,(a,b)->sg(a-b)),mid=(E[1]+E[2])/2);
  Mat(vector(2,j,sg(mid-es[j])<0)~)*Mod(1,2);
};
split(z,p,pi,bid)=concat(localclass(z[1],p,pi,bid),localclass(z[2],p,pi,bid));
checkselmer()={
  my(U=squareclasses(S),H=matid(2*#U),p,pi,bid,M,N,G);
  foreach(S,p,pi=nfbasistoalg(nf,p[2]);bid=localdata(p);
    M=localimage(es,p,pi,bid);N=matker(M~)~;
    if(N*Mat(apply(z->split(z,p,pi,bid),rows))!=0,error("local Selmer representatives"));
    G=Mat(apply(z->localclass(z,p,pi,bid),U));H=H*matker(N*matconcat([G,0;0,G])*H));
  for(v=1,nf.sign[1],isolationIndex=v;N=matker(realimage()~)~;
    if(N*Mat(apply(z->vector(2,k,sg(z[k])<0)~,rows))!=0,error("real representatives"));
    G=Mat(apply(z->sg(z)<0,U))*Mod(1,2);H=H*matker(N*matconcat([G,0;0,G])*H));
  if(matsize(H)[2]!=5,error("Selmer dimension"));
};
checkselmer();
print("Sextic elliptic cover: Selmer dimension 5, rank at most 3.");

\\ The final abscissa in each row is the reference lift on the 2-cover.
localpoints=eval(concat(readstr("scripts/type1_descent/sextic_local_points.txt")));
localpoints=apply(row->apply(z->Mod(Polrev(z,y),nf.pol),row),localpoints);
places2=[select(p->p[3]==2,idealprimedec(nf,2))[1],select(p->p[4]==3,idealprimedec(nf,2))[1]];
scales=apply(p->2^(-idealval(nf,r,p)/p[3]),places2);
localbasis(i)={
  my(p=places2[i],pi=nfbasistoalg(nf,p[2]),bid=localdata(p),roots=scales[i]^2*es);
  my(xs=localpoints[i],M=matrix(0,0));
  for(k=1,#xs-1,my(x=xs[k],rhs=prod(j=1,3,x-roots[j]));
    if(rhs && !nfislocalpower(nf,p,rhs,2),error("basis abscissa"));
    my(s=vector(2,j,if(roots[j]==x,prod(h=1,3,if(j==h,1,x-roots[h])),x-roots[j])));
    M=matconcat([M,split(s,p,pi,bid)]));
  if(matrank(M)!=2+p[3]*p[4],error("local basis incomplete"));
  M;
};
LB=vector(2,i,localbasis(i));
\\ Reduced log coordinates, checked in sextic_logs.m.
logs=[[0,1,0,0;0,0,1,0],[0,0,1,0,0;0,0,0,1,0;0,0,0,0,1]];
refbits=[[1,0]~,[1,1,1]~];
logimage(i)={
  my(p=places2[i],pi=nfbasistoalg(nf,p[2]),bid=localdata(p));
  logs[i]*matinverseimage(LB[i],Mat(vector(5,k,split(rows[k],p,pi,bid))));
};
joint=matconcat([logimage(1);logimage(2)]);
ray=vector(5,k,1)~*Mod(1,2);
if(matrank(joint)!=2 || matrank(matconcat([joint,ray]))!=3,error("first logarithm test"));
kernelMasks=select(i->joint*mask(i)==0,[0..31]);
if(kernelMasks!=Set(concat(concat([0],tors),apply(t->bitxor(t,13),concat([0],tors)))),error("primitive kernel"));

tang=vector(3,j,my(a=j%3+1,b=a%3+1,co=[ds[a],-ds[b],es[a]-es[b]],P=apply(z->Mod(z,nf.pol),sexticTangents[j]));if(sum(k=1,3,co[k]*P[k]^2),error("tangent conic"));vector(3,k,co[k]*P[k]));
support=Set(S);
{for(j=1,3,support=setunion(support,Set(primesof(ds[j])));
  support=setunion(support,Set(primesof(idealadd(nf,idealadd(nf,tang[j][1],tang[j][2]),tang[j][3])))));}
affine()={
  my(U=squareclasses(support),n=#U,eq=List(),fix=vector(2),p,pi,bid,md,jj,x,z,L,M,G,F,err,prec=96);
  foreach(support,p,pi=nfbasistoalg(nf,p[2]);bid=localdata(p);md=nfmodprinit(nf,p);
    jj=select(i->places2[i]==p,[1,2]);
    x=if(#jj,localpoints[jj[1]][#localpoints[jj[1]]]/scales[jj[1]]^2,localx(es,ds,p,pi));
    z=vector(3,k,localroot((x-es[k])/ds[k],p,pi,md,prec));
    L=vector(3,k,tang[k][1]*z[k%3+1]+tang[k][2]*z[(k+1)%3+1]+tang[k][3]);
    for(k=1,3,err=min(idealval(nf,tang[k][1]*z[k%3+1],p),idealval(nf,tang[k][2]*z[(k+1)%3+1],p))+prec-(p[1]==2)*p[3];
      if(err<=idealval(nf,L[k],p)+(p[1]==2)*2*p[3],error("tangent precision")));
    F=split([L[2]*L[3],L[1]*L[3]],p,pi,bid);
    M=if(#jj,LB[jj[1]],localimage(es,p,pi,bid));
    G=Mat(apply(u->localclass(u,p,pi,bid),U));G=matconcat([matconcat([G,0;0,G]),F]);
    listput(eq,matker(M~)~*G);if(#jj,fix[jj[1]]=G));
  for(v=1,nf.sign[1],isolationIndex=v;L=realrow(es,ds,tang);
    F=[L[2]*L[3]<0,L[1]*L[3]<0]~;G=Mat(apply(z->sg(z)<0,U))*Mod(1,2);
    listput(eq,matker(realimage()~)~*matconcat([matconcat([G,0;0,G]),F])));
  my(A=matker(matconcat(Col(eq))),last=2*n+1);
  my(k=select(i->A[last,i]!=0,[1..matsize(A)[2]])[1],base=A[,k],diff=A*matker(Mat(A[last,])));
  my(liftbits=vector(2,i,logs[i]*matinverseimage(LB[i],fix[i]*matconcat([base,diff]))));
  my(allbits=matconcat(Col(liftbits)),offset=allbits[,1]+concat(refbits),dirs=allbits[,2..matsize(allbits)[2]]);
  my(test=[1,0,0,1,0]*Mod(1,2));
  if(test*offset!=1 || test*dirs!=0 || test*ray!=0,error("four-cover obstruction"));
  print("The sum of the ramified first and unramified second coordinates is 1 for every primitive lift and 0 on the center direction.");
};
affine();
