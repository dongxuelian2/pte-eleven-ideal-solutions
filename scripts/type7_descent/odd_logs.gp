\\ Formal parameters of multiples, evaluated in the ring modulo p^2 by
\\ the standard division polynomials of the short Weierstrass equation.
divpol(n)={
  my(v,m=n\2);
  if(mapisdefined(cache,n,&v),return(v));
  if(n%2,v=divpol(m+2)*divpol(m)^3-divpol(m-1)*divpol(m+1)^3,
    v=divpol(m)/(2*ordinate)*(divpol(m+2)*divpol(m-1)^2-divpol(m-2)*divpol(m+1)^2));
  mapput(cache,n,v);v;
};
formal(n,X,Y,a,b,c)={
  my(A=b-a^2/3,B=c-a*b/3+2*a^3/27,u=X+a/3);
  if(Y^2!=u^3+A*u+B,error("point off curve"));
  ordinate=Y;cache=Map();
  my(initial=[0,1,2*Y,3*u^4+6*A*u^2+12*B*u-A^2,
    4*Y*(u^6+5*A*u^4+20*B*u^3-5*A^2*u^2-4*A*B*u-8*B^2-A^3)]);
  for(i=0,4,mapput(cache,i,initial[i+1]));
  my(v=divpol(n),phi=u*v^2-divpol(n-1)*divpol(n+1));
  my(omega=(divpol(n+2)*divpol(n-1)^2-divpol(n-2)*divpol(n+1)^2)/(4*Y));
  -phi*v/omega;
};
basis3=[[225*y^2-630*y-90,-6750*y^2+17550*y+5400],[225*y^2-180*y+1260,-6750*y^2+5400*y-37800]];
basis10=[[-185*y^2-50*y+1340,-11150*y^2+1000*y+85250],[100*y^2-80*y-844,1900*y^2-1520*y-16036]];
inert(d,t,p,pts,det)={
  my(f=y^3-8*y-2);if(!polisirreducible(Mod(1,p)*f),error("inert prime"));
  my(ff=ffgen(Mod(1,p)*f,'u),l=Mod(Mod(1,p^2)*y,f),e=Mod((-5*y^2+4*y+53)/27,f));
  my(k=135*d*(t-e),a=k*(e-79/27),b=k^2*(e^2-79/27*e+57/27),A=subst(lift(a),y,l),B=subst(lift(b),y,l));
  my(E=ellinit([0,subst(lift(a),y,ff),0,subst(lift(b),y,ff),0]),n=ellcard(E),rows=List());
  for(i=1,#pts,
    my(X=Mod(pts[i][1],f),Y=Mod(pts[i][2],f));
    if(Y^2!=X*(X^2+a*X+b),error("global generator"));
    my(z=liftall(formal(n,subst(lift(X),y,l),subst(lift(Y),y,l),A,B,0)));
    if(content(z)%p,error("formal reduction"));
    listput(rows,vector(3,j,polcoef(z,j-1)/p%p)));
  my(X=k*t,Y=nfroots(nfinit(f),x^2-X*(X^2+a*X+b))[1],dr=liftall(subst(lift(k/(2*Y)),y,l)));
  listput(rows,vector(3,j,polcoef(dr,j-1)%p));
  my(v=matdet(Mat(vector(3,j,Col(rows[j]))))%p);
  if(n%p==0 || (v!=det && v!=p-det),error("inert logarithm determinant"));
};
base=ellinit([0,1,0,-13533,-323937]);
if(elllocalred(base,2)!=[2,-4,[1,0,0,0],3],error("base reduction at 2"));
inert(3,3,61,basis3,31);
inert(10,2/5,389,basis10,216);
splitpacket()={
  my(p=83,f=y^3-5*y^2-27*y+81,rt=polrootsmod(f,p),logs=List(),slopes=List(),z0=1/3);
  if(#rt!=3,error("splitting prime"));
  for(i=1,3,my(h=Mod(lift(rt[i]),p^2));h-=subst(f,y,h)/subst(deriv(f,y),y,h);
    my(r=h/9,k=(z0-r)*(9*z0^3+5*z0^2-3*z0-1),m=k*(9*r^3+5*r^2-3*r-1));
    my(a=k*(27*r^2+10*r-3),b=k*m*(27*r+5),c=9*k*m^2,X=m/(z0-r),Y=m*k/(z0-r)^2);
    my(n=ellcard(ellinit(liftall([0,a,0,b,c]),p)),z=lift(formal(n,X,Y,a,b,c)));
    if(z%p,error("formal reduction"));
    listput(logs,z/p%p);listput(slopes,lift(-m/(z0-r)^2/(2*Y))*n%p));
  my(minors=vector(3,i,(logs[i]*slopes[i%3+1]-logs[i%3+1]*slopes[i])%p));
  if(minors!=[67,70,22],error("split logarithm minors"));
};
splitpacket();
print("Verified the two inert determinants and the three split minors.");
