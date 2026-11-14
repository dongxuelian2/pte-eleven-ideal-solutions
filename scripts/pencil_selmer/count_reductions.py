from math import gcd
coeff=[15173754195676653158400,-27536385315056910336,3787455316296960,930770808176,6403320,-729]
D=2**77*3**33*5**4*7**4*61**3*773**2*3087467**2*35498020282877**3
def count(p):
 def fv(x):
  y=0
  for c in coeff[::-1]: y=(y*x+c)%p
  return y
 def chi(x): return 0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1)
 n1=1+sum(1+chi(fv(x)) for x in range(p))
 d=next(d for d in range(2,p) if chi(d)==-1)
 def mul(x,y): return ((x[0]*y[0]+d*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0])%p)
 def power(x,k):
  y=(1,0)
  while k:
   if k&1:y=mul(x,y)
   x=mul(x,x); k//=2
  return y
 n2=1
 for a in range(p):
  for b in range(p):
   y=(0,0)
   for c in coeff[::-1]:
    y=mul(y,(a,b)); y=((y[0]+c)%p,y[1])
   q=power(y,(p*p-1)//2)
   n2+=1+(0 if y==(0,0) else (1 if q==(1,0) else -1))
 return n1,n2,(n1*n1+n2)//2-p
G=0
for p in [11,13,17,19,23,29,31]:
 assert D%p
 a,b,c=count(p);G=gcd(G,c);print(p,a,b,c,'gcd',G)
