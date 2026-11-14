import json
from fractions import Fraction
from pathlib import Path
P=Path(__file__).parent
o=P.parents[1]/'build'/'pencil_selmer';o.mkdir(parents=True,exist_ok=True)
D=json.load(open(o/'global_space.json'))
g=[-102060,70983,0,-660,0,1]
th=[Fraction(-750870,2),Fraction(-255150,2),Fraction(70983,2),0,Fraction(-165,2)]
def trim(a):
 while a and not a[-1]:a.pop()
 return a
def add(a,b):
 c=[Fraction(0)]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return trim(c)
def mul(a,b):
 c=[Fraction(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def modpoly(a,h,m):
 a=[x%m for x in a]
 while len(a)>=len(h):
  c=a[-1]; k=len(a)-len(h)
  for j,z in enumerate(h):a[k+j]=(a[k+j]-c*z)%m
  trim(a)
 return a

def powpoly(a,n,h,p):
 r=[1]
 while n:
  if n&1:r=modpoly(mul(r,a),h,p)
  a=modpoly(mul(a,a),h,p);n//=2
 return r

def ev(c,r,m):
 a=0
 for z in reversed(c):a=(a*r+z)%m
 return a

def vp(n,p):
 if n==0:return 10**9
 v=0
 while n%p==0:n//=p;v+=1
 return v

def char57(c,p,N=32):
 m=p**N;r=0
 # Newton in Z/p^N; g'(0) is a unit for p=5,7.
 dg=[i*g[i] for i in range(1,len(g))]
 for _ in range(N.bit_length()+1):r=(r-ev(g,r,m)*pow(ev(dg,r,m),-1,m))%m
 assert ev(g,r,m)==0
 h=[0]*5;h[4]=1
 for j in range(3,-1,-1):h[j]=(g[j+1]+r*h[j+1])%m
 assert (-r*h[0]-g[0])%m==0
 cc=[int(Fraction(z).numerator)*pow(int(Fraction(z).denominator),-1,m)%m for z in c]
 out=[]
 for hc in [[-r%m,1],h]:
  rem=modpoly(cc,hc,m);v=min([vp(int(a),p) for a in rem],default=10**9);assert v<N
  res=[int(a)//p**v%p for a in rem]; hp=[int(a)%p for a in hc];f=len(hp)-1
  leg=powpoly(res,(p**f-1)//2,hp,p)
  assert leg in [[1],[p-1]]
  out += [v%2,0 if leg==[1] else 1]
 return out

out={}
for p,b,c in [(5,85,23),(7,22,78271)]:
 beta=add(add(mul(th,th),[b*z for z in th]),[c])
 rows=[char57(cs,p) for cs in D['gens']]; br=char57(beta,p)
 assert rows==[char57(cs,p,40) for cs in D['gens']] and br==char57(beta,p,40)
 assert br!=[0]*4
 allowed=[]
 for mask in D['allowed_masks']:
  a=[sum(rows[j][k] for j in range(14) if mask>>j&1)%2 for k in range(4)]
  if a in [[0]*4,br]:allowed.append(mask)
 out[str(p)]={'rows':rows,'beta':br,'allowed_masks':allowed}
 print('p',p,'beta',br,'rows',rows,'survivors',len(allowed))
both=sorted(set(out['5']['allowed_masks'])&set(out['7']['allowed_masks']))
print('both survivors',len(both));out['both']=both
json.dump(out,open(o/'audit_local57.json','w'))
