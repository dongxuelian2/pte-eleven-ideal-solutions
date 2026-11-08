from math import gcd
A,B=-560,72900
Disc=16*B*B*(A*A-4*B)
def add(S,T,p):
 if S is None:return T
 if T is None:return S
 x,y=S;xx,yy=T
 if x==xx and (y+yy)%p==0:return None
 v=((yy-y)*pow((xx-x)%p,-1,p))%p if x!=xx else ((3*x*x+2*A*x+B)*pow(2*y,-1,p))%p
 z=(v*v-A-x-xx)%p
 return z,(v*(x-z)-y)%p
M=3120
sur=[{n for n in range(M) if n%240 in [0,1,120]},{n for n in range(M) if n%240 in [0,120,239]}]
for p in [61,101,107]:
 assert Disc%p!=0
 P=(180%p,900%p);S=None;cycle=[]
 while True:
  cycle.append(S);S=add(S,P,p)
  if S is None:break
 order=len(cycle)
 assert order=={61:26,101:13,107:104}[p]
 assert len(set(cycle))==order
 assert all(S is None or (S[1]**2-S[0]*(S[0]**2+A*S[0]+B))%p==0 for S in cycle)
 assert order%13==0
 squares={j*j%p for j in range(p)}
 def possible(S):
  if S is None or S[0]==0:return True
  X,Y=S;z=Y*pow(2*X,-1,p)%p
  qs=[(z*z+5)%p,(4*z*z+30*z+80)%p,(z*z+36*z+185)%p,(z*z+275)%p]
  if any(q==0 for q in qs):return True
  return len({q in squares for q in qs})==1
 g=gcd(M,order)
 allowed=[set(),set()]
 for n,S in enumerate(cycle):
  for e in (0,1):
   SS=add(S,(0,0),p) if e else S
   if possible(SS):allowed[e].add(n%g)
 old=list(map(len,sur));sur=[{n for n in sur[e] if n%g in allowed[e]} for e in (0,1)]
 if old!=list(map(len,sur)):print(p,order,g,old,list(map(len,sur)),[sorted({n%13 for n in z}) for z in sur],flush=True)
print('FINAL',list(map(len,sur)),[sorted(z) for z in sur])
assert sur==[{0,1,1560},{0,1560,3119}]
assert {n%13 for n in sur[0] if n%2==0}=={0}
assert {n%13 for n in sur[1] if n%2==1}=={12}
print('PASS: sixth-mode B1/B5 elliptic index residues must be0/12 modulo13.')
# Independent elliptic formal logarithm calibration, using rational group law.
from fractions import Fraction as F
Pq=(F(180),F(900))
def qadd(S,T):
 if S is None:return T
 if T is None:return S
 x,y=S;xx,yy=T
 if x==xx and y==-yy:return None
 slope=(yy-y)/(xx-x) if x!=xx else (3*x*x+2*A*x+B)/(2*y)
 z=slope*slope-A-x-xx
 return z,slope*(x-z)-y
six=None
for _ in range(6):six=qadd(six,Pq)
t6=-six[0]/six[1]
t6mod=int(t6.numerator)*pow(int(t6.denominator),-1,169)%169
assert t6mod==78
# E integral at13, t(6P) in13Z13. Integral invariant differential has
# leading coefficient1; its formal logarithm is t+O(13²) on this disc.
assert (t6mod//13)*pow(6,-1,13)%13==1
print('PASS: independent formal log(P)/13=1 mod13 from t(6P)=78 mod169.')
