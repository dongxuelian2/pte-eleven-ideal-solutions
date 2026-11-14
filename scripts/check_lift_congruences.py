"""Uniform necessary index congruences, not a complete PTE classification."""
from fractions import Fraction as Q
from math import gcd,lcm
A,B=-560,72900
Delta=16*B*B*(A*A-4*B)
def add(S,T,p):
 if S is None:return T
 if T is None:return S
 x,y=S;xx,yy=T
 if x==xx and (y+yy)%p==0:return None
 slope=((yy-y)*pow((xx-x)%p,-1,p))%p if x!=xx else ((3*x*x+2*A*x+B)*pow(2*y,-1,p))%p
 z=(slope*slope-A-x-xx)%p
 return z,(slope*(x-z)-y)%p

def possible(S,p):
 # Keep every point where this affine reduction of the lift is undefined.
 if S is None or S[0]==0:return True
 X,Y=S;z=Y*pow(2*X,-1,p)%p;t=(z+5)*pow(3,-1,p)%p
 d=(3*t*t-10*t+10)%p
 if d==0:return True
 a=2*(6*t*t-5*t+5)*pow(d,-1,p)%p
 b=(3*t*t+26*t+10)*pow(4*d,-1,p)%p
 squares={j*j%p for j in range(p)}
 return a in squares and b in squares

local=[]
for p,order in [(13,6),(37,40),(59,16),(113,120)]:
 assert Delta%p!=0 and p>3
 P=(180%p,900%p);S=None;cycle=[]
 for _ in range(order):
  cycle.append(S);S=add(S,P,p)
 assert S is None and len(set(cycle))==order
 assert all(T is None or (T[1]*T[1]-T[0]*(T[0]*T[0]+A*T[0]+B))%p==0 for T in cycle)
 sets=[{n for n,T in enumerate(cycle) if possible(add(T,(0,0),p) if e else T,p)} for e in (0,1)]
 local.append((order,sets))
assert lcm(*(order for order,_ in local))==240
survivors=[{n for n in range(240) if all(n%order in sets[e] for order,sets in local)} for e in (0,1)]
assert survivors==[{0,1,120},{0,120,239}]
# The only allowed |n|<120 are the two torsion points and P,-P+T.
assert [{n for n in range(-119,120) if n%240 in survivors[e]} for e in (0,1)]==[{0,1},{-1,0}]
# At t=5/2, the recovered B length-three chain meets central A at 3/2.
t=Q(5,2);d=3*t*t-10*t+10
assert (2*(6*t*t-5*t+5)/d,(3*t*t+26*t+10)/(4*d),(3*t*t-10*t+100)/(4*d))==(16,Q(25,4),Q(25,4))
assert Q(5,2)-1==Q(3,2)
# At t=infinity, (a,b,c)=(4,1/4,1/4), also a central collision.
assert Q(1,2) in {Q(-3,2),Q(-1,2),Q(1,2),Q(3,2)}

def rational_add(S,T):
 if S is None:return T
 if T is None:return S
 x,y=S;xx,yy=T
 if x==xx and y==-yy:return None
 slope=(yy-y)/(xx-x) if x!=xx else (3*x*x+2*A*x+B)/(2*y)
 z=slope*slope-A-x-xx
 return z,slope*(x-z)-y
P=(Q(180),Q(900));T=(Q(0),Q(0))
assert rational_add((P[0],-P[1]),T)==(Q(405),Q(2025))
for R in [P,(Q(405),Q(2025))]:assert (R[1]/(2*R[0])+5)/3==Q(5,2)
PP=rational_add(P,P);PPPP=rational_add(PP,PP)
u=PPPP[0]/270
assert max(abs(u.numerator),u.denominator)==77763
# Exact logarithmic inequality proving scaled elliptic height >10^3298.
assert 120**2//16==900
assert 77763**2700*11**900 > 10**9894*51840**900*440
print('Uniform residues mod 240, all |n|<120 collisions, and elliptic height lower bound verified.')
print('This is not a complete simultaneous-lift classification or an original PTE height bound.')
