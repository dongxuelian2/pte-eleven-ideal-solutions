"""Dependency-free exact certificate; coefficients in increasing order."""
def add(a,b,scale=1):
 out=[0]*max(len(a),len(b))
 for i,c in enumerate(a):out[i]+=c
 for i,c in enumerate(b):out[i]+=scale*c
 while len(out)>1 and out[-1]==0:out.pop()
 return out
def mul(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,c in enumerate(a):
  for j,d in enumerate(b):out[i+j]+=c*d
 return out
A=[162,-159,118,-36,8]
B=[54,-27,54,-20,8]
fplus=mul(mul([0,1],[3,1,2]),mul(A,A))
fminus=mul(mul([-1,2],[8,-3,1]),mul(B,B))
assert add(fplus,fminus,-1)==[23328]
# D=Dn/4; S=Sn/256; T=Tn/64. Clear the common denominator.
Dn=mul(mul([0,1],[-1,2]),mul([8,-3,1],[3,1,2]))
Sn=add(fplus,fminus)
Tn=mul(A,B)
assert add(mul(Sn,Sn),mul(Dn,mul(Tn,Tn)),-4)==[729**2*1024]
assert len(Sn)-1==11 and len(Tn)-1==8 and len(Dn)-1==6
print('PASS: Fplus-Fminus=23328 and degree-eleven sextic Pell identity, using integer arithmetic only.')

from fractions import Fraction
def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def rem(a,b):
 a=list(map(Fraction,a));b=list(map(Fraction,b))
 while len(a)>=len(b) and a!=[0]:
  k=len(a)-len(b);c=a[-1]/b[-1]
  for j,v in enumerate(b):a[k+j]-=c*v
  trim(a)
 return a
d=[i*Dn[i] for i in range(1,len(Dn))]
a,b=Dn,d
while b!=[0]:a,b=b,rem(a,b)
assert len(a)==1 and a[0]!=0
assert Dn[-1]==4 and Sn[-1]==256 and Tn[-1]==64
assert (-3)**2-4*8==-23 and 1-4*2*3==-23
print('PASS: D squarefree; D,S,T monic; both quadratic discriminants equal -23.')
