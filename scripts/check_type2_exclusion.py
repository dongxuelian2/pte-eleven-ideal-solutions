"""Independent exact quotient maps, point list and collisions for PROOF 5.45.
The elementary two-isogeny residue tables prove the unconditional rank bound.
"""
from fractions import Fraction as F
import sympy as s
v,y,X,Y=s.symbols('v y X Y')
f=-240*v**4+760*v**2+81
H=lambda a,b:8*a*a-40*a*b-58*a+80*b*b+50*b+95
a=(y+20*v*v+29)/8
assert s.expand(8*H(a,v*v)-(y*y-f))==0
ex=F(760,81)+(F(2,9)*y+2)/v**2
ey=2*ex/v
E=lambda xx,yy:yy**2-(xx**3-F(1520,81)*xx**2+F(655360,6561)*xx)
assert s.rem(s.together(E(ex,ey)).as_numer_denom()[0],y*y-f,y)==0
iv=2*X/Y
iy=F(9,2)*(X-F(760,81))*iv**2-9
assert s.factor(ex.subs({v:iv,y:iy})-X)==0
assert s.factor(ey.subs({v:iv,y:iy})-Y)==0
points=[(F(0),F(9)),(F(0),F(-9))]
points += [(vv,yy) for vv in [F(1,2),F(-1,2)] for yy in [F(16),F(-16)]]
points += [(vv,yy) for vv in [F(3,2),F(-3,2)] for yy in [F(24),F(-24)]]
assert len(set(points))==10
for vv,yy in points:
 assert yy*yy==f.subs(v,vv)
 aa=F(yy+20*vv*vv+29,8)
 if vv==0:
  assert aa in [F(19,4),F(5,2)]
  assert not s.sqrt(s.Rational(aa)).is_Rational
 else:
  bchain={vv+dd for dd in [-1,0,1]}|{-vv+dd for dd in [-1,0,1]}
  assert bchain&{F(-1,2),F(1,2)}
# Independent torsion upper bound on the integral minimal model.
b2,b4,b6,b8=1,-56,1088,272-784
disc=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
assert disc==-31104000
for p in [7,11]:
 assert disc%p
 assert 1+sum((yy*yy+xx*yy-xx**3+28*xx-272)%p==0
              for xx in range(p) for yy in range(p))==10
print('PASS: type2 quotient maps, all ten points, good reductions and exclusions')
# Elementary two-isogeny descent, independent of the Magma rank bound.
rows=[(-95,2560,2,3,3,{2}),(-95,2560,5,3,3,{2}),
 (190,-1215,-1,2,16,{7,12,15}),
 (190,-1215,3,2,16,{3,11,12}),
 (190,-1215,-5,2,16,{3,11,12}),
 (190,-1215,15,2,16,{7,12,15}),
 (190,-1215,-3,2,128,{5,21,29,69,80,93,117,125}),
 (190,-1215,5,2,128,{5,13,21,77,80,85,101,109})]
for A,B,d,p,n,expected in rows:
 charts=[(u,1) for u in range(n)]+[(1,vv) for vv in range(0,n,p)]
 residues={(d*u**4+A*u*u*vv*vv+(B//d)*vv**4)%n for u,vv in charts}
 assert residues==expected
 assert not residues&{w*w%n for w in range(n)}
# Scaling the quotient E to the integral two-isogeny model.
XX,YY=F(81,16)*ex,F(729,64)*ey
res=YY**2-(XX**3-95*XX**2+2560*XX)
assert s.rem(s.together(res).as_numer_denom()[0],y*y-f,y)==0
integral_points={(0,0)}
for vv,yy in points:
 if not vv:continue
 pp=(s.factor(XX.subs({v:vv,y:yy})),s.factor(YY.subs({v:vv,y:yy})))
 integral_points.add(pp)
assert integral_points=={(0,0),(16,144),(16,-144),(40,120),(40,-120),
 (64,192),(64,-192),(160,1440),(160,-1440)}
for xx,yy in integral_points:assert yy*yy==xx**3-95*xx**2+2560*xx
# Integral E has discriminant 16*B^2*(A^2-4B), and reduction counts 10.
assert 16*2560**2*((-95)**2-4*2560)==-(2**22)*(3**5)*(5**3)
for p in [7,11]:
 assert 1+sum((yy*yy-xx**3+95*xx**2-2560*xx)%p==0
              for xx in range(p) for yy in range(p))==10
print('PASS: elementary two-isogeny rank-zero proof and ten integral E points')

# The exact integral-coordinate maps printed in the proof.
iv_int=F(9,2)*X/Y
iy_int=F(8,9)*(X-F(95,2))*iv_int**2-9
assert s.factor(XX.subs({v:iv_int,y:iy_int})-X)==0
assert s.factor(YY.subs({v:iv_int,y:iy_int})-Y)==0
print('PASS: integral-coordinate forward and inverse maps')
