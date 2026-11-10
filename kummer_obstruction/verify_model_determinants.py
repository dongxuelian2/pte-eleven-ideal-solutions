import json,re,pathlib,sympy as s
p=pathlib.Path(__file__).parent
models=json.loads((p/'all_seven_models.json').read_text());z=s.symbols('z')
g=4185074574376548950859684249600*z**6+41672499468952215879696384*z**5+31450141261690198740*z**4-42408244947519*z**3+1600830*z**2+z
f6=s.Poly(g,z).LC();G=s.zeros(6)
for i in range(6):G[i,5-i]=1
for obj in models:
 la,*a=map(int,re.findall(r'-?\d+',obj['model']));assert len(a)==21
 H=s.zeros(6);k=0
 for i in range(6):
  for j in range(i,6):
   H[i,j]=H[j,i]=a[k]*(2 if i==j else 1);k+=1
 for x in range(7):
  assert (la*x*G-H).det()==-s.Rational(la**6,f6)*g.subs(z,x), (obj['vector'],x)
 print('det polynomial verified at 7 distinct integers; degree <= 6:',obj['vector'])
print('ALL_SEVEN_DETERMINANT_IDENTITIES_VERIFIED_INDEPENDENTLY')
