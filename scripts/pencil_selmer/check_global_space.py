import json,re
import sympy as s
from pathlib import Path
p=Path(__file__).parent
o=p.parents[1]/'build'/'pencil_selmer';o.mkdir(parents=True,exist_ok=True)
raw=dict(success=True,truncated=False,stdout=(p/'selected.txt').read_text())
assert raw['success'] and not raw['truncated']
txt=raw['stdout'].replace('\\\n','')
gens=[];vals=[]
for i,g,v in re.findall(r'GEN\s+(\d+)\s+\[([^]]+)\]\s+VAL\s+\[([^]]+)\]',txt):
 assert int(i)==len(gens)+1
 gens.append([s.Rational(z.strip()) for z in g.split(',')]);vals.append([int(z.strip()) for z in v.split(',')])
assert len(gens)==14
x=s.symbols('x');g=x**5-660*x**3+70983*x-102060
intervals=[ab for ab,m in s.Poly(g,x).intervals(eps=s.Rational(1,10**150))]
assert len(intervals)==5

def valint(c,ab):
 lo=hi=s.Rational(0)
 for z in reversed(c):
  v=[lo*ab[0],lo*ab[1],hi*ab[0],hi*ab[1]];lo=min(v)+z;hi=max(v)+z
 return lo,hi
signs=[]
for c in gens:
 row=[]
 for ab in intervals:
  a,b=valint(c,ab);assert a>0 or b<0
  row.append((1 if b<0 else 0))
 signs.append(row)
# Necessary real image in increasing a order: span 10100,01010.
rows=[]
for i in range(3):
 rows.append(sum((((v[0]^v[1]^v[2]) if i==0 else (v[5]^v[6]) if i==1 else (v[7]^v[8]))<<j) for j,v in enumerate(vals)))
for i in range(3):
 rows.append(sum(((v[4] if i==0 else v[0]^v[2] if i==1 else v[1]^v[3])<<j) for j,v in enumerate(signs)))
def rank(rr):
 b={}
 for z in rr:
  while z:
   i=z.bit_length()-1
   if i in b:z^=b[i]
   else:b[i]=z;break
 return len(b)
allowed=[z for z in range(1<<14) if all((z&r).bit_count()%2==0 for r in rows)]
assert len(allowed)==2**(14-rank(rows))
print('Matrix rank',rank(rows),'candidate dimension',14-rank(rows))
print('Sign vectors',signs)
print('Constraint rows',rows)
out={'gens':[[str(x)for x in row]for row in gens],'vals':vals,'signs':signs,'constraints':rows,'allowed_masks':allowed}
json.dump(out,open(o/'global_space.json','w'))
