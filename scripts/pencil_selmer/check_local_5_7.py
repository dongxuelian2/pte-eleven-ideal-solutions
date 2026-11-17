import json,re
from pathlib import Path
p=Path(__file__).parent
o=p.parents[1]/'build'/'pencil_selmer';o.mkdir(parents=True,exist_ok=True)
d=json.load(open(o/'global_space.json'))
blocks=(p/'local_5_7.txt').read_text().split('PRIME ')[1:]
primes=[int(block.splitlines()[0]) for block in blocks]
assert len(primes)==2 and set(primes)=={5,7}, primes
a=d['allowed_masks']
for block in blocks:
 prime=int(block.splitlines()[0]);v=[[int(x.strip())for x in z.split(',')] for z in re.findall(r'\[([^]]+)\]',block)];assert len(v)==15
 assert all(len(row)==4 and set(row)<={0,1} for row in v)
 def pack(row):return sum(x<<i for i,x in enumerate(row))
 cols=list(map(pack,v[:14]));beta=pack(v[-1]);assert beta
 def img(z):
  r=0
  for j,c in enumerate(cols):
   if z>>j&1:r^=c
  return r
 a=[z for z in a if img(z) in (0,beta)]
 print('After prime',prime,'candidates',len(a),'dimension',len(a).bit_length()-1)
 assert len(a)&(len(a)-1)==0
 d[f'local_{prime}']={'columns':cols,'beta':beta}
d['allowed_masks']=a
json.dump(d,open(o/'space_after_5_7.json','w'))
