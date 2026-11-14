import json,re
from pathlib import Path
p=Path(__file__).parent
o=p.parents[1]/'build'/'pencil_selmer';o.mkdir(parents=True,exist_ok=True)
d=json.load(open(o/'global_space.json'));raw=dict(success=True,truncated=False,stdout=(p/'local_5_7.txt').read_text());assert raw['success'] and not raw['truncated']
a=d['allowed_masks']
for block in raw['stdout'].split('PRIME ')[1:]:
 prime=int(block.splitlines()[0]);v=[[int(x.strip())for x in z.split(',')] for z in re.findall(r'\[([^]]+)\]',block)];assert len(v)==15
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
