import json,re
from pathlib import Path
p=Path(__file__).parent
o=p.parents[1]/'build'/'pencil_selmer';o.mkdir(parents=True,exist_ok=True)
d=json.load(open(o/'space_after_5_7.json'));raw=dict(success=True,truncated=False,stdout=(p/'local_odd.txt').read_text());assert raw['success'] and not raw['truncated']
a=d['allowed_masks']
def pack(row):return sum(x<<i for i,x in enumerate(row))
def rows(s):return [[int(x.strip())for x in z.split(',')]for z in re.findall(r'\[([^]]+)\]',s)]
for block in raw['stdout'].split('PRIME ')[1:]:
 prime=int(block.splitlines()[0]);vstr,tstr=block.split('VECTORS')[1].split('TORSION');v=rows(vstr);ts=rows(tstr)
 cols=list(map(pack,v[:14]));loc=list(map(pack,ts+v[14:]));span={0}
 for z in loc:span|={t^z for t in list(span)}
 assert len(span)=={3:8,773:4,3087467:8}[prime]
 def img(z):
  r=0
  for j,c in enumerate(cols):
   if z>>j&1:r^=c
  return r
 a=[z for z in a if img(z) in span]
 print('After prime',prime,'candidates',len(a),'dimension',len(a).bit_length()-1)
 assert len(a)&(len(a)-1)==0
 d[f'local_{prime}']={'columns':cols,'span':sorted(span)}
d['allowed_masks']=a
json.dump(d,open(o/'space_after_odd.json','w'))
