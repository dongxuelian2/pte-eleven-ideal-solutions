"""Check the global certificate's fixed inputs, rank argument, and receipts.
This does not reimplement Magma's ideal arithmetic or FactorBasisVerify.
"""
import ast
import json
import re
from pathlib import Path
import sympy as s
p=Path(__file__).resolve().parent
d=json.loads((p.parents[1]/'build/pencil_selmer/global_space.json').read_text())
def rank(rows):
    pivots={}
    for row in rows:
        v=sum(c<<i for i,c in enumerate(row))
        while v:
            j=v.bit_length()-1
            if j in pivots: v^=pivots[j]
            else:
                pivots[j]=v
                break
    return len(pivots)
assert all(not any(row) for row in d['vals'][:5])
assert rank(d['signs'][:5]) == 5
assert rank(d['vals']) == 9
assert rank([v+sig for v,sig in zip(d['vals'],d['signs'])]) == 14
code=(p/'certify_global_basis.m').read_text()
expected=ast.literal_eval(re.search(r'Expected := (.*?);',code,re.S).group(1))
assert expected==d['vals']
x=s.symbols('x');g=x**5-660*x**3+70983*x-102060
terms=re.search(r'G:=\[(.*?)\];',code).group(1).split(',K!')
facts=re.search(r'Facts := \[\* (.*?) \*\];',code,re.S).group(1)
facts=ast.literal_eval('['+facts.replace('<','(').replace('>',')')+']')
assert len(terms)==len(facts)==14
for row,term,ff in zip(d['gens'],terms,facts):
    gen=sum(s.Rational(c)*x**i for i,c in enumerate(row))
    assert s.expand(s.sympify(term.removeprefix('K!').replace('^','**'),locals={'a':x})-gen)==0
    assert all(s.isprime(q) and e>0 for q,e in ff)
    assert len(set(q for q,e in ff))==len(ff)
    assert abs(s.resultant(g,gen,x))==s.prod(q**e for q,e in ff)
# Successful outputs cover exactly the declared intervals. No failed or
# incomplete response was admitted when preparing this mathematical receipt.
text=(p/'global_basis.txt').read_text()
intervals=[tuple(map(int,z)) for z in re.findall(r'VERIFIED_INTERVAL\s+(\d+)\s+(\d+)',text)]
expected=[(318,50000),(50000,1000000),(1000000,3000000),(3000000,6000000),
          (6000000,10000000),(10000000,14000000),(14000000,17507890)]
assert intervals==expected
assert all(b==c for (a,b),(c,e) in zip(intervals,intervals[1:]))
assert re.findall(r'EXACT_RELATIONS_VERIFIED\s+(\d+)',text)==['80']+['85']*6
assert re.findall(r'LOW_PRIME_COVERAGE\s+(\d+)\s+(\d+)',text)==[('63','318')]*7
assert re.findall(r'MINKOWSKI_BOUND\s+(\d+)',text)==['17507890']*7
assert list(map(int,re.findall(r'EXACT_GLOBAL_IDEAL_PARITY\s+(\d+)\s+VERIFIED',text)))==list(range(1,15))
assert 'ALL_14_GLOBAL_PARITIES_CERTIFIED_WITHOUT_CLASSGROUP_ASSUMPTIONS' in text
# Minkowski's totally-real degree-five bound is (5!/5^5)*sqrt(D).
D=207876406776527712; B=17507890
assert B**2*3125**2 <= 120**2*D < (B+1)**2*3125**2
print('Exact norm factorizations and basis bindings passed.')
print('Sign rank 5, ideal-parity rank 9, combined rank 14.')
print('Successful interval receipts cover the entire Minkowski bound.')
