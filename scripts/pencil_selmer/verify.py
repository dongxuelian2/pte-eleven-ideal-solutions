"""Exact audit of the fixed pencil and recorded Magma 2-descent certificates."""
import json
import re
import runpy
from pathlib import Path
import sympy as s

p = Path(__file__).resolve().parent
o = p.parents[1] / 'build' / 'pencil_selmer'
for script in ['check_pencil', 'count_reductions', 'check_global_space',
               'check_local_5_7', 'check_local_odd', 'audit_local57', 'audit_witnesses']:
    print('\nCHECK', script)
    result = runpy.run_path(str(p / (script + '.py')))
    if script == 'count_reductions':
        assert result['count'](11) == (3, 123, 55)
        assert result['count'](13) == (12, 206, 162)

d = json.loads((o / 'space_after_odd.json').read_text())
audit = json.loads((o / 'audit_local57.json').read_text())
assert audit['both'] == json.loads((o / 'space_after_5_7.json').read_text())['allowed_masks']
for prime in [5, 7]:
    rows = audit[str(prime)]['rows']
    assert [sum(v << i for i, v in enumerate(row)) for row in rows] == d[f'local_{prime}']['columns']
    assert sum(v << i for i, v in enumerate(audit[str(prime)]['beta'])) == d[f'local_{prime}']['beta']
assert len(d['allowed_masks']) == 16
# Bind the dyadic input masks and all local scripts to the exported basis.
x, t = s.symbols('x t')
g = x**5 - 660*x**3 + 70983*x - 102060
th = (-165*x**4 + 70983*x**2 - 255150*x - 750870)/2
gens = [sum(s.Rational(c)*x**i for i, c in enumerate(row)) for row in d['gens']]
for name in ['local_5_7', 'local_odd', 'local_two', 'local_remaining']:
    text = (p / (name + '.m')).read_text()
    terms = re.search(r'G:=\[(.*?)\];', text).group(1).split(',K!')
    for term, gen in zip(terms, gens):
        term = term.removeprefix('K!')
        assert s.expand(s.sympify(term.replace('^', '**'), locals={'a': x}) - gen) == 0
    assert len(terms) == 14
masks = [int(v) for v in re.search(r'Masks:=\[([^]]+)\]', (p/'local_two.m').read_text()).group(1).split(',')]
assert masks == d['allowed_masks']
raw = (p / 'local_two.txt').read_text()
assert 'LOCAL_CLASS_COUNT 16' in raw
survivors = [int(v) for v in re.search(r'SURVIVORS\s*\[([^]]+)\]', raw).group(1).split(',')]
assert survivors == [0, 4, 3169, 3173, 12298, 12302, 15467, 15471]
assert set(survivors) <= set(masks)
assert sorted({(4 if i&1 else 0) ^ (3169 if i&2 else 0) ^ (12298 if i&4 else 0) for i in range(8)}) == survivors

def rows(text):
    return [[int(v) for v in row.split(',')] for row in re.findall(r'\[([^]]+)\]', text)]
def pack(row):
    assert all(v in [0, 1] for v in row)
    return sum(v << i for i, v in enumerate(row))
for block in (p/'local_remaining.txt').read_text().split('PRIME ')[1:]:
    prime = int(block.splitlines()[0])
    vtext, ttext = block.split('VECTORS')[1].split('TORSION')
    v, torsion = rows(vtext), rows(ttext)
    assert len(v) == 15 and len(torsion) == 3
    # All local-image generators have even valuations at these omitted primes.
    assert all(all(z == 0 for z in row[::2]) for row in torsion + v[14:])
    span = {0}
    for z in map(pack, torsion + v[14:]):
        span |= {a ^ z for a in list(span)}
    assert len(span) == 4
    for mask in survivors:
        image = 0
        for j, row in enumerate(v[:14]):
            if mask >> j & 1:
                image ^= pack(row)
        assert image in span
    print('Full local image at', prime, 'has dimension 2; all eight classes survive.')
F = x**5+3478005*x**4-38594700111069*x**3-14431460273343469305*x**2+47793145649188822875249984*x+22058431254036587705175782036640
f = -729*t**5+6403320*t**4+930770808176*t**3+3787455316296960*t**2-27536385315056910336*t+15173754195676653158400
assert s.expand(F.subs(x, -729*t/4-375435) - s.Rational(531441,32)**2*f) == 0
assert s.rem(F.subs(x, th), g, x) == 0
assert s.Poly(g, x).is_irreducible
assert s.discriminant(f, t) == 2**77*3**33*5**4*7**4*61**3*773**2*3087467**2*35498020282877**3
rcoeff = [s.Rational(v) for v in re.search(r'\[([^]]+)\]', (p/'known_root.txt').read_text()).group(1).split(',')]
r = sum(c*x**i for i, c in enumerate(rcoeff))
assert s.rem(r*r - (-375435-th)*gens[2], g, x) == 0
# Certify the theta ordering by rational interval arithmetic, not midpoints.
ints = s.Poly(g, x).intervals(eps=s.Rational(1,10**30))
def interval(poly, ab):
    lo = hi = s.Rational(0)
    for c in s.Poly(poly, x).all_coeffs():
        vals = [lo*ab[0], lo*ab[1], hi*ab[0], hi*ab[1]]
        lo, hi = min(vals)+c, max(vals)+c
    return lo, hi
bounds = [interval(th, ab) for ab, mult in ints]
order = [4, 0, 2, 3, 1]
assert all(bounds[i][1] < bounds[j][0] for i,j in zip(order,order[1:]))
real1 = sum(1 << j for j in order[1:])
real2 = sum(1 << j for j in order[3:])
assert {0,real1,real2,real1^real2} == {0,5,10,15}
print('Exact real Kummer image: span{10100,01010}, increasing-a order.')
print('Known point has mask 4; exact square identity passed.')
print('All certificate checks passed. Conditional Selmer dimension 3; 1 <= rank <= 3 under GRH.')
