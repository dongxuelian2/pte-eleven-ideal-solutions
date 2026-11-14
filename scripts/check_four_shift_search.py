"""Rebuild the four reductions and independently check the bounded-search hits.

Usage: python scripts/check_four_shift_search.py [build/four_types_1000000.txt]
With no argument, checks identities, the six controls, and a direct height-100
enumeration against the C++ executable, if it has been compiled in build/.
Requires the existing SymPy package; no network access or software installation.
"""
from collections import Counter
from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path
import re
import subprocess
import sys

import sympy as s

c = s.symbols("c")
FAMILIES = {
    4: ((20*c+25)/6, (240*c*c+40*c+167)/48,
        (44*c+91)/12, (256*c*c+97*c+466)/36),
    5: ((20*c-85)/16, (80*c*c-840*c+589)/256,
        (12*c-83)/16, (4*c-49)*(4*c-9)/256),
    15: ((20*c+15)/8, (160*c*c-20*c+49)/96,
         (69-4*c)/8, (2*c-23)*(4*c-31)/48),
    16: ((195-20*c)/12, (160*c*c-1335*c+2304)/36,
         4*(9-c)/3, 2*(c-6)*(4*c-9)/3),
}
DATA = {4:(1,2,1,3), 5:(1,2,2,1), 15:(3,1,1,3), 16:(3,1,2,1)}
EXPECTED = {(4,1,2), (5,7,2), (5,9,2), (15,1,2), (15,7,2), (16,3,2)}


def offsets(length):
    return [s.Rational(2*j-length+1, 2) for j in range(length)]


def moment(length, k, square):
    return sum(s.binomial(2*k,2*j)*square**(k-j)
               *sum(x**(2*j) for x in offsets(length)) for j in range(k+1))


def paired_moment(length, k, total, product):
    powers = [s.Integer(2), total]
    for j in range(2, k+1):
        powers.append(s.expand(total*powers[-1]-product*powers[-2]))
    z = s.symbols("z")
    return sum(coeff*powers[mon[0]] for mon,coeff in s.Poly(moment(length,k,z),z).terms())


for typ, (S,T,U,V) in FAMILIES.items():
    m, la, lb, lc = DATA[typ]
    for k in range(1,5):
        residual = (sum(x**(2*k) for x in offsets(2*m))/2
                    +paired_moment(la,k,S,T)-paired_moment(lb,k,U,V)-moment(lc,k,c))
        assert s.cancel(residual) == 0
print("PASS: four exact moment identities rebuilt for each type.")

# Check the actual integer tables consumed by the C++ executable against the
# independently reconstructed rational center quadratics, including scale.
cpp = Path(__file__).with_name("search_four_shift_types.cpp").read_text()
rows = re.findall(r"\{(4|5|15|16),\s*(\d+),\s*\{([^}]+)\},\s*\{([^}]+)\},\s*\{([^}]+)\},\s*\{([^}]+)\}\}", cpp)
assert len(rows) == 4
for typ,D,ss,uu,tt,vv in rows:
    typ,D = int(typ),int(D)
    poly = lambda raw: s.Poly.from_list(list(map(int,raw.split(","))),c).as_expr()
    S,T,U,V = FAMILIES[typ]
    assert s.expand(poly(ss)-D*S) == s.expand(poly(uu)-D*U) == 0
    assert s.expand(poly(tt)-D*D*T) == s.expand(poly(vv)-D*D*V) == 0
print("PASS: every C++ integer coefficient and scaling matches the moment reduction.")


def sqrtq(x):
    if x < 0:
        return None
    n,d = isqrt(x.numerator), isqrt(x.denominator)
    return Q(n,d) if n*n == x.numerator and d*d == x.denominator else None


# Convert from SymPy once; the direct small search uses Fraction arithmetic.
COEFFICIENTS = {typ:[[Q(a) for a in s.Poly(f,c).all_coeffs()] for f in fs]
                for typ,fs in FAMILIES.items()}


def eval_poly(coeffs, value):
    out = Q(0)
    for x in coeffs:
        out = out*value+x
    return out


def centers(typ, value):
    S,T,U,V = [eval_poly(f,value) for f in COEFFICIENTS[typ]]
    out = []
    for total,product in ((S,T),(U,V)):
        r = sqrtq(total*total-4*product)
        if r is None:
            return None
        pair = [sqrtq((total+r)/2), sqrtq((total-r)/2)]
        if None in pair:
            return None
        out.append(pair)
    return out


def classify(hit):
    typ,n,d = hit
    value = Q(n,d)**2
    ca,cb = centers(typ,value)
    m,la,lb,lc = DATA[typ]
    def chains(cs,length):
        return [sign*u+Q(v) for u in cs for sign in (-1,1) for v in offsets(length)]
    A = list(map(Q,offsets(2*m)))+chains(ca,la)
    B = chains(cb,lb)+chains([Q(n,d)],lc)
    assert len(A) == len(B) == 10
    for k in range(10):
        assert sum(x**k for x in A) == sum(x**k for x in B)
    problems = []
    if len(set(A)) < 10: problems.append("repeated A entries")
    if len(set(B)) < 10: problems.append("repeated B entries")
    if set(A) & set(B): problems.append("A/B overlap")
    if problems:
        return ", ".join(problems)
    left = Counter(A)+Counter(b+1 for b in B)
    right = Counter(B)+Counter(a+1 for a in A)
    shared = left & right
    left -= shared
    right -= shared
    assert sum(left.values()) == sum(right.values())
    if sum(left.values()) != 11 or max(left.values()) != 1 or max(right.values()) != 1:
        return "not a distinct-entry size-eleven output"
    for k in range(11):
        assert sum(x**k for x in left) == sum(x**k for x in right)
    raise AssertionError(f"VALID WITNESS FOUND: {hit}, {sorted(left)}, {sorted(right)}")


def parse(text):
    hits = [tuple(map(int,line.split()[1:])) for line in text.splitlines() if line.startswith("HIT ")]
    assert len(hits) == len(set(hits))
    match = re.search(r"^COMPLETE candidates=\d+ modular_survivors=\d+ exact_hits=(\d+) seconds=", text, re.M)
    assert match and int(match[1]) == len(hits), "incomplete or inconsistent search"
    return set(hits)


if len(sys.argv) > 1:
    text = Path(sys.argv[1]).read_text()
    bound = re.search(r"^BOUND (\d+) DENOMINATORS (\d+) (\d+)$", text, re.M)
    assert bound and int(bound[2]) == 1 and int(bound[3]) == int(bound[1])
    hits = parse(text)
    assert all(0 <= n <= int(bound[1]) and 1 <= d <= int(bound[1]) and gcd(n,d)==1 for _,n,d in hits)
else:
    hits = {(typ,n,d) for typ in FAMILIES for d in range(1,101)
            for n in range(101) if gcd(n,d)==1 and centers(typ,Q(n,d)**2) is not None}
    binary = Path(__file__).resolve().parents[1]/"build"/"search_four_shift_types"
    if binary.exists():
        text = subprocess.check_output([str(binary),"100","2"],text=True)
        assert parse(text) == hits
        print("PASS: independent Fraction enumeration and C++ sieve agree through height 100.")

assert hits == EXPECTED
for hit in sorted(hits):
    print(hit, classify(hit))
print("PASS: all six center-square configurations fail the distinct/disjoint seed conditions.")
print("No unbounded rational-point or unrestricted PTE exclusion is asserted.")
