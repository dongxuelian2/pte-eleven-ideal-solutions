#!/usr/bin/env python3
"""Exact local certificate for the five-common-center-pair symmetric PTE ansatz.

F(x)=(x+10)*(H(x-1)+lambda), H(t)=product(t^2-k^2), k=1,3,5,7,9.
F(0)=10*lambda. For each prime p>5, test all nonzero fibers of H over F_p.
S consists of tested primes with no completely split nonzero fiber. Every
prime in S divides the constant of any integer solution in this ansatz.
Two exact methods independently test splitting, allowing repeated roots.
No third-party dependencies. Run with --limit 1000 for the smaller bound.
"""
import argparse
from collections import defaultdict
from math import comb, isqrt, prod
from time import perf_counter

H = [1]  # ascending coefficients
for k in (1,3,5,7,9):
    out = [0]*(len(H)+2)
    for i,v in enumerate(H):
        out[i] -= k*k*v
        out[i+2] += v
    H = out
TAYLOR = [[comb(k,j)*H[k] for k in range(j,11)] for j in range(11)]


def evaluate(coefficients,x,p):
    value = 0
    for coefficient in coefficients[::-1]:
        value = (value*x+coefficient) % p
    return value


def allowed_lambdas(p):
    """Method 1: total exact Taylor multiplicities in each nonzero fiber."""
    multiplicity = defaultdict(int)
    for x in range(p):
        lam = -evaluate(H,x,p) % p
        if lam:
            m = next(j for j in range(1,11) if evaluate(TAYLOR[j],x,p))
            multiplicity[lam] += m
    return sorted(lam for lam,m in multiplicity.items() if m == 10)


def independent_check(p):
    """Method 2: divide H+lambda repeatedly by every root in its fiber.

    Each possible root is enumerated. A constant quotient proves complete
    splitting; positive residual degree proves failure. This method uses no
    Taylor coefficients, derivatives, or root-multiplicity counts.
    """
    fibers = defaultdict(list)
    for x in range(p):
        lam = -evaluate(H,x,p) % p
        if lam:
            fibers[lam].append(x)
    result = []
    for lam, roots in fibers.items():
        poly = [a % p for a in H[::-1]]  # descending coefficients
        poly[-1] = (poly[-1]+lam) % p
        for x in roots:
            while len(poly) > 1:
                quotient = [poly[0]]
                for coefficient in poly[1:-1]:
                    quotient.append((coefficient+x*quotient[-1]) % p)
                if (poly[-1]+x*quotient[-1]) % p:
                    break
                poly = quotient
        if len(poly) == 1:
            result.append(lam)
    return sorted(result)


def ceiling_nth_root(value,n):
    lo, hi = 0, 10**((len(str(value))+n-1)//n)
    while lo+1 < hi:
        middle = (lo+hi)//2
        if middle**n >= value:
            hi = middle
        else:
            lo = middle
    assert (hi-1)**n < value <= hi**n
    return hi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--limit',type=int,default=5000)
    args = parser.parse_args()
    if args.limit < 7:
        parser.error('--limit must be at least 7')
    started = perf_counter()
    primes = [p for p in range(7,args.limit+1)
              if all(p%k for k in range(2,isqrt(p)+1))]
    obstruction_primes = []
    for p in primes:
        values = allowed_lambdas(p)
        other_values = independent_check(p)
        assert values == other_values, (p,values,other_values)
        if not values:
            obstruction_primes.append(p)
        if p == 349:
            assert values == [283]
            print('Positive control: p=349, admissible nonzero lambda=[283].')
    divisor = prod(obstruction_primes)
    height_bound = ceiling_nth_root(divisor,11)
    digits = len(str(divisor))
    if args.limit >= 1000:
        lower_primes = [p for p in obstruction_primes if p <= 1000]
        assert len(lower_primes) == 149
        assert len(str(prod(lower_primes))) == 371
        assert ceiling_nth_root(prod(lower_primes),11) == 4501915390982750858550567111089423
        print('Cutoff 1000 verified: 149 obstruction primes; product has 371 digits.')
    if args.limit == 5000:
        assert len(obstruction_primes) == 393
        assert digits == 1195
        assert height_bound == 4264914784302653186685253811259318130882488669797076480249850785413789679541836572698823410736766760419030769
    print(f'Cutoff {args.limit}: {len(primes)} primes tested by both exact methods.')
    print(f'Obstruction primes: {len(obstruction_primes)}; product D has {digits} digits.')
    print(f'Conditional integer height bound ceil(D^(1/11)): {height_bound}')
    print(f'All checks passed in {perf_counter()-started:.3f} seconds.')

if __name__ == '__main__':
    main()
