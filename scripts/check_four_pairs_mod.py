#!/usr/bin/env python3
"""Exact finite-field certificate for symmetric size-11 PTE, four paired roots.
No third-party dependencies. Repeated roots are explicitly allowed.
Normal form: F(x)=Q((x-1)^2)*P(x-1), P(t)=t^3+11*t^2+b*t+d.
For integer common center c, or half-integer c, normalize by c modulo odd p.
If c=0 modulo p, any common-center pair already forces the PTE constant=0.
Zero survivors therefore imply p divides the integer PTE constant.
This verifies a conditional obstruction, not nonexistence of PTE solutions.
"""
from collections import Counter
from math import comb, isqrt


def splits(poly, p):
    """Monic polynomial, descending coefficients; divide every linear factor."""
    poly = [a % p for a in poly]
    for a in range(p):
        while len(poly) > 1:
            quotient = [poly[0]]
            for coefficient in poly[1:-1]:
                quotient.append((coefficient + a * quotient[-1]) % p)
            if (poly[-1] + a * quotient[-1]) % p:
                break
            poly = quotient
    return len(poly) == 1


def q_coefficients(b, d, p):
    q3 = (9*b-d-330)*pow(2, -1, p) % p
    q2 = ((7*b-d-224)*q3+84*b-28*d-1848)*pow(4, -1, p) % p
    q1 = ((5*b-d-130)*q2+(35*b-15*d-644)*q3+126*b-70*d-1980)*pow(6, -1, p) % p
    q0 = ((3*b-d-56)*q1+(10*b-6*d-144)*q2+(21*b-15*d-272)*q3+36*b-28*d-440)*pow(8, -1, p) % p
    return [1, 0, q3, 0, q2, 0, q1, 0, q0]


def nonzero_constant(b, d, p):
    # 384 F(0) = -product[k*b-d-(11*k*k-k**3)], k=1,3,5,7,9.
    return all((k*b-d-(11*k*k-k**3)) % p for k in (1, 3, 5, 7, 9))


def survivors(p):
    return [(b, d) for b in range(p) for d in range(p)
            if nonzero_constant(b, d, p)
            and splits([1, 11, b, d], p)
            and splits(q_coefficients(b, d, p), p)]


def independent_check(p):
    """Enumerate cubic root multisets; count octic roots via Taylor coefficients.
    Independent of both the b,d enumeration and synthetic division test.
    """
    cubics = set()
    for r in range(p):
        for s in range(r, p):
            t = (-11-r-s) % p
            if t >= s:
                cubics.add(((r*s+r*t+s*t) % p, (-r*s*t) % p))
    result = []
    for b, d in sorted(cubics):
        if not nonzero_constant(b, d, p):
            continue
        coefficients = q_coefficients(b, d, p)[::-1]
        multiplicities = []
        for x in range(p):
            multiplicities.append(next(j for j in range(9)
                if sum(comb(k,j)*coefficients[k]*pow(x,k-j,p)
                       for k in range(j,9)) % p))
        if sum(multiplicities) == 8:
            result.append((b,d))
    return result


def pair_maximum(pattern, p):
    counts = Counter(pattern)
    def pairs(c):
        return (sum(min(counts[x], counts[(2*c-x) % p])
                    for x in range(p) if x != c)//2 + counts[c]//2)
    return max(pairs(c) for c in range(1,p))


def main():
    primes = [p for p in range(5,102)
              if all(p % k for k in range(2,isqrt(p)+1))]
    zero = []
    for p in primes:
        found = survivors(p)
        if not found:
            zero.append(p)
        print(f'p={p:3d}: {len(found):2d} surviving (b,d) pairs; first={found[:1]}')
        if p in (23,29,37,61):
            assert found == independent_check(p), ('Independent check failed',p)
    assert zero == [5,7,13,17,19,23,29,31,37,41,43,47,53,59,67,73]
    assert survivors(61) == [(58,16)]
    # These local pattern lists are from CMSV, arXiv:2304.11254v1, section 3.
    # Their completeness is a literature input, not established by this script.
    patterns = [(23,[1,2,3,4,6,8,9,12,13,16,18],3),
                (23,[1,2,2,3,3,5,5,7,7,17,17],3),
                (29,[1,1,1,15,19,21,21,22,23,25,25],2)]
    for p,pattern,expected in patterns:
        assert set(pattern).isdisjoint({-x % p for x in pattern})
        assert all(sum(pow(x,j,p) for x in pattern) % p == 0
                   for j in (1,3,5,7,9))
        actual = pair_maximum(pattern,p)
        assert actual == expected
        print(f'CMSV pattern mod {p}: maximum common-center pairs = {actual}')
    print('Zero-survivor primes:',zero)
    print('All independent checks passed.')

if __name__ == '__main__':
    main()
