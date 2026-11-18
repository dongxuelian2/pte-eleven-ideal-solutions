#!/usr/bin/env python3
"""Independently verify a Type18 rigid-prime CSV and its denominator bound.

Usage: python scripts/verify_type18_rigid_primes.py certificate.csv
       python scripts/verify_type18_rigid_primes.py certificate.csv \
           --denominator-bound 1000000
       python scripts/verify_type18_rigid_primes.py --self-test

The CSV has one column, headed q (an unheaded one-column list also works).
A prime q>5 is rigid when c=0 is the only c for which P(X)-c splits
completely with all roots square residues, where P=prod(X-a^2),
a=1,3,5,7,9. Repeated roots count with their polynomial multiplicities.
Multiplicities here use explicit formal derivatives, not synthetic division.
Only the Python standard library is required.
"""

import argparse
import csv
from itertools import combinations_with_replacement
from math import isqrt, prod
from pathlib import Path
import sys


# Coefficients in increasing order, independently expanded from the five roots.
SEED = (1, 3, 5, 7, 9)
P = (-893025, 1057221, -172810, 8778, -165, 1)
HEIGHT_CONSTANT = 33**5


def derivative(coefficients):
    return tuple(i * a for i, a in enumerate(coefficients) if i)


DERIVATIVES = []
_derivative = P
for _ in range(5):
    _derivative = derivative(_derivative)
    DERIVATIVES.append(_derivative)


def evaluate(coefficients, x, modulus):
    value = 0
    for a in reversed(coefficients):
        value = (value * x + a) % modulus
    return value


def is_prime(q):
    """Exact trial division; the O(q) fiber check dominates this O(sqrt(q)) test."""
    if q < 2:
        return False
    if q % 2 == 0:
        return q == 2
    return all(q % d for d in range(3, isqrt(q) + 1, 2))


def split_square_fibers(q):
    """Return {c: [(root, multiplicity), ...]} for all full square fibers."""
    residues = {x * x % q for x in range((q + 1) // 2)}
    fibers = {}
    for root in residues:
        c = evaluate(P, root, q)
        # q>5, so the fifth derivative 5! is nonzero and some order occurs.
        multiplicity = next(
            order for order, d in enumerate(DERIVATIVES, 1)
            if evaluate(d, root, q)
        )
        fibers.setdefault(c, []).append((root, multiplicity))
    return {
        c: roots for c, roots in fibers.items()
        if sum(m for _, m in roots) == 5
    }


def integer_nth_root(n, exponent):
    """Floor root, using only integer comparisons and a binary interval."""
    if n < 0 or exponent < 1:
        raise ValueError("root requires n>=0 and exponent>=1")
    if n < 2:
        return n
    low = 0
    high = 1 << ((n.bit_length() + exponent - 1) // exponent)
    while high - low > 1:
        middle = (low + high) // 2
        if middle**exponent <= n:
            low = middle
        else:
            high = middle
    return low


def strict_denominator_bound(modulus, exponent):
    # C*D^exponent < modulus iff D^exponent <= floor((modulus-1)/C).
    return integer_nth_root((modulus - 1) // HEIGHT_CONSTANT, exponent)


def read_primes(path):
    primes = []
    seen = set()
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = csv.reader(stream)
        first = True
        for line, row in enumerate(rows, 1):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) != 1:
                raise ValueError(f"line {line}: expected exactly one q column")
            value = row[0].strip()
            if first and value == "q":
                first = False
                continue
            first = False
            try:
                q = int(value)
            except ValueError as exc:
                raise ValueError(f"line {line}: invalid integer {value!r}") from exc
            if q in seen:
                raise ValueError(f"line {line}: duplicate value {q}")
            if q <= 5:
                raise ValueError(f"line {line}: q must exceed 5, got {q}")
            if not is_prime(q):
                raise ValueError(f"line {line}: {q} is not prime")
            seen.add(q)
            primes.append(q)
    if not primes:
        raise ValueError("empty certificate")
    return sorted(primes)


def verify(primes, verbose=False):
    for q in primes:
        fibers = split_square_fibers(q)
        if 0 not in fibers:
            raise ValueError(f"internal polynomial error at q={q}: missing seed fiber")
        witnesses = [c for c in fibers if c != 0]
        if witnesses:
            c = min(witnesses)
            raise ValueError(
                f"q={q} is not rigid: c={c}, square roots/multiplicities="
                f"{sorted(fibers[c])}"
            )
        if verbose:
            print(f"PASS q={q} class={q % 10}")
    m0 = prod(q for q in primes if q % 10 != 1)
    m1 = prod(q for q in primes if q % 10 == 1)
    modulus = m0 * m1
    d10 = strict_denominator_bound(m0, 10)
    d11 = strict_denominator_bound(modulus, 11)
    for m, e, d in ((m0, 10, d10), (modulus, 11, d11)):
        assert HEIGHT_CONSTANT * d**e < m
        assert HEIGHT_CONSTANT * (d + 1)**e >= m
    return m0, m1, modulus, d10, d11


def self_test():
    coefficients = [1]
    for a in SEED:
        new = [0] * (len(coefficients) + 1)
        for i, c in enumerate(coefficients):
            new[i] -= a * a * c
            new[i + 1] += c
        coefficients = new
    assert tuple(coefficients) == P
    assert P[0] == -(prod(SEED) ** 2) == -(945**2)
    # Independent exhaustive five-root multisets: Newton's identities say
    # the first four root power sums fix precisely the first four coefficients.
    for q in (7, 11, 13):
        residues = sorted({x * x % q for x in range(q)})
        target = tuple(sum(pow(a * a, j, q) for a in SEED) % q
                       for j in range(1, 5))
        brute = set()
        for roots in combinations_with_replacement(residues, 5):
            moments = tuple(sum(pow(x, j, q) for x in roots) % q
                            for j in range(1, 5))
            if moments == target:
                brute.add(evaluate(P, roots[0], q))
        assert set(split_square_fibers(q)) == brute
    for exponent in (1, 2, 10, 11):
        for root in range(1, 20):
            assert integer_nth_root(root**exponent, exponent) == root
            if exponent > 1:
                assert integer_nth_root(root**exponent - 1, exponent) == root - 1
    for exponent in (10, 11):
        for d in (1, 2, 10):
            assert strict_denominator_bound(HEIGHT_CONSTANT * d**exponent,
                                            exponent) == d - 1
            assert strict_denominator_bound(HEIGHT_CONSTANT * d**exponent + 1,
                                            exponent) == d
    print("PASS: seed constant, derivative fibers versus root multisets, exact strict roots")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, nargs="?")
    parser.add_argument("--denominator-bound", type=int)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    if args.csv is None:
        if args.self_test:
            return 0
        parser.error("a CSV certificate is required unless --self-test is used")
    if args.denominator_bound is not None and args.denominator_bound < 1:
        parser.error("--denominator-bound must be positive")
    primes = read_primes(args.csv)
    m0, m1, modulus, d10, d11 = verify(primes, args.verbose)
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(0)
    print(f"PASS: {len(primes)} listed primes are >5, prime, distinct, and rigid")
    print(f"count0={sum(q % 10 != 1 for q in primes)}")
    print(f"count1={sum(q % 10 == 1 for q in primes)}")
    print(f"M0={m0}")
    print(f"M1={m1}")
    print(f"M={modulus}")
    print(f"Dmax10={d10} (M0 > 33^5*D^10)")
    print(f"Dmax11={d11} (M > 33^5*D^11)")
    print(f"Dmax={max(d10, d11)}")
    if args.denominator_bound is not None:
        if args.denominator_bound > max(d10, d11):
            raise ValueError("requested denominator bound is not certified")
        print(f"PASS: nontrivial Type18 tuples with D<={args.denominator_bound} excluded")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
