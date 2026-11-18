"""Prove the norm bounds in §§5.55 and 5.60 by prime counts and rational intervals."""
from argparse import ArgumentParser
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
from functools import cache
from pathlib import Path
import subprocess

SUPPORT = 2**31
BOUND = 48_000_001
DISC = 172707117717588564754287819295488000000
SCALE = 2**64


@cache
def log_unit_interval(value):
    """Outward fixed-point bounds for log(value), for 1 <= value <= 2."""
    assert 1 <= value <= 2
    z = (value - 1) / (value + 1)
    terms = 24
    s = sum(2 * z**(2*j+1) / (2*j+1) for j in range(terms))
    tail = 2 * z**(2*terms+1) / ((2*terms+1) * (1-z*z))
    return int(s*SCALE), int((s+tail)*SCALE) + 1


def log_interval(value):
    value = Q(value)
    assert value >= 1
    exp = value.numerator.bit_length() - value.denominator.bit_length()
    if value < 2**exp:
        exp -= 1
    lo, hi = log_unit_interval(value / 2**exp)
    log2lo, log2hi = log_unit_interval(Q(2))
    return exp*log2lo + lo, exp*log2hi + hi


def read_histogram(lines):
    hist = Counter()
    for line in lines:
        e, b, m, count = map(int, line.split())
        assert 0 <= e <= 31 and 0 <= b < 256 and 1 <= m <= 31 and count > 0
        hist[e, b, m] += count
    return hist


def prove(hist, support=SUPPORT, bound=BOUND, disc=DISC, real_places=2):
    length_lo, length_hi = log_interval(support)
    prime_sum = 0
    for (e, b, m), count in hist.items():
        lower = Q(2**e * (256+b), 256)
        upper = min(Q(support), Q(2**e * (257+b), 256))
        log_lower = log_interval(lower)[0]
        log_upper = log_interval(upper)[1]
        if upper <= bound:
            term = (4*count*log_lower * Q(max(0, length_lo-log_upper), length_lo)
                    / (m*(upper+1)))
        else:
            term = (-4*count*log_upper * Q(length_hi-log_lower, length_hi)
                    / (m*(lower+1)))
        prime_sum += term.numerator // term.denominator

    # H_n - log(n) - 1/(2n) < Euler's constant by the trapezoid
    # bound for the convex decreasing function 1/x.
    gamma_lo = sum(Q(1, j) for j in range(1, 257)) - Q(log_interval(256)[1], SCALE) - Q(1, 512)
    pi_lo, pi_hi = Q(314159, 100000), Q(314160, 100000)
    ell = Q(length_lo, SCALE)
    arch = 12*(gamma_lo + Q(log_interval(4*pi_lo)[0], SCALE) - pi_hi**2/(4*ell))
    arch += real_places*(1 - 2*Q(log_interval(2)[1], SCALE)/ell)
    margin = Q(prime_sum, SCALE) + arch - Q(log_interval(disc)[1], SCALE)
    if disc == DISC:
        assert Q(prime_sum, SCALE) > Q(50333, 1000) and margin > Q(8, 100)
    else:
        assert margin > Q(26, 100)
    print(f"Norm bound {bound}, explicit-formula margin > {int(margin*1000)}/1000")


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('counter', type=Path, help='compiled prime_counts executable')
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--type17', action='store_true', help='verify the totally real Type 17 field')
    args = parser.parse_args()
    support = 2**20 if args.type17 else SUPPORT
    block = 10_000_000

    def subtotal(lo):
        hi = min(support, lo+block-1)
        command = [args.counter.resolve(), str(lo), str(hi)]
        if args.type17:
            command.append('type17')
        result = subprocess.run(command,
                                check=True, capture_output=True, text=True)
        return read_histogram(result.stdout.splitlines())

    hist = Counter()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for result in pool.map(subtotal, range(2, support+1, block)):
            hist.update(result)
    if args.type17:
        prove(hist, support, support, 58160039859394517078535559761100800000, 12)
    else:
        prove(hist)


if __name__ == '__main__':
    main()
