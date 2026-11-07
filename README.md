# Ideal Prouhet–Tarry–Escott Solutions of Size 11

**Research target:** The existence or obstruction of a nontrivial ideal Prouhet–Tarry–Escott (PTE) solution with 11 terms on each side.

Find integer multisets \(A=\{a_1,\ldots,a_{11}\}\) and \(B=\{b_1,\ldots,b_{11}\}\) such that
\[
\sum_{i=1}^{11}a_i^k=\sum_{i=1}^{11}b_i^k \qquad (k=1,\ldots,10),
\]
with \(A\neq B\). Equivalently, the associated monic root polynomials differ by a **nonzero constant**:
\[
\prod_{i=1}^{11}(x-a_i)-\prod_{i=1}^{11}(x-b_i)=C\neq0.
\]

## Questions

1. **Existence:** Construct an explicit nontrivial ideal size-11 pair, or prove a rigorous obstruction under useful hypotheses.
2. **Structured families:** Develop or rule out common-center / paired constructions without confusing results for a restricted family with the full problem.
3. **Arithmetic constraints:** Derive effective divisibility, congruence, and \(p\)-adic restrictions on putative solutions.
4. **Polynomial approach:** Use coefficient cancellation and the constant-difference identity to seek decisive general structure.

## Work

Research arguments go in [`PROOF.md`](PROOF.md). Introduce programs only for mathematical exploration or verifiable finite statements, and preserve the distinction between finite searches and general proofs.

**Status:** Open research target; no size-11 ideal solution or general impossibility theorem is claimed.

## Established conditional results

The complete proofs from
[ddy314/pte-order11-research, draft PR #1](https://github.com/ddy314/pte-order11-research/pull/1)
are consolidated in [PROOF.md, Section 3](PROOF.md#3-common-center-pair-obstructions).
They concern **symmetric solutions with distinct entries**, with the specified
disjoint pairs on one side sharing a common center:

- **Five pairs and a singleton:** an exact one-parameter normal form and
  height lower bound of approximately 4.2649 × 10^108, measured after centering
  the reflection symmetry. The exact integer bound is in Theorem 3.
- **Four pairs:** an exact two-parameter normal form and the forced divisors
  23, 29, 37, 41, 43, 47, 53, 59, 67, 73 of the constant C.
- **Three pairs:** 29 divides C, using the published CMSV local classification.

These are necessary conditions for restricted subclasses, not an integer
construction or a nonexistence proof for the general multiset problem.
This transfer makes no new novelty or record claim.

Run the exact certificates with Python 3.8 or newer, without dependencies:

```sh
python3 scripts/check_four_pairs_mod.py
python3 scripts/check_five_pairs_mod.py
```

The five-pair certificate checks all 666 primes from 7 through 5000 by two
exact splitting methods, yielding 393 obstruction primes. Both certificates
allow repeated roots modulo a prime, even though the integer entries in these
theorems are distinct.
