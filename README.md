# Ideal Prouhet–Tarry–Escott Solutions of Size 11

**Research target:** Construct two disjoint sets of eleven distinct integers
forming an ideal Prouhet–Tarry–Escott (PTE) solution, or prove general
nonexistence.

Find integer sets \(A=\{a_1,\ldots,a_{11}\}\) and \(B=\{b_1,\ldots,b_{11}\}\),
with all eleven entries distinct on each side and \(A\cap B=\varnothing\), such that
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

**Scope:** Section 3 gives conditional results for distinct-entry symmetric
families. Section 4 studies a multiset relaxation. Its repeated-root and
support results do not solve or exclude the original distinct-entry task.

## Distinct-entry shift-lift route

[PROOF.md, Section 5](PROOF.md#5-a-distinct-entry-construction-route-shift-chains)
gives an exact criterion for lifting a symmetric size-ten solution to
eleven distinct entries. Nine same-side cancellations and no cross-side
unit differences yield eighteen possible normalized chain types.

For four types, exact elimination leaves rational conics and a cubic
whose three roots must be rational squares. Their discriminants give
four explicit smooth genus-two curves. Faltings's theorem implies
finitely many normalized valid configurations in those four types,
but does not enumerate them or prove their absence. The known control
points are trivial A=B and are rejected.

A fifth type reduces to an explicit smooth plane quartic of genus three,
with a displayed elliptic quotient. It also has finitely many normalized
valid configurations, without a complete rational-point classification.
A sixth type also maps to a smooth genus-three quartic. Its explicit
nonexceptional rational point recovers a nonsquare center, showing why
a rational quartic point alone is insufficient. This gives normalized
finiteness, not a solution or an effective point classification.
The sixth-mode Jacobian has certified rational rank seven: its three
elliptic factors have ranks 2, 2, 3 by exact 2-isogeny descent. Real and
modulo-256 exclusions give upper bounds, and explicit points attain every
remaining descent class. Thus ordinary Chabauty's rank < genus hypothesis
fails here. This does not determine the rational points or settle PTE
existence. The other twelve types and other construction routes are not
settled by these reductions.

The polynomial certificates require SymPy; the exact-rank descent
certificate is dependency-free:

```sh
python3 scripts/check_shift_chains.py
python3 scripts/check_fifth_mode.py
python3 scripts/check_sixth_mode.py
python3 scripts/check_sixth_rank_maps.py
python3 scripts/check_sixth_exact_ranks.py
```

## Ancillary multiset results

[PROOF.md, Section 4](PROOF.md#4-general-multiset-restrictions) proves, without
symmetry or distinct-entry assumptions, that every nontrivial rational
size-11 ideal pair has:

- Disjoint supports, with multiplicity at most two for every entry.
- Ten distinct real zeros of the common derivative of its root polynomials.
- At least fourteen distinct values in the two supports combined.
  Six/six forces an irrational Chebyshev factor; six/seven forces a
  rational elliptic point of order 11, excluded by Mazur's theorem.
- At least eight distinct values per side if the pair is symmetric.
  The seven/seven pattern descends explicitly to a quartic polynomial
  Pell equation and yields the same elliptic torsion contradiction.
- An exact higher-genus Pell/Jacobian reduction for the remaining patterns.
  The smallest remaining cases reach genus two, where the elliptic
  obstruction alone no longer applies.
- At this genus-two boundary, a Jacobian simple over Q with rational
  torsion containing [2,2,2,22] (order 176) is necessary. Rational
  elliptic quotients, including centrally paired sextic branch sets,
  are excluded.
- An exact degree-eleven sextic Pell example shows why Pell solvability
  alone is insufficient: its sextic does not split over Q, so it is
  not a PTE construction.
- For symmetric eight-value supports, the three double entries cannot
  sum to zero; a second explicit exceptional parameter locus is also
  excluded. Theorem 12 gives the unique quintic factor and an exact
  two-parameter rational splitting problem after scaling.
- A fully real eight-value example with double entries 1, 5, 9 shows
  that real-root arguments cannot finish the exclusion. Its quintic
  fails complete rational splitting already modulo 5.
- The fourteen-value general boundary is exactly equivalent to a
  degree-eleven Pell identity with a completely rationally split
  squarefree sextic D and a completely rationally split degree-eight T.
  This is an equivalence with simultaneous splitting, not with torsion
  alone, and supplies no integer solution.

These are necessary conditions or exact reformulations within the multiset
relaxation. In particular, the repeated-root boundary cases do not meet
the original requirement of eleven distinct entries on each side.
They do not prove existence, global nonexistence, or attainment of the
support bound. The Chebyshev and Pell/torsion mechanisms are classical; Mazur's
established torsion theorem is cited explicitly. No novelty claim is made.

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
construction or a nonexistence proof for the original distinct-entry problem.
This transfer makes no new novelty or record claim.

Run the exact certificates with Python 3.8 or newer, without dependencies:

```sh
python3 scripts/check_four_pairs_mod.py
python3 scripts/check_five_pairs_mod.py
python3 scripts/check_sextic_pell.py
```

The five-pair certificate checks all 666 primes from 7 through 5000 by two
exact splitting methods, yielding 393 obstruction primes. Both certificates
allow repeated roots modulo a prime, even though the integer entries in these
theorems are distinct.
