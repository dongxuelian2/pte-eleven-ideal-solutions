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

## Splitting-pencil Selmer computation

[PENCIL_SELMER.md](PENCIL_SELMER.md) gives a complete staged 2-descent for
the genus-two discriminant curve of an explicit quintic–sextic splitting
pencil. Its 2-Selmer dimension is unconditionally three, so 1≤rank J(Q)≤3.
Exact principal relations, verification through the Minkowski bound, and
global ideal-parity/sign checks remove the earlier GRH assumption. This is not a rank-one theorem or a valid
size-eleven splitting fiber. Run `python scripts/pencil_selmer/verify.py`
for the exact certificates; the Magma inputs and successful mathematical
outputs are included.

## Sixth-mode exclusion

[The complete sixth-mode proof](SIXTH_MODE_QC_EXCLUSION.md) now closes the local QC calculation and candidate sieve for the sixth shift-chain type. Certified moving normalization and analytic tails isolate48 local QC roots; exact elliptic tests at61,101,107 exclude every candidate except the eight known colliding orbit points. Thus this type has no valid distinct-entry lift. The new steps were independently reviewed, with upstream arithmetic dependencies listed explicitly. Reproduce them with scripts/check_sixth_mode_qc_roots13.py and scripts/check_sixth_mode_candidate_sieve.py. This is not a full rational-point classification or a solution of general size-eleven PTE; the other17 chain types are not excluded.

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
existence.

Adding the third center-square condition gives an explicit genus-nine
cover with a rank-one elliptic projection. A birational map and exact
2-isogeny descent certify the projection. The simultaneous-square lifts,
remaining singleton roots and collision tests are not classified.

The elliptic projection's full group is now certified as
Z·(180,900) plus the order-two point (0,0). An explicit height inequality
forces any missing generator below scaled height 24, where an exhaustive
exact check rules it out. Thus all elliptic candidates are integer
multiples of the certified generator, with an optional torsion translate.
Which candidates satisfy all lifting conditions remains unresolved.

Four good-prime reductions give uniform index restrictions modulo 240.
Every surviving candidate with |n|<120 has a collision; any valid lift
therefore has scaled elliptic-coordinate height H(x/270)>10^3298.
This is not a bound on the original integer entries or a global exclusion.
A finite covering descent gives eight explicit genus-two quotient models,
with squareclass 1 or 5 fixed by the parity of n. All eight Jacobians
have certified positive rank and are Q-simple. Exact ranks for the six
curves with j other than 2 and all complete rational-point sets remain
unproved.

The genus-nine center cover itself has Jacobian rank at least twelve,
so ordinary Chabauty does not apply there either. For the j=2 genus-two
pair, explicit Richelot duals have torsion exactly (Z/2)^3 and known
infinite-order classes. Full two-descent now computes both dual 2-Selmer
groups with dimension five; explicit independent divisor classes prove
both ranks exactly two. The same exact ranks hold for the original j=2
Jacobians, so ordinary genus-two Chabauty fails its rank condition here.
For the duals, Sha[2]=0 and the displayed subgroup has finite odd index.
Saturation at 13 is now certified below; saturation at the other odd
primes and complete original-curve point sets are not certified. An isogeny of Jacobians does not identify the two curves'
rational-point sets.
The two selected genus-two Jacobians have rational Neron-Severi rank one,
so the standard quadratic-Chabauty inequality does not hold on them.
However, all sixth-mode center lifts lie on two explicit genus-three covers
with Jacobian rank three and rational Neron-Severi rank two. Their
Mordell-Weil-restricted quadratic-Chabauty finiteness criterion does hold.
The models, maps, rational base points, trace-zero correspondence and
ordinary good prime 13 are certified. This establishes an applicable route,
not a completed computation: p-adic heights, sufficient Mordell-Weil data,
rigorous root isolation and the final point-completeness sieve remain.
The covers need not contain every point of the original genus-two curves,
and their points must still pass the remaining PTE conditions.

Exact de Rham inputs are now explicit: a full rational basis, signed cup
matrix, holomorphic-first symplectic normalization and the genuine
trace-zero correspondence. A symbolic Frobenius-projector identity and all
twelve reduction discs at 13 are checked.

Exact arithmetic modulo 169 and tiny-integral error bounds now certify
13-adic logarithm rank three for both genus-three covers. On each dual
genus-two factor the normalized logarithm determinant is 6 modulo 13;
the elliptic factor supplies the third independent logarithm.
The known subgroups are 13-saturated, proved both from the logarithm
lattice and independently from reductions at 181 and 199. This does
not provide explicit transported divisor coordinates or a full
target-basis logarithm matrix. No full genus-three height, isolated QC root
or completed point-completeness sieve has been computed.

Section 5.17 now computes the elliptic-factor Coleman–Gross height for the
specified quotient splitting: h(P,P)=59106874840+O(13^10). Exact formal
series and rigorous tail bounds agree through 6P, 12P and 24P. The nonzero
away-from-13 contribution and all rational elliptic-point component
corrections are explicit. The genus-two height block, full genus-three
local contributions, QC roots and completeness sieve remain uncomputed.

Section 5.18 replaces the missing genus-two self-pairings with three mixed
pairings of explicit disjoint divisors. It gives residue-correct third-kind
differentials and computes every finite local contribution outside the bad
prime set and 13, including nonzero terms at 7,17,29,31,313,593. The vertical
corrections of these representatives vanish at 11,19,53. Section 5.19 now
finishes their horizontal intersections using explicit blowups: 11 and 53
contribute zero, while each curve has one nonzero mixed term +2log13(19).
Together these certify every finite contribution outside {2,3,5,13}.
Section 5.20 now also certifies the terms at 3 and 5 using explicit
quadratic-extension charts and local-height descent: one mixed term is
-log13(3), another is -2log13(5), and the others vanish. Thus every finite
contribution outside {2,13} is known for these representatives. Section 5.21
now resolves the dyadic terms too: all six vanish, using explicit smooth
charts and same/opposite-sheet cancellation over a ramified quadratic field.
All away-from-13 mixed-height inputs are therefore complete, with certified
13-adic residues.

Section 5.22 now makes the Richelot splitting transport explicit, including
the full de Rham matrix, exact primitives, cup scaling and compatible
height weighting. Section 5.23 computes both genus-two Frobenius matrices
modulo 13^3 with a proved infinite-tail bound. Independent exact-rational
and modular algorithms agree, as do cup and quadratic-twist checks.
The compatible splitting is therefore known.

Section 5.24 now evaluates all six Coleman basis-integral vectors modulo
13^3, with retained primitives, certified endpoint arithmetic and a separate
principal-divisor reciprocity check. Their sign-correct third-kind
normalizations are explicit. Section 5.25 computes three actual normalized
mixed heights via good even-degree models, independently reproduced with
proved error bounds. Combining the finite contributions gives h(P,P)=923
and h(P,Q)=1586 modulo 13^3 for the k=-15 Jacobian, and h(P,P)-h(P,Q)=1417
for k=-3. Sections 5.26--5.27 complete both global genus-two height
matrices modulo13^3, including the remaining Q/R integrals and an exact
moved divisor T=4P-Q that avoids the shared-disc integration. Its new finite
corrections are recomputed in full. The matrices are
[[923,1586],[1586,689]] for k=-15 and [[2093,676],[676,1027]] for k=-3.
They yield certified global height quadratic forms in normalized logarithms.
Independent implementations and precision proofs agree. These are completed
global height blocks, not the local genus-three QC functions: basepoint
normalization, allowed finite local-value sets, root isolation and the final
rational-point sieve remain unfinished.

Sections 5.28--5.29 supply the corrected curve-level QC divisor, including
the basepoint-dependent linear term, and an explicit rational differential
with the required infinity cases. Exact regular-model and graph-symmetry
arguments show that BOTH curves have local target {0} at5,11,19,53;
all other good primes also contribute zero. Only the finite target sets
at2 and3 remain. These conclusions concern the correctly normalized
curve-level height, not the earlier generator-pair local contributions.
The moving13-adic de Rham projection, certified QC series, roots and
final sieve remain uncomputed. No size-eleven solution is claimed.

The other twelve types and other construction routes are not settled
by these reductions.

The polynomial and de Rham certificates require SymPy; the exact-rank descent
and generator-saturation certificates are dependency-free:

```sh
python3 scripts/check_shift_chains.py
python3 scripts/check_fifth_mode.py
python3 scripts/check_sixth_mode.py
python3 scripts/check_sixth_rank_maps.py
python3 scripts/check_sixth_exact_ranks.py
python3 scripts/check_sixth_center_lift.py
python3 scripts/check_rank_one_generator.py
python3 scripts/check_lift_congruences.py
python3 scripts/check_sixth_genus_two_covers.py
python3 scripts/check_center_cover_rank_obstruction.py
python3 scripts/check_genus_two_richelot_certificate.py
python3 scripts/check_richelot_two_selmer.py
python3 scripts/check_genus_three_qc_targets.py
python3 scripts/check_genus_three_qc_de_rham.py
python3 scripts/check_thirteen_logarithms.py
python3 scripts/check_thirteen_saturation.py
python3 scripts/check_elliptic_height13.py
python3 scripts/check_genus_two_mixed_heights.py
python3 scripts/check_genus_two_semistable_heights.py
python3 scripts/check_genus_two_odd_local_heights.py
python3 scripts/check_genus_two_dyadic_heights.py
python3 scripts/check_richelot_de_rham_transport.py
python3 scripts/check_genus_two_frobenius13.py
python3 scripts/check_coleman_basis13.py
python3 scripts/check_normalized_mixed_heights13.py
python3 scripts/check_remaining_qr_heights13.py
python3 scripts/check_complete_height_blocks13.py
python3 scripts/check_qc_normalization_and_local_targets.py
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

## Latest local arithmetic certificates

[Complete away-from-13 targets](COMPLETE_AWAY_TARGETS.md) proves dyadic vanishing and the exact 3-adic local sets for the two sixth-mode covers. [Certified infinity projection](INFINITY_PROJECTION.md) supplies precision-three de Rham and holomorphic normalization constants. Both have independent review and executable certificates. The subsequent sixth-mode proof above now completes moving QC functions, root isolation and the candidate sieve for the necessary center lifts. The original distinct-integer problem remains unresolved.

