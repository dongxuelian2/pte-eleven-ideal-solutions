# Ideal Prouhet–Tarry–Escott Solutions of Size 11

Fourteen of the eighteen normalized types in the symmetric shift construction
are excluded at every rational height. In type 17, elliptic Chabauty excludes
the odd multiples. On the even branch, two descents put every candidate on
an elliptic curve over a quartic field. At its unique dyadic place the
direction required by the chain lies outside the global 2-Selmer image.
Types 6, 7, 8 and 18 remain open.
[Theorem 32, §5.60](PROOF.md#560-complete-exclusion-of-type-17-by-a-dyadic-selmer-obstruction)

In type 12, both parity classes lead to rank-one elliptic curves over a
cubic field. Elliptic Chabauty at 907 forces the squared center to be
1 or 9. Each value makes the two chains share a vertex.
[Theorem 34, §5.62](PROOF.md#562-complete-exclusion-of-type-12-by-two-rank-one-covers)

For type 1 the two disks modulo 73 contain only the collision centers
0 and 1. A rank-one elliptic cover treats 1. A descent over the sextic
discriminant field and a second division at two dyadic places treat 0.
[Theorem 37, §5.65](PROOF.md#565-complete-exclusion-of-type-1-by-descent-over-the-discriminant-field)

In type 3, two squareclass pairs survive local tests. A rank-zero elliptic
curve eliminates one. For the other, elliptic Chabauty over the field of
discriminant −15 leaves one finite parameter, which makes the two input sets
share a vertex.
[Theorem 30, §5.58](PROOF.md#558-complete-exclusion-of-type-3-by-elliptic-chabauty)

For type 9, comparing three collinear points on a rank-one elliptic curve
determines all split fibers. The associated genus-four curve of ordered
roots has exactly nine rational points.
[Theorem 28, §5.56](PROOF.md#556-complete-exclusion-of-type-9-by-the-heights-of-three-collinear-points)

For fifteen types, including every open type, each finite set of primes
admits simultaneous local lifts with a common rational center of unbounded
height. [Proposition 29, §5.57](PROOF.md#557-simultaneous-local-lifts-in-fifteen-types)

In §5.32 all five nontrivial
doubled B centers in type 18 have a common reduced denominator D, supported
at primes `1 mod 10`. Each positive denominator valuation at those primes
is locally possible. A [certificate of 521 rigid primes](TYPE18_RIGID_PRIMES.md)
forces `D > 2.07 × 10^153`.

The original target remains open: construct two disjoint sets of eleven
distinct integers with equal power sums through degree ten, or prove that
none exist. Their monic root polynomials must differ by a nonzero constant:

\[
\prod_{a\in A}(x-a)-\prod_{b\in B}(x-b)=C\ne0.
\]

## Complete finiteness theorem for the symmetric shift construction

[Theorem 14, PROOF.md §5.32](PROOF.md#532-finiteness-for-all-eighteen-symmetric-shift-chain-types)
proves that **all eighteen** normalized chain types in the symmetric
size-ten shift-lift construction have only finitely many valid rational
configurations. This completes the twelve types previously left untreated.
Every type maps with finite fibers to an explicit smooth curve of genus
2, 3, or 4. Exact moment elimination, exceptional-branch analysis and
genus certificates establish these maps; Faltings's theorem supplies
finiteness.

Consequently this entire construction can produce only finitely many
successful affine-equivalence classes. An infinite family of inequivalent
size-eleven solutions cannot come from this route. The theorem is
ineffective: it supplies no rational-point list or height bound, and does
not show that the finite set is empty. It assumes neither GRH nor a
conjecture about a Tate–Shafarevich group.

[The sixth-mode exclusion](SIXTH_MODE_QC_EXCLUSION.md) separately eliminates
one type completely. The historical “sixth mode” is **type 10** in the
new exhaustive table. Its certified quadratic-Chabauty calculation isolates
48 local roots; the exact elliptic-image sieve leaves only eight known
colliding orbit points. The upstream arithmetic dependencies are listed
in that proof.

[Theorem 15, §5.40](PROOF.md#540-complete-exclusion-of-type-16-by-an-individual-square-condition)
also excludes **type 16 at every rational height**. Keeping two individual
square conditions gives a quartic birational to
`Y²=X³+44X²+100X`. An explicit two-isogeny descent proves rank zero;
good reduction determines its four rational points. Back-substitution
leaves only an irrational-center obstruction or repeated entries.
The proof and SymPy checker need no external arithmetic service or GRH.

[The exact bounded search in §5.33](PROOF.md#533-exact-bounded-center-search-in-four-remaining-types)
excludes types 4, 5, 15 and 16 when their specified rational center
`n/d` has `0<=n<=1000000`, `1<=d<=1000000` in lowest terms. The C++
multithreaded sieve retains only six colliding configurations, checked
independently with rational arithmetic. This finite box does not
classify the unbounded rational-point sets; the later theorems below
exclude these four types at all heights by separate arguments.

[Theorem 16, §5.42](PROOF.md#542-complete-exclusion-of-type-5-by-normalized-two-cover-descent)
also excludes **type 5 at every rational height**. An elementary
scalar-square lemma repairs the fake-cover ambiguity using the original
center conditions. Two rank-one elliptic Chabauty computations, over
Q(sqrt(505)) and Q(sqrt(-95)), exhaust both normalized branches and
leave only collisions or a negative center-square sum. All class-group
proof bounds use actual Minkowski bounds.

[Theorem 17, §5.43](PROOF.md#543-complete-exclusion-of-type-15-by-a-quadratic-unit-normalization)
excludes **type 15 at every rational height**. An elementary ideal and
unit argument in Q(sqrt(15)) restricts all genuine configurations to two
explicit squareclasses. Two rank-one elliptic Chabauty certificates
exhaust their rational images, leaving only nonsquare center parameters
or collisions.

[Theorem 18, §5.44](PROOF.md#544-complete-exclusion-of-type-4-by-individual-component-norms)
excludes **type 4 at every rational height**. Individual component norms
remove two of three fake two-cover classes. Local scalar-square
normalization and elliptic Chabauty over Q(sqrt(5)) exhaust the last
class and leave only a repeated-entry seed. The certificate transports
the original center coordinate through the exact minimal-model
isomorphism and uses unconditional class-group proof bounds.

[Theorem 19, §5.45](PROOF.md#545-complete-exclusion-of-type-2-by-an-elementary-elliptic-descent)
excludes **type 2 at every rational height**. An elementary two-isogeny
descent proves rank zero for its elliptic quotient, and good reduction
exhausts its ten rational points. Each point forces an irrational center
or a collision. Its independent SymPy checker needs no arithmetic service
or GRH.

[Theorem 20, §5.46](PROOF.md#546-complete-exclusion-of-type-14-by-chabauty-at-seven)
excludes **type 14 at every rational height**. A certified rank-one bound
and Chabauty at seven give exactly `∞` and `(9,±3240)` on its necessary
genus-two curve. The finite value `b=9` makes the two sets coincide.

[Theorem 23, §5.49](PROOF.md#549-complete-exclusion-of-type-11-by-chabauty-at-37)
excludes **type 11 at every rational height**. The squareclass
calculation leaves three relevant residue disks at 37 on a genus-three
cover. Two independent divisor logarithms prove rank two for its
Jacobian and give an annihilating differential that is a unit in all
three disks.
The three rational points there force a collision or a zero chain center.
This completes the exclusion of type 11.

[Proposition 22, §5.48](PROOF.md#548-exact-product-rank-and-two-cubic-field-covers-for-type-13)
proves rank two for the type-13 product Jacobian and identifies the
only two Kummer classes compatible with a square center. The resulting
genus-one covers are defined over a cubic field. One has certified rank
upper bound one.

[Propositions 24 and 25, §§5.50–5.51](PROOF.md#550-complete-exclusion-of-the-type-13-alpha-class)
exclude the entire alpha class and the part of the beta class where
the length-two center has 2-adic valuation at least two. Formal
logarithms and exact local halves give Selmer obstructions at two.
This reduces the problem to the beta class with center valuation one.
[Sections 5.52–5.53](PROOF.md#552-two-odd-prime-restrictions-and-a-dyadic-local-family-for-type-13)
restrict that class to the beta residue at 11 and 23 and give an explicit
genus-two Prym from the three individual singleton squares.
[Propositions 26 and 27, §§5.54–5.55](PROOF.md#554-six-independent-prym-classes-forced-by-a-type-13-lift)
exclude this last branch. A valid lift would force six independent Prym
classes. An unconditional class-group certificate for the degree-twelve
field makes the two-descent rank upper bound five. The exact Prym rank
lies between three and five.

Sections 5.37 and 5.39 combine software-assisted two-Selmer upper bounds
with exact Kummer certificates to give ranks **3,3,4,2** for types
**4,5,15,16**, respectively, and Sha[2]=0 for all four. The new type-4
certificate consists of four explicit rational degree-two divisors,
an existing rational point, and a five-by-five character matrix of
determinant one. Its discovery used Fisher–Yan's third-party software;
the lower-bound certificate is independently checked by SymPy.
Section 5.36 proves rational Néron–Severi rank one for all four
Jacobians, so both direct Chabauty rank criteria fail for all four
necessary curves. These rank obstructions alone do not exclude the chain
types; the separate individual-square and normalized-cover arguments
above exclude all four.
Reproduce the upper bounds, with unconditional quadratic class bounds, using
`scripts/check_shift_curve_selmer.m` in Magma; no rational-point
completeness or saturated basis is asserted.

Section 5.38 proves that the type-4 rational Richelot/double-Richelot
graph has exactly two Jacobian vertices, with no split quotient or
additional compatible 2-power-isogenous Jacobian. An exact Galois-action
argument closes this route to further models. The separate rank-three
result does not classify the rational points. The graph proof has
independent SymPy and Magma checks.

## Splitting-pencil arithmetic

[PENCIL_SELMER.md](PENCIL_SELMER.md) studies an explicit quintic–sextic
splitting pencil whose completely rationally split nondegenerate fiber
would give an original distinct-entry size-eleven solution. The necessary
genus-two discriminant curve has unconditionally
`dim Sel²(J/Q)=3`, `1≤rank J(Q)≤3`, and trivial rational torsion. Exact
principal-ideal relations and a complete Minkowski-bound verification
remove the former GRH hypothesis.

[The Kummer obstruction](kummer_obstruction/KUMMER_FIBRATION_OBSTRUCTION.md)
proves that the least extension degree permitting a genus-one fibration
over a split projective line on each of the seven nonzero Selmer twisted
Kummers is ten. All odd-degree extensions are still excluded. The result does not exclude rational
points or determine the Jacobian rank. A nonzero Cassels–Tate pairing or
two further independent rational divisor classes would distinguish the
current rank-one and rank-three routes. Neither is supplied here.

## Other proved restrictions and their scope

[PROOF.md §3](PROOF.md#3-common-center-pair-obstructions) contains the
imported common-center obstructions for symmetric distinct-entry pairs,
including the conditional five-pair height bound and forced divisors in
the three- and four-pair cases. Their stated hypotheses and literature
inputs remain essential.

[PROOF.md §4](PROOF.md#4-general-multiset-restrictions) concerns the multiset
relaxation. It proves multiplicity at most two, total support at least
fourteen, symmetric support at least eight per side, and explicit
Pell/Jacobian and simultaneous-splitting conditions at the boundary.
Repeated-root constructions in that section do not solve the original
distinct-entry problem.

## Necessary next mathematical work

The complete shift-route theorem replaces an unfinished type classification
with fixed higher-genus rational-point problems. An effective classification
must compute their rational points and retain only points satisfying all
center-square, complete-splitting and noncollision conditions. A square
discriminant alone is insufficient. The pencil route instead requires a
sharper Mordell–Weil result followed by a complete point computation.
Neither route currently gives an original size-eleven witness or a general
nonexistence proof.

## Exact verification

The all-type certificate rebuilds the polynomial identities and genus
checks from the four even moments. It is not a bounded point search.
Python 3 and the existing SymPy package suffice:

```sh
python scripts/check_all_shift_chain_types.py
python scripts/check_type11_thirty_seven_adic.py
python scripts/check_type4_richelot_graph.py
python scripts/check_type4_rank_three.py
python scripts/check_sixth_mode_qc_roots13.py
python scripts/check_sixth_mode_candidate_sieve.py
python scripts/pencil_selmer/verify.py
```

The pencil verifier includes recorded successful Magma arithmetic receipts;
it does not independently reimplement Magma's number-field and ramified
local arithmetic. The individual proofs identify their exact dependencies.
Research arguments and full reconstruction details remain in
[PROOF.md](PROOF.md); generated exploratory files belong in ignored `build/`.
