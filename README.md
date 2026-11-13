# Ideal Prouhet–Tarry–Escott Solutions of Size 11

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
The remaining fifteen types' rational-point sets and lifting conditions
have not been completely determined.

[The exact bounded search in §5.33](PROOF.md#533-exact-bounded-center-search-in-four-remaining-types)
excludes types 4, 5, 15 and 16 when their specified rational center
`n/d` has `0<=n<=1000000`, `1<=d<=1000000` in lowest terms. The C++
multithreaded sieve retains only six colliding configurations, checked
independently with rational arithmetic. This finite box does not
classify the unbounded rational-point sets.

[Theorem 16, §5.42](PROOF.md#542-complete-exclusion-of-type-5-by-normalized-two-cover-descent)
also excludes **type 5 at every rational height**. An elementary
scalar-square lemma repairs the fake-cover ambiguity using the original
center conditions. Two rank-one elliptic Chabauty computations, over
Q(sqrt(505)) and Q(sqrt(-95)), exhaust both normalized branches and
leave only collisions or a negative center-square sum. All class-group
proof bounds use actual Minkowski bounds. Thus types 5, 10 and 16 are
completely excluded; the remaining fifteen types and the original
size-eleven problem remain open.

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
types; the separate individual-square argument above excludes type 16.
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
proves that the seven nonzero Selmer twisted Kummers admit no genus-one
fibration over Q, any finite odd-degree extension, or any quadratic
extension. This excludes that shortcut; it does not exclude rational
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
