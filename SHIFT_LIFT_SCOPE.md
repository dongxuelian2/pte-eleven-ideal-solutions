# What the shift exclusion establishes

Theorem 49 completes the eighteen-type exclusion. Two consequences give
different ways to state the result. No ten-element distinct-entry precursor
can produce a symmetric ideal size-eleven pair by one rational shift. Also,
a rational arithmetic progression of length ten is determined among rational
multisets of that size by its first nine moments. The second theorem has no
lift or symmetry hypothesis on the competing multiset.
[Proofs of Theorem 49 and Corollaries 50–51](PROOF.md#576-complete-exclusion-of-type-18-and-rigidity-of-a-ten-term-progression)

## Why try this construction?

Choudhry's *A new approach to the Tarry–Escott problem*, §2, Lemma 1 and the
examples on manuscript pages 3–4, uses

\[
 A\cup(B+h)\;=_{k+1}\;B\cup(A+h)
\]

when \(A=_k B\), followed by cancellation of common entries. His arithmetic
progressions are chosen so that long interiors cancel and only endpoints
remain. The paper applies this method to obtain parametric ideal solutions of
degrees four through seven. These successes explain why trying a more demanding
cancellation at size eleven is mathematically motivated. The paper
allows longer precursors and repeated applications of the lemma. It supplies
no assertion that a single symmetric ideal size-ten precursor should succeed
at size eleven. [Choudhry, arXiv:1603.00206v2, §2](https://arxiv.org/pdf/1603.00206)

Wright gives a sequence of shifts immediately after equation (11). Its
first four lifts take size two successively to ideal sizes three, four,
five and six. For the step from three to four terms, start
with the ideal size-three pair \(\{0,4,5\}=_2\{1,2,6\}\) and use shift
\(h=5\). Cancelling the common entries 5 and 6 gives
\(\{0,4,7,11\}=_3\{1,2,9,10\}\), a distinct ideal size-four pair.
Both identities follow by direct expansion of the first two or three moments. [E. M. Wright, *Equal sums of like powers*, equation (11) and its example](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/2F4A6DC62111B56B3B6891C4EE316917/S000843950005462Xa.pdf/equal-sums-of-like-powers.pdf)

The size-eleven existence question has an independent literature. Coppersmith,
Mossinghoff, Scheinerman and VanderKam report integer ideal solutions for sizes
at most ten and size twelve, and symmetric searches through height 3500 at
size eleven. Their introduction also recalls the infinite family of symmetric
size-ten solutions obtained by Smyth using Letac's construction and an elliptic
curve. Availability of size-ten inputs and success of shift cancellation at
lower degrees support exploring the mechanism. They give no completeness claim
for it at size eleven. [Coppersmith et al., introduction and §4.1](https://arxiv.org/html/2304.11254v1)

The specific single-shift question at size eleven was selected in this
project. The checked sources establish the construction and its earlier
successes. They leave open which other mechanisms might yield a size-eleven
solution. The progression theorem has a separate formulation: the polynomial
\(T(T-1)\cdots(T-9)-c\) splits completely over the rationals exactly when
\(c=0\). Banks and van der Poorten give a related obstruction when the
*union* of the two sets is an arithmetic progression of more than four
terms. Corollary 50 instead fixes just one side to be a ten-term progression
and allows every rational multiset on the other side.
[Banks and van der Poorten, *Squares from products of integers*, §3](https://austms.org.au/wp-content/uploads/Gazette/2004/Mar04/Mar2004.pdf)

The arithmetic-progression rigidity is sensitive to size and to the
ground field. For four terms, the progression \(\{\pm1,\pm3\}\) and
\(\{\pm9/5,\pm13/5\}\) have the same first three moments, since
\(1^2+3^2=(9/5)^2+(13/5)^2\). The conic \(u^2+v^2=10\) supplies
infinitely many rational examples. For ten terms, moving the constant
slightly also preserves ten simple *real* roots, by continuity near
the ten distinct roots at zero shift. Corollary 50 singles out complete
rational splitting.

## A general inverse theorem

For a finitely supported integer-valued function \(u\) on \(\mathbf Q\), let
\(T_hu(x)=u(x-h)\). For a pair of disjoint finite sets write
\(\nu=\mathbf1_A-\mathbf1_B\). Equal moments through degree \(k\), including
cardinality at degree zero, mean that
\(\sum_x x^j\nu(x)=0\) for \(0\le j\le k\).

**Proposition.** Fix \(h\ne0\). There is at most one finitely supported
function \(u\) satisfying \(\nu=(1-T_h)u\). On each coset of
\(h\mathbf Z\), it exists precisely when the signed sum of \(\nu\) on that
coset is zero. For \(h>0\), it is given by the cumulative sums

\[
 u(r+jh)=\sum_{i\le j}\nu(r+ih).
\]

It is the signed indicator of two disjoint finite sets precisely when all
these cumulative sums belong to \(\{-1,0,1\}\). If \(\nu\) has zero moments
through degree \(k\ge1\), then \(u\) has zero moments through degree
\(k-1\). Its positive and negative masses are equal. For a signed
indicator this gives equal support sizes.

**Proof.** The equation gives the recurrence
\(u(r+jh)-u(r+(j-1)h)=\nu(r+jh)\). Finite support fixes the initial value
as zero. Its eventual value is the sum on that coset, which must vanish.
For moments put \(M_j=\sum_x x^ju(x)\). Then

\[
 \sum_x x^j\nu(x)
 =-\sum_{i=0}^{j-1}\binom ji h^{j-i}M_i.
\]

For \(j=1\) this gives \(M_0=0\). Induction using the nonzero coefficient
\(-jh\) gives \(M_{j-1}=0\). The statement about indicators and equal
cardinalities follows. ∎

The moment-divisibility aspect belongs to the classical generating-polynomial
formalism. Borwein, Lisonek and Percival, Lemma 1 on journal page 2065,
state its equivalence with divisibility of
\(\sum_{a\in A}z^a-\sum_{b\in B}z^b\) by a power of \(z-1\).
The proposition records the exact finite-support and unit-coefficient
conditions needed when reversing a construction. The support conditions follow directly from the recurrence.
[Borwein et al., *Math. Comp.* 72 (2003), 2063–2070](https://www.daemonology.net/papers/pte.pdf)

## Rational shifts into an integer ideal target

**Theorem.** Let \(A,B\subset\mathbf Z\) be a nontrivial distinct-entry
ideal PTE pair of size \(n+1\). Suppose it is obtained by cancelling common
terms from a single shift of disjoint rational sets \(X,Y\), each of size
\(n\). Then the shift is an integer and \(X,Y\subset\mathbf Z\).
Consequently, for any fixed integer target, all such precursors can be tested
using finitely many integer shifts. The statement allows nonsymmetric targets
and precursors.

**Proof.** Change the orientation to assume \(h=p/q>0\), with \(p,q\) coprime
positive integers. A rational translation orbit containing an integer is in
\(q^{-1}\mathbf Z\). Orbits not meeting the target have zero cumulative
inverse, so the precursor is supported in \(q^{-1}\mathbf Z\). Let

\[
 P(z)=\sum_{a\in A}z^a-\sum_{b\in B}z^b.
\]

Laurent polynomials suffice if an entry is negative. Finite support gives
\(U(z)=P(z^q)/(1-z^p)\). Since exponentiation by \(q\) permutes the
\(p\)-th roots of unity, \(1-z^p\) divides \(P(z)\) as well. Put
\(V(z)=P(z)/(1-z^p)\). Then

\[
 U(z)=\bigl(1+z^p+\cdots+z^{(q-1)p}\bigr)V(z^q).
\]

The \(q\) summands have disjoint supports modulo \(q\). Hence \(V\) has
coefficients in \(\{-1,0,1\}\), and its positive and negative supports each
have size \(n/q\). The inverse proposition makes these supports a nontrivial
PTE pair through degree \(n-1\). Newton's identities require at least \(n\)
entries per side for such a pair. Thus \(n/q\ge n\), so \(q=1\).
Finite-support uniqueness then gives integer precursor entries.

For an integer shift \(p\), balance in every residue class pairs target
entries on the two sides modulo \(p\). At least one nonzero cross difference
is a multiple of \(p\), so \(p\le\max(A\cup B)-\min(A\cup B)\). Testing
the cumulative inverse for each such positive integer is exact and finite. ∎

In fact, for a fixed numerator \(p\), passage from shift \(p\) to \(p/q\)
multiplies precursor cardinalities by \(q\). Increasing the denominator of the
shift cannot rescue an ideal-to-ideal single lift of a fixed integer target.
The theorem follows from the classical divisibility representation and the
minimal-cardinality bound. Its provenance should be checked further before
claiming priority for this formulation.

## Symmetry is forced backwards

For \(c\in\mathbf Q\), put \(R_cu(x)=u(2c-x)\).

**Theorem.** If \(\nu=(1-T_h)u\), with both functions finitely supported, and
\(R_c\nu=\epsilon\nu\), \(\epsilon\in\{1,-1\}\), then

\[
 R_{c-h/2}u=-\epsilon u.
\]

**Proof.** Reflection satisfies \(R_cT_h=T_{-h}R_c\), and
\(1-T_{-h}=-T_{-h}(1-T_h)\). Thus the reflection identity gives

\[
 (1-T_h)(R_{c-h/2}u+\epsilon u)=0.
\]

The finite-support kernel of \(1-T_h\) is zero, by the inverse proposition.
This proves the assertion. ∎

For a symmetric odd-size target \(B=2c-A\), one has \(R_c\nu=-\nu\).
Its precursor therefore has both sides separately symmetric about
\(c-h/2\). This holds even if no symmetry was imposed on the proposed input.
Disjoint symmetric sets with equal cardinalities have even cardinality:
if both had odd size, both would contain the reflection center.

**Corollary 51 of the main proof.** A distinct symmetric ideal
size-eleven pair cannot be obtained from any distinct ideal size-ten pair by
a single rational shift and cancellation. Every disjoint distinct precursor
for such a target has at least twelve entries per side. This includes a final
step inside a longer sequence of shifts whenever that step has a simple
precursor.

**Proof.** The inverse proposition forces equal moments through degree nine,
so a nontrivial precursor has at least ten entries per side. Backward symmetry
makes both precursor sides separately symmetric and makes their cardinality
even. At cardinality ten, §5.1 places the precursor in one of the eighteen excluded classes. The next possibility is twelve. ∎

This extends the input hypothesis in the construction statement and gives an
explicit lower bound on precursor size. It still concerns symmetric
size-eleven targets. General nonsymmetric size-ten inputs can only produce
nonsymmetric size-eleven targets.
The existence of those targets remains an independent problem.

## A transport interpretation

For an integer ideal pair, the classical interlacing theorem implies that
\(\#\{a\in A:a\le t\}-\#\{b\in B:b\le t\}\) lies in
\(\{-1,0,1\}\) for every real \(t\). One proof uses the common derivative
of the two monic root polynomials: its critical points divide the line into
monotonic intervals, each containing exactly one root from each polynomial.
Within successive intervals their order reverses. Thus the unit-shift inverse
always consists of disjoint sets with no repetitions.

If \(a_1<\cdots<a_N\) and \(b_1<\cdots<b_N\), the cardinality on each side
of that inverse is

\[
 r=\frac12\sum_{i=1}^N|a_i-b_i|.
\]

Indeed, the sum of the absolute cumulative differences is the integral
transport distance between the two ordered sets. Its positive and negative
sums agree because their first moments agree. For a symmetric size-eleven
target, Corollary 51 gives \(r\ge12\).
This is a necessary condition on every symmetric target in its given integer coordinates.

The interlacing theorem is recorded with proof as Theorem 2.1.3 in Caley's
thesis. The support and transport deductions here follow from it and the
inverse theorem.
[Caley, *The Prouhet–Tarry–Escott problem*, pp. 6–7](https://dspacemainprd01.lib.uwaterloo.ca/server/api/core/bitstreams/15979d95-3eb1-42e3-8101-ef89e4e9bbbe/content)

## Minimum ordered distance forces symmetry

**Theorem 52.** Let \(A,B\) be a disjoint pair of integer sets of size
\(n\ge2\) with equal power sums through degree \(n-1\). Order the
entries of each set increasingly and put \(L=\sum_{i=1}^n|a_i-b_i|\).
Then \(L\ge2(n-1)\). Equality forces standard PTE symmetry: a common
reflection exchanges the two sets for odd \(n\), and preserves each of
them for even \(n\). Equivalently, a nonsymmetric pair satisfies
\(L\ge2n\). For example, the size-four pair displayed above attains six.

**Proof.** Translate the minimum to zero, and let \(u\) be the
unit-shift inverse of the target. The common-derivative argument above
shows that the nonzero values of \(u\) form \(n\) successive blocks
of alternating signs. The blocks have lengths \(|a_i-b_i|\). A block
starts at the smaller of the paired roots and ends immediately before
the larger. Consecutive blocks have at least one unused integer between
them. Each sign occurs at \(L/2\) integers, and these two sets have
equal moments through degree \(n-2\). Newton's identities give
\(L/2\ge n-1\). At equality the inverse itself is ideal; its signed
cumulative count also stays between \(-1\) and \(1\). Each interior
block therefore has at most two entries, and the first and last blocks
have at most one each. Their upper bounds sum to \(2(n-1)\), so all
bounds are attained.

Write \(D\) for the largest inverse entry, and choose the sign so the
singleton at zero is positive. Let \(A_0(x),B_0(x)\) be the monic
polynomials with roots at the starting integers of the positive and
negative *interior* doubletons. Their roots are disjoint. A doubleton
beginning at \(s\) contributes \((x-s)(x-s-1)\) to the appropriate
root polynomial. Write \(R=B_0/A_0\). Since the inverse is ideal,
its two monic root polynomials differ by a constant.

For \(n=2m+1\), the endpoint \(D\) is positive. The two polynomials
are \(x(x-D)A_0(x)A_0(x-1)\) and \(B_0(x)B_0(x-1)\), with
\(\deg A_0=m-1,\ \deg B_0=m\). In Laurent series at infinity,
\[
 R(x)R(x-1)=x(x-D)+O(x^{-2m+2}).
\]
The function \(R^*(x)=-R(D-1-x)\) satisfies the same estimate, by
substitution of \(D-x\) in the product. Both start with \(x\).
If \(R-R^*\) has leading term \(c x^e\ne0\), its denominator has
degree at most \(2m-2\), so \(e\ge-2m+2\). The difference of the
two products then starts with \(2c x^{e+1}\), contradicting the
estimate. Thus \(R=R^*\). The zeros and poles of this reduced function
are each invariant under \(s\mapsto D-1-s\). Reflection about \(D/2\)
therefore preserves each color of inverse entries, including the two
endpoints.

For \(n=2m\), the singleton at \(D\) is negative. Now both starting
polynomials have degree \(m-1\), and the root polynomials are
\(x A_0(x)A_0(x-1)\) and \((x-D)B_0(x)B_0(x-1)\). Hence
\[
 R(x)R(x-1)=\frac{x}{x-D}+O(x^{-2m+1}).
\]
Here \(R^*(x)=1/R(D-1-x)\) satisfies the same estimate: substitute
\(D-x\) and take the reciprocal of a series with constant term one.
Both functions start with one. A nonzero difference with leading term
\(c x^e\) again has \(e\ge-2m+2\), whereas the product difference
starts with \(2c x^e\). This is impossible. Thus the starting points
of the two colors are reflected into each other. Reflection about
\(D/2\) exchanges the colors of the inverse, including the endpoints.

Taking the unit difference reverses the reflection sign and moves the
center to \((D+1)/2\). This gives the asserted symmetry of the target
in both cases. The distance sum is even, which gives the bound for a
nonsymmetric pair. ∎

**Corollary 53.** For *every* disjoint distinct-entry integer ideal
size-eleven pair, \(\sum_{i=1}^{11}|a_i-b_i|\ge22\). For a symmetric
pair the lower bound is 24. Both bounds scale by the common lattice
spacing if all the entries belong to a translate of \(d\mathbf Z\).

**Proof.** At equality in the elementary bound 20, Theorem 52 gives
a symmetric target with a ten-entry unit-shift predecessor. Corollary 51
excludes this. In fact only the Type 5 exclusion is needed here. After
centering the predecessor, its negative side consists of five doubletons,
one fixed by reflection and two reflected pairs of doubletons. Its other
side consists of four doubletons and the two end singletons. This is
exactly Type 5, with \((m,aa,bb)=(1,(2,2),(1,2,2))\). The ordered
distance sum is even, and for symmetric targets its half is even by
backward symmetry. Translation and division by \(d\) prove the lattice
version. ∎ [Complete exclusion of Type 5, §5.42](PROOF.md#542-complete-exclusion-of-type-5-by-normalized-two-cover-descent)

A focused paper could lead with the ten-term progression theorem. Its
statement only involves power sums or the completely split fibers of
an explicit polynomial. The proof uses the single vertex descent and
rank-one calculation and is independent of the other seventeen exclusions.
Theorem 52 gives a separate general mechanism, and Corollary 53 needs only the Type 5 exclusion in addition. A paper on the full single-shift obstruction would also need the other
seventeen proofs and their arithmetic certificates. The unrestricted
size-eleven existence question still allows other constructions and
nonsymmetric targets. Priority for the progression theorem and the inverse
statements, and independent review of the arithmetic certificates, remain
open.
