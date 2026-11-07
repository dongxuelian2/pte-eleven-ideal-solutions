# Working Mathematics

## 1. Problem

Seek integer multisets \(A,B\) of size 11, not equal as multisets, with equal power sums of orders 1 through 10. State explicitly if a subproblem additionally requires all entries on each side to be distinct.

Equivalently, for monic polynomials
\[
P_A(x)=\prod_{a\in A}(x-a), \qquad P_B(x)=\prod_{b\in B}(x-b),
\]
Newton's identities give \(P_A(x)-P_B(x)=C\), a constant, with \(C\neq0\) for a nontrivial pair.

Translation and nonzero integer scaling preserve the equal-power-sum property; record normalization carefully.

## 2. Principal target

Settle existence of an ideal size-11 integer PTE pair, or obtain a genuinely new general obstruction. If the general problem resists solution, prioritize proved structural results with clear hypotheses.

### Routes

- Study factorization and discriminant constraints implied by two split monic degree-11 polynomials differing by a nonzero constant.
- Develop prime-adic and modular restrictions; distinguish necessary conditions from sufficient ones.
- Reexamine common-center paired constructions and the known conditional obstructions; independently rederive any imported claim.
- Use symbolic cancellation identities and constrained searches to find candidate patterns, then verify exactly.
- Compare with prior PTE literature before claiming novelty or a new record.

## 3. Common-center pair obstructions

These proofs and exact certificates are transferred from
[ddy314/pte-order11-research, draft PR #1](https://github.com/ddy314/pte-order11-research/pull/1)
([source revision](https://github.com/ddy314/pte-order11-research/tree/011cf4b8681817bacc1f889b1f8fd8cf9b83123d)).
This is a consolidation of that work, not a new result. The source's
literature attribution and theorem hypotheses are retained. Unlike the
general multiset target above, every result below assumes distinct entries
and symmetry B = -A after an integer translation, together with the stated
number of common-center pairs.

### 3.1. Problem and status

For this section only, require disjoint sets A,B of eleven distinct integers such that

\[
\sum_{a\in A}a^j=\sum_{b\in B}b^j \qquad(1\le j\le10).
\]

Write f_A(x)=∏_{a∈A}(x−a), and similarly for B. Newton's identities give

\[
f_A(x)-f_B(x)=C\ne0.
\]

Here we study the **symmetric subclass** B=−A. Since the degree is odd,

\[
f_B(x)=-f_A(-x),\qquad
f_A(x)+f_A(-x)=C=2f_A(0)=-2\prod_{a\in A}a.\tag{1}
\]

The height is H=max_{a∈A}|a| in these centered integer coordinates. If a
solution is instead presented as B=2h−A, then h is an integer: reflection
makes 2h integral, and the first moment gives ∑A=11h integral. Thus an integer
translation puts it in the form above. The height is then max|a−h|.

A **common-center pair** is a pair of different entries {c−u,c+u} in A.
Several pairs are disjoint as sets of entries and use the same c. For integer
entries, c and u lie in Z[1/2]. This terminology is internal to A; it must not
be confused with the reflection B=−A.

The complete size-eleven problem remains open in the literature examined.
Coppersmith, Mossinghoff, Scheinerman, and VanderKam [1, §4.1 and Table 2]
exclude symmetric integer solutions only through height 3500. Their
quadratic-field near-solutions are not integer solutions. The published
constant-divisor baseline is recorded in their Table 1. A November 2026
preprint [2] claims a stronger symmetric 2-adic bound and local sharpness;
it explicitly reports no professional mathematical review. Nothing here
uses that claim. Neither these references nor our calculations justify
an exclusion of arbitrary size-eleven solutions.

The results below are independently derived conditional statements. A
novelty claim, a publication claim, and a complete classification of their
remaining rational points are deliberately not made.

### 3.2. Five pairs: an exact one-parameter normal form

**Theorem 1.** Suppose A consists of five common-center pairs and a singleton:

\[
A=\{s,c-u_1,c+u_1,\ldots,c-u_5,c+u_5\},
\]

and (A,−A) is an ideal size-eleven solution. Then c≠0, s=−10c, and for some
λ∈Z[1/2] with λ≠0,

\[
f_A(x)=(x+10c)\left[
\prod_{j=1}^{5}\big((x-c)^2-(2j-1)^2c^2\big)+\lambda
\right].\tag{2}
\]

In particular,

\[
C=20c\lambda.\tag{3}
\]

Here λ is the unnormalized constant in (2). If μ=λ/c^10 is used after
scaling the roots by c, then

\[
c^{-11}f_A(cx)=(x+10)(G(x-1)+\mu),\qquad C=20c^{11}\mu,
\]

where G is defined in (4). These two parameter conventions must not be mixed.

Conversely, every polynomial (2) satisfies f(x)+f(−x)=20cλ. If all its roots
are distinct integers and cλ≠0, its roots and their negatives form a valid
symmetric ideal solution. Integer splitting is an essential remaining
condition; the formula alone does not construct a solution.

**Proof.** The first moment gives s=−10c. If c=0, the five internal opposite
pairs already make A intersect −A, so c≠0. Put

\[
E(y)=\prod_{i=1}^{5}(y-u_i^2),\qquad
E_0(y)=\prod_{j=1}^{5}(y-(2j-1)^2c^2).
\]

Both are monic of degree five. The polynomial

\[
(x+10c)E_0((x-c)^2)=\prod_{k=-5}^{5}(x-2kc)
\]

is odd. Therefore (1) says that the even part of
(x+10c)(E−E_0)((x−c)^2) is constant.

If E−E_0 has degree r with 1≤r≤4 and leading coefficient d≠0, the coefficient
of x^(2r) in this last polynomial is (10−2r)cd≠0. Lower-degree terms cannot
contribute to that coefficient. This is a contradiction. Hence E−E_0=λ is
constant, which proves (2). Its coefficients lie in Z[1/2]. Evaluating at
zero, one factor of E_0(c²) is zero, so f_A(0)=10cλ, proving (3).

Conversely, the grid product is odd and the remaining summand is λ(x+10c),
so the constant identity holds. If a root a and its negative were both
roots, the constant identity would force C=0. Thus nonzero C supplies the
required disjointness once the eleven roots are distinct integers. ∎

**Important boundary.** The moment equations do not force λ=0. They leave
one free constant, and a real-root argument cannot simply eliminate it.
The obstruction below addresses integer splitting.

### 3.3. Exact finite-field fibers and a global conditional height bound

Define the fixed polynomial

\[
G(t)=\prod_{k\in\{1,3,5,7,9\}}(t^2-k^2)\in\mathbb Z[t].\tag{4}
\]

For a prime p>5 let L_p be the set of nonzero λ∈F_p for which G(t)+λ splits
completely into linear factors in F_p[t], **with multiplicities allowed**.

**Lemma 2 (local obstruction).** If L_p is empty, then p divides the integer
constant C of every solution covered by Theorem 1.

**Proof.** Work in Z[1/2] reduced modulo p. If c≡0, any one of the five
internal pairs reduces to {a,−a}. Thus f_A and f_{−A} share a root over F_p;
their difference C is zero modulo p.

If c≠0 modulo p, substitute x=c(t+1) in (2). After removal of the linear
factor and division by c^10, complete splitting of f_A implies complete
splitting of G(t)+λ/c^10. If C≠0 modulo p, equation (3) and p>5 imply
λ/c^10≠0. This contradicts L_p=∅. ∎

#### Finite certificate

For every prime 7≤p≤5000, `scripts/check_five_pairs_mod.py` computes L_p exactly.
The nonempty cases in the smaller, easily displayed range p≤1000 are:

| p | L_p |
|---|---|
| 11 | {1} |
| 349 | {283} |
| 353 | {332} |
| 601 | {139,215} |
| 647 | {307} |
| 673 | {501} |
| 719 | {686} |
| 727 | {146} |
| 743 | {392} |
| 787 | {157,302} |
| 809 | {491} |
| 821 | {110} |
| 853 | {353} |
| 857 | {449} |
| 877 | {25} |
| 919 | {683} |

Every other prime between 7 and 1000 has empty L_p, giving 149 obstruction
primes in that smaller interval. In particular, every prime from 13 through
347 is an obstruction. Extending the same exhaustive calculation through
5000 gives 393 obstruction primes. No conjecture or extrapolation beyond
the tested prime range is used.

This finite computation is exhaustive for each field: for every t∈F_p,
there is exactly one possible fiber λ=−G(t). The multiplicity of t in that
fiber is the least j≥1 for which the coefficient of X^j in G(t+X) is nonzero.
It is computed by

\[
\sum_{k=j}^{10}\binom{k}{j}g_k t^{k-j}\pmod p,
\quad G(X)=\sum_{k=0}^{10}g_kX^k.\tag{5}
\]

The expression with j=10 equals 1, so the multiplicity is always found.
Accumulating these multiplicities over a fiber gives 10 if and only if the
monic degree-ten polynomial splits completely. Nonzero λ only are retained.
Thus no root, repeated root, or fiber is omitted. The calculation uses only
integer arithmetic. An independent algorithm groups points by their fiber values, then removes
each possible linear factor by synthetic division, repeatedly where needed.
The two methods agree at every prime through 5000. The second method does not
use the Taylor multiplicity computation.

For a positive control, at p=349 and λ=283 the roots of G(t)+λ, counting
multiplicity, are

\[
0,0,14,30,96,132,217,253,319,335.
\]

The double root at zero is valid for the local test. This witness also
shows why the finite-field obstruction must not be declared universal.

Let

\[
S=\{p:7\le p\le5000,\ p\text{ prime},\ L_p=\varnothing\},
\qquad D=\prod_{p\in S}p.\tag{6}
\]

The integer D has 1195 decimal digits. Definitions (4)–(6) and the exhaustive
certificate specify it exactly, without a large decimal transcription. The
certificate asserts both the count |S|=393 and the height calculation below.

**Theorem 3 (conditional global height bound).** Every symmetric integer
ideal size-eleven solution having five common-center pairs on one side
satisfies

\[
H\ge\left\lceil D^{1/11}\right\rceil
=4264914784302653186685253811259318130882488669797076480249850785413789679541836572698823410736766760419030769.
\tag{7}
\]

**Proof.** By Lemma 2 and the finite certificate, D|C. All primes in D are
odd, so (1) gives D|∏a. This product is a nonzero integer, whence
D≤|∏a|≤H^11. The exact integer ceiling in (7) is checked by comparing its
11th power and its predecessor's 11th power with D. Its size is approximately
4.2649 × 10^108. ∎

The finite certificate concerns all residues of parameters modulo finitely
many primes. Its consequence (7) concerns all integer parameters of the
specified ansatz. It is not an exhaustive integer search and does not
assert global nonexistence. No primitivity assumption is needed.

### 3.4. Four pairs: a two-parameter normal form

Suppose eight entries of A form four common-center pairs, leaving three
singletons. Again c≠0. Normalize x by c over Q. Then

\[
F(x)=c^{-11}f_A(cx)=Q((x-1)^2)P(x-1),
\]

where P(t)=t³+11t²+bt+d and
Q(z)=z⁴+q_3z³+q_2z²+q_1z+q_0. The coefficient 11 in P follows from the
vanishing coefficient of x^10 in F. The coefficients of x^8,x^6,x^4,x²
then vanish exactly when

\[
\begin{aligned}
q_3&=(9b-d-330)/2,\\
q_2&=((7b-d-224)q_3+84b-28d-1848)/4,\\
q_1&=((5b-d-130)q_2+(35b-15d-644)q_3+126b-70d-1980)/6,\\
q_0&=((3b-d-56)q_1+(10b-6d-144)q_2\\
&\qquad +(21b-15d-272)q_3+36b-28d-440)/8.
\end{aligned}\tag{8}
\]

These equations are triangular and have a unique solution for Q. Direct
substitution also gives the factored constant

\[
F(0)=\frac{P(-1)P(-3)P(-5)P(-7)P(-9)}{384}\\
=-\frac{(b-d-10)(3b-d-72)(5b-d-150)(7b-d-196)(9b-d-162)}{384}.
\tag{9}
\]

Equations (8)–(9) are polynomial identities over Q; they can be checked
by multiplying the displayed cubic and quartic composition and comparing
coefficients. They remain valid modulo every prime p>3.

This gives an exact rational-point formulation of the four-pair ansatz:
P must have three distinct rational roots, Q must have four distinct
nonzero rational-square roots, the two factors must have no common root,
and (9) must be nonzero. Conversely these conditions give eleven distinct
rational roots and the constant identity; clearing denominators gives an
integer ideal solution. No such rational point is exhibited here.

**Theorem 4 (additional conditional divisors).** If a symmetric integer
size-eleven solution contains four common-center pairs in A, then

\[
23\cdot29\cdot37\cdot41\cdot43\cdot47\cdot53\cdot59\cdot67\cdot73\mid C.
\tag{10}
\]

**Proof.** For any of the displayed primes, if c≡0 modulo p the shared-root
argument of Lemma 2 applies. Otherwise normalize c=1 modulo p. The cubic P
and the octic Q(t²) must both split, allowing repetitions. Formula (9) must
be nonzero if C is nonzero. `scripts/check_four_pairs_mod.py` exhausts every pair
(b,d)∈F_p², constructs Q by (8), and tests both factors by repeated synthetic
division by every linear factor. There are zero survivors at each displayed
prime. Hence C cannot be nonzero modulo any of them. ∎

The complete counts through p=101 are printed by the certificate. At p=61
there is exactly one surviving pair (b,d)=(58,16), another positive control
against an unwarranted all-prime conclusion. A second check enumerates all
split cubics by their root multisets and tests the octic with Taylor
multiplicities, instead of enumerating coefficients and dividing factors.

### 3.5. Three pairs already force 29

This section uses the completeness of the local classifications in [1, §3]
as an explicit literature input. The preceding certificates for Theorems
3 and 4 do not depend on that input.

For a multiset with multiplicity function m on F_p, the maximum number of
disjoint pairs about a fixed center c is

\[
\frac12\sum_{x\ne c}\min\{m(x),m(2c-x)\}
+\left\lfloor\frac{m(c)}2\right\rfloor.\tag{11}
\]

The floor term is necessary because distinct integer entries may coincide
with the center after reduction.

For nonzero C, [1] gives, up to nonzero scaling, two symmetric local
size-eleven patterns modulo 23:

- {1,2,3,4,6,8,9,12,13,16,18};
- {1,2,2,3,3,5,5,7,7,17,17};

and one modulo 29:

- {1,1,1,15,19,21,21,22,23,25,25}.

Evaluating (11) at every nonzero center gives maxima 3, 3, and 2,
respectively. The short `pair_maximum` routine reproduces these counts.
Nonzero scaling preserves the maximum. A zero center would already force
C≡0. It follows that four common-center pairs force 23|C, while **three
common-center pairs already force 29|C**.

### 3.6. What remains

The general integer size-eleven existence problem is not solved. Even in
the four-pair and five-pair subclasses, these results give necessary
conditions, not a full exclusion or an integer construction. The exact
forms (2) and (8) leave rational splitting problems that require further
arithmetic work. The large bound (7) applies only to the five-pair subclass.

The useful outcome of this round is a proved normal form, explicit
unbounded structural restrictions, and an effective conditional lower
bound supported by small exact finite-field certificates. There is no
large integer search, no floating-point candidate, and no claim that
finite computations alone settle the unrestricted problem.

## 4. Open boundary

The general size-11 existence problem remains the target; restrictions to symmetric constructions do not settle it.

## References

[1] D. Coppersmith, M. J. Mossinghoff, D. Scheinerman, and J. M. VanderKam,
*Ideal solutions in the Prouhet–Tarry–Escott problem*, Mathematics of
Computation 93 (2024), 2473–2501.
[DOI](https://doi.org/10.1090/mcom/3917),
[arXiv version](https://arxiv.org/abs/2304.11254),
[full text](https://arxiv.org/html/2304.11254v1).

[2] D. Niedbala Giraudin, *The 2-adic valuation of the constant for symmetric
ideal Prouhet–Tarry–Escott solutions*, Hexagon preprint 2611.00057v1,
November 4, 2026.
[Record and provenance](https://hexagonmath.org/2611.00057v1).
This is an explicitly AI-generated, not professionally reviewed preprint;
its claimed bounds are not used in our proofs.
