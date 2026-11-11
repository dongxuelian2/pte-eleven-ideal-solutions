# Working Mathematics

## 1. Problem

The principal task is to construct two disjoint sets \(A,B\subset\mathbb Z\),
each containing eleven distinct integers, with equal power sums of orders
1 through 10, or prove that no such pair exists.

Section 4 also studies the broader multiset relaxation, where repeated
entries are allowed. Its results on repeated-root patterns are ancillary:
they do not solve or exclude the principal distinct-entry problem.

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
literature attribution and theorem hypotheses are retained. These results
address the distinct-entry target only under the additional assumption
of symmetry B = -A after an integer translation, together with the stated
number of common-center pairs.

### 3.1. Problem and status

In this section, let A,B be disjoint sets of eleven distinct integers such that

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

## 4. General multiset restrictions

This section studies a multiset relaxation, allowing repeated entries.
All statements apply to rational multisets and therefore to integer multisets.
The multiplicity and support exclusions do not exclude the principal case
of eleven distinct entries per side, and the repeated-root constructions
here would not satisfy that task. They are ancillary structural results,
not a solution of the principal existence problem.
The extremal argument is the classical polynomial Chebyshev/Pell mechanism;
no novelty claim is made.

### 4.1. Repetitions and critical points

**Theorem 5.** Let A and B be nontrivial ideal size-eleven rational
multisets. Their supports are disjoint, no entry has multiplicity greater
than two, and their common derivative has ten distinct real zeros.
Consequently each side has at least six distinct entries.

**Proof.** Write P=P_A and Q=P_B=P-C, where C is nonzero. A common
root would give C=0. Both polynomials split over the reals.
For either root polynomial R and any real r outside its roots,

\[
\left(\frac{R'}R\right)'(r)
=-\sum_{z\in R^{-1}(0)}\frac{m_z}{(r-z)^2}<0,
\tag{12}
\]

where m_z is the root multiplicity. If P has a root r of multiplicity
at least three, then Q(r)=-C is nonzero but Q'(r)=Q''(r)=0,
contradicting (12). The same argument applies to Q.

Rolle's theorem, including multiplicities, shows that all ten zeros of
P' are real. At any such zero r, at least one of P(r),Q(r) is nonzero.
Apply (12) to that polynomial: since its first derivative is zero,
its second derivative is nonzero. Thus every zero of P' is simple.
The support bound follows from eleven entries of multiplicity at most two. ∎

### 4.2. The six/six case is impossible over the rationals

**Theorem 6.** In Theorem 5 the two supports cannot both have size six.
Thus the disjoint union of the supports has at least thirteen elements.
In particular, a symmetric rational solution B=2h-A must have at least
seven distinct entries on each side.

**Proof.** If each support has size six, the multiplicity bound forces
five double roots and one simple root on each side. Hence

\[
P(x)=(x-a)U(x)^2,\qquad Q(x)=(x-b)V(x)^2,
\tag{13}
\]

where a,b are rational and U,V are monic degree-five polynomials with
rational roots. Disjointness gives a≠b and gcd(U,V)=1.
Each of U,V divides P'=Q'. Degree and leading coefficient then give
P'=11UV, and therefore

\[
(x-a)(x-b)P'(x)^2=121P(x)(P(x)-C).
\tag{14}
\]

Set m=(a+b)/2, d=(b-a)/2, x=m+dt, and
F(t)=(2P(m+dt)-C)/C. Then d≠0, and (14) becomes

\[
(t^2-1)F'(t)^2=121(F(t)^2-1).
\tag{15}
\]

Differentiate this polynomial identity and cancel the nonzero polynomial
2F'. This gives

\[
(t^2-1)F''+tF'=121F.
\tag{16}
\]

The polynomials satisfying (16) form a one-dimensional space spanned by
the Chebyshev polynomial T_11. Indeed, T_11 satisfies (16), as follows
from T_11(cos θ)=cos(11θ). If another solution has degree at most eleven,
subtract a multiple of T_11 to cancel its leading coefficient. A nonzero
remainder of degree j<11 would have leading coefficient multiplied by
j²-121 in (16), a contradiction. Thus F=kT_11. Since x=b corresponds
to t=1 and P(b)=C, we have F(1)=1, so k=1.

The exact factorization

\[
T_{11}(t)+1=(t+1)
(32t^5-16t^4-32t^3+12t^2+6t-1)^2
\tag{17}
\]

now contradicts rational splitting of P. To see this without any
algebraic-number classification, let
h(t)=32t^5-16t^4-32t^3+12t^2+6t-1.
A rational root in lowest terms has denominator dividing 32, by the
rational-root theorem, and can therefore be reduced modulo 3. But

\[
h(0)\equiv2,\qquad h(1)\equiv1,\qquad h(2)\equiv1\pmod3.
\]

Thus h has no rational root. Since F+1=2P(m+dt)/C and the affine change
has rational coefficients with nonzero slope, rational splitting of P
would force h to split over the rationals, which is impossible.

Each support has at least six elements by Theorem 5, so excluding six/six
gives the bound thirteen. Reflection preserves support size; hence a
symmetric pair must have at least seven on each side. ∎

The rationality hypothesis in Theorem 6 is essential. The two real
fibers T_11(t)+1 and T_11(t)-1, after division by their common leading
coefficient 1024, are monic degree-eleven polynomials differing by a
nonzero constant. Both split over the reals, and each has precisely
six distinct roots. They furnish a real extremal example, not a rational
or integer PTE construction. Theorem 6 does not assert that the bound
thirteen is attained over the rationals.

### 4.3. An elliptic obstruction: fourteen values in total

We use the following established consequence of Mazur's torsion theorem
[3, Introduction, Theorem (8)]: an elliptic curve over Q has no rational
point of order 11. This is a theorem, not an additional conjectural
hypothesis. The polynomial Pell/divisor connection below is classical;
see [4, Section 2] for the general framework. We give the specific
argument needed here in full.

**Lemma 7 (quartic Pell obstruction).** Let D∈Q[x] be monic, squarefree,
and of degree four. There are no S,T∈Q[x] and k∈Q\{0} satisfying

\[
S^2-DT^2=k^2,\qquad \deg S=11.
\tag{18}
\]

**Proof.** Degree comparison gives deg T=9 and equality, up to sign,
of the leading coefficients of S and T. Change the sign of T if necessary
to make these coefficients equal, say to s≠0.
The smooth projective curve E with affine equation y²=D(x) has genus one.
It has two distinct rational points O,R at infinity, characterized by
y/x² tending to +1 and -1. More explicitly, t=1/x and w=y/x² give
w²=t⁴D(1/t), whose two points at t=0 are w=±1; they are smooth and t is
a local parameter at each.

The function f=(S+Ty)/k and its conjugate (S-Ty)/k are mutual inverses
and are regular on the affine curve. Thus f has no finite zero or pole.
At O its leading term is (2s/k)t^(-11), so it has a pole of exact order 11.
Conjugation exchanges O,R and sends f to its inverse. Consequently

\[
\operatorname{div}(f)=11R-11O.
\tag{19}
\]

Use O as the rational origin of E. The genus-one isomorphism
E→Pic⁰(E), X↦[X-O], identifies R with a nonzero rational point killed
by 11. Since 11 is prime and R≠O, its order is exactly 11, contradicting
Mazur's theorem. ∎

**Theorem 8.** The supports of any nontrivial ideal size-eleven rational
pair contain at least fourteen distinct values in total.

**Proof.** Theorem 6 excludes totals below thirteen. A total of thirteen
must have support sizes six and seven. Interchange the sides if necessary.
By Theorem 5 their root multiplicities force

\[
P=(x-a)U^2,\qquad Q=(x-b)(x-c)(x-d)V^2,
\quad \deg U=5,\quad\deg V=4.
\tag{20}
\]

The four simple roots a,b,c,d are distinct rational numbers, since the
supports are disjoint. Hence D=(x-a)(x-b)(x-c)(x-d) is a squarefree
monic quartic. With S=2P-C and T=2UV, the identity P-Q=C gives

\[
S^2-DT^2=(2P-C)^2-4PQ=C^2.
\tag{21}
\]

Here deg S=11, contradicting Lemma 7. ∎

### 4.4. Symmetry strengthens the bound to eight values per side

**Theorem 9.** If a nontrivial ideal size-eleven rational pair satisfies
B=2h-A, each side has at least eight distinct values. Equivalently,
there are at most three double entries on each side.

**Proof.** The first moment gives h=(∑A)/11∈Q, so translate by h and
assume B=-A. Write P=P_A and Q=P_B=-P(-x). Then P+P(-x)=C.
Theorem 6 already gives at least seven distinct values per side.
Suppose there are exactly seven. By Theorem 5,

\[
P(x)=(x-a)(x-b)(x-c)U(x)^2,\qquad \deg U=4.
\tag{22}
\]

Disjointness of A and -A ensures abc≠0 and that a²,b²,c² are pairwise
distinct. Define

\[
S(x)=2P(x)-C,\quad T(x)=2U(x)U(-x),\quad
D(x)=(x^2-a^2)(x^2-b^2)(x^2-c^2).
\]

Then S is odd of degree eleven, T is even of degree eight, and
S²-DT²=C². Write S=xR(x²), T=T₀(x²). The degrees of R,T₀ are five
and four, respectively, and both have leading coefficient 2.
Putting z=x² and D₀(z)=(z-a²)(z-b²)(z-c²) gives

\[
zR(z)^2-D_0(z)T_0(z)^2=C^2.
\tag{23}
\]

Now set

\[
D_1(z)=zD_0(z),\quad
A_1(z)=2zR(z)^2-C^2,\quad B_1(z)=2R(z)T_0(z).
\]

Direct multiplication, using (23), gives

\[
A_1^2-D_1B_1^2=C^4.
\tag{24}
\]

The polynomial D₁ is a squarefree monic quartic with distinct rational
roots 0,a²,b²,c². Moreover, deg A₁=11 and deg B₁=9, with leading
coefficients both equal to 8. Equation (24) contradicts Lemma 7 with
k=C². This excludes seven values per side and proves the result. ∎

### 4.5. The remaining arithmetic boundary

The elliptic obstruction does not extend automatically to higher-genus
Jacobians. To make the precise remaining problem explicit, let r be the
total support size of a rational solution. Let U,V be the monic products
of its double roots on the respective sides, and let D be the monic
product of all its simple roots across both sides. Then

\[
\deg D=2r-22,\quad\deg(UV)=22-r,\quad
(2P-C)^2-D(2UV)^2=C^2.
\tag{25}
\]

D is squarefree and completely split over Q. The smooth projective curve
y²=D has genus g=r-12. For the remaining range r≥14 this gives g≥2.
The same infinity-divisor calculation as in Lemma 7 gives a rational
divisor class [R-O] killed by 11 in its Jacobian. This class is nonzero:
a principal divisor R-O would give a degree-one map to the projective
line, impossible on a curve of positive genus. Thus its order is exactly
11. Mazur's theorem concerns elliptic curves and does not rule out such
a class in a genus-two or higher Jacobian.

For a symmetric pair with s distinct values per side, the same explicit
descent as in Section 4.4 gives a monic squarefree polynomial

\[
D_1(z)=z\prod_{i=1}^{2s-11}(z-a_i^2)
\]

and a degree-eleven Pell equation on a curve of genus s-6, where the a_i
are the simple roots on one side. In detail, deg R=5 and deg T₀=11-s
in (23), while deg A₁=11 and deg B₁=16-s in (24). At the first remaining
case s=8, D₁ has degree six and the obstruction again lies in a genus-two
Jacobian. Neither this torsion condition nor complete splitting of D₁
is proved sufficient for PTE splitting. No admissible genus-two point
or ideal integer pair is constructed here.

### 4.6. The genus-two Jacobian must be simple over Q

**Theorem 10.** In either of the following cases, the genus-two Jacobian
arising in Section 4.5 is simple over Q and its rational torsion contains

\[
(\mathbb Z/2\mathbb Z)^4\times\mathbb Z/11\mathbb Z
\simeq [2,2,2,22],
\tag{26}
\]

a subgroup of order 176:

- a rational ideal size-eleven pair with total support size fourteen;
- a symmetric rational ideal size-eleven pair with eight distinct entries
  per side.

Here simplicity over Q excludes an elliptic subvariety defined over Q;
it does not assert geometric simplicity.

**Proof.** In both cases Section 4.5 gives a monic squarefree sextic D
split completely over Q, and a rational class δ of exact order eleven
in J=Jac(y²=D). Its six Weierstrass points are rational. The fifteen
nonzero two-torsion classes are represented by W_i+W_j-K, for the fifteen
unordered pairs of distinct Weierstrass points and a canonical divisor K.
They are all rational, so J[2](Q)≅(Z/2Z)^4. This standard description is
also recorded in [5, Section 1.2]. The subgroup they generate with δ is
the direct product (26), because the orders are coprime.

If J contained an elliptic subvariety E defined over Q, the quotient
F=J/E would be an elliptic curve over Q with scheme-theoretic kernel E.
The image of δ in F(Q) is killed by eleven, so Mazur's theorem makes it
zero. Thus δ belongs to E(Q), where the same theorem makes δ zero,
a contradiction. This argument uses only exactness at J(Q), not
surjectivity onto F(Q). It also excludes every Q-isogeny to a product
of elliptic curves, regardless of its degree: such an isogeny implies
an elliptic subvariety over Q. ∎

**Corollary 11.** The genus-two curve in Theorem 10 cannot admit a
nonconstant map over Q to an elliptic curve. In particular, its six
branch values cannot form three opposite pairs about a rational center
in the given affine coordinate.

**Proof.** A nonconstant map induces a nonzero, hence surjective,
homomorphism from its Jacobian to the elliptic curve. The identity
component of the kernel is an elliptic subvariety over Q, contradicting
Theorem 10. For the stated branch configuration, translate the center
to zero and write

\[
y^2=(x^2-u^2)(x^2-v^2)(x^2-w^2),
\tag{27}
\]

where u,v,w are nonzero rational numbers with distinct squares.
The map (x,y)↦(x²,y) has target
Y²=(X-u²)(X-v²)(X-w²), a smooth elliptic curve over Q. ∎

This strengthens the conditions on the remaining boundary, but does not
raise the support bounds. In particular, a table of known genus-two
torsion groups is not an upper bound. The introduction of [5] explicitly
distinguishes its realized groups from a complete classification.

### 4.7. A sextic Pell counterexample that is not a PTE solution

The degree-eleven quartic obstruction in Lemma 7 cannot be extended to
arbitrary squarefree sextics. Here is an exact counterexample, using the
curve in [5, Section 4.4, equation (10)] after completing the y-square:

\[
D(x)=\frac{x(2x-1)(x^2-3x+8)(2x^2+x+3)}4.
\tag{28}
\]

Define

\[
\begin{aligned}
A(x)&=8x^4-36x^3+118x^2-159x+162,\\
B(x)&=8x^4-20x^3+54x^2-27x+54,\\
F_+(x)&=x(2x^2+x+3)A(x)^2,\\
F_-(x)&=(2x-1)(x^2-3x+8)B(x)^2.
\end{aligned}
\]

Exact multiplication gives

\[
F_+-F_-=23328.
\tag{29}
\]

Consequently

\[
S=\frac{F_++F_-}{256},\qquad T=\frac{AB}{64},\qquad
S^2-DT^2=\left(\frac{729}{8}\right)^2.
\tag{30}
\]

D is monic and squarefree, and deg S=11, deg T=8. The two quadratics in
(28) both have discriminant -23; they are irreducible over Q and have
no real roots. Thus D is not completely split over Q, and (30) supplies
no rational PTE pair. The monic polynomials F₊/128 and F₋/128 differ
by a constant, but neither splits over Q. This shows concretely why
Pell solvability alone does not settle the remaining problem.

`scripts/check_sextic_pell.py` verifies (29)–(30), the degrees and
squarefreeness by exact rational/integer coefficient arithmetic without
third-party dependencies. The source curve and the Pell/torsion method
are not claimed as new.


### 4.8. Three double roots: exceptional-locus exclusion and exact normal form

**Theorem 12.** Suppose a symmetric rational ideal size-eleven pair has
eight distinct values per side, and translate so B=-A. Its three double
entries on side A are u,v,w. Put

\[
p=u+v+w,\quad q=uv+uw+vw,\quad r=uvw,\quad
H=p^3r+p^2q^2-3pqr+3r^2.
\tag{31}
\]

Then pH≠0. In particular, the three double entries cannot sum to zero.
More generally the entire exceptional locus H=0 is excluded.

The root polynomial is uniquely determined by u,v,w:

\[
P(x)=(x^3-px^2+qx-r)^2Q(x),\qquad
Q=x^5+ax^4+bx^3+cx^2+dx+e,
\tag{32}
\]

where

\[
\begin{aligned}
a&=2p,\\
b&=p^2+2q-\frac{r^2(p^3+r)}{pH},\\
c&=2pb-2(p^3+pq-r),\\
d&=q^2+2pr-\frac{r^2(p^5+2p^3q+qr)}{pH},\\
e&=\frac{2r(pq-r)^3}{pH}.
\end{aligned}
\tag{33}
\]

The PTE constant is

\[
C=2P(0)=\frac{4r^3(pq-r)^3}{pH}.
\tag{34}
\]

Conversely, choose distinct nonzero rational u,v,w with no opposite pair
and with pH≠0. Form Q by (33). If Q has five distinct rational roots,
none equal to u,v,w, then the roots of P in (32), with their stated
multiplicities, and their negatives give a symmetric rational ideal
size-eleven pair with eight distinct values per side. Clearing
denominators gives an integer pair. Complete rational splitting of Q
is the essential unresolved condition, not a consequence of (33).

**Proof.** The multiplicity bound in Theorem 5 gives exactly three
double entries and five simple entries on each side. Disjointness of
A and -A gives

\[
r\ne0,\qquad pq-r=(u+v)(u+w)(v+w)\ne0.
\]

The coefficients of x¹⁰,x⁸,x⁶,x⁴,x² in P must vanish. Expansion of
(32) gives the linear system M(a,b,c,d,e)^t=y, with

\[
M=\begin{pmatrix}
1&0&0&0&0\\
p^2+2q&-2p&1&0&0\\
2pr+q^2&-2pq-2r&p^2+2q&-2p&1\\
r^2&-2qr&2pr+q^2&-2pq-2r&p^2+2q\\
0&0&r^2&-2qr&2pr+q^2
\end{pmatrix},\qquad
y=\begin{pmatrix}2p\\2pq+2r\\2qr\\0\\0\end{pmatrix}.
\tag{35}
\]

To handle singular parameters without dividing by a determinant, use the
polynomial row vector

\[
L=\bigl(qr(q^3+q^2p^2-3qrp+3r^2),\ -r(q^3+r^2),\
r(q^2+rp),\ -r(q+p^2),\ r+p^3\bigr).
\]

Direct polynomial multiplication yields

\[
LM=(0,0,0,0,pH),\qquad Ly=2r(pq-r)^3.
\tag{36}
\]

Thus every solution of the even-coefficient equations satisfies
pHe=2r(pq-r)³≠0. Both p=0 and H=0 are impossible. Moreover

\[
\det M=4p(pq-r)H\ne0.
\]

The system has a unique solution; substitution verifies exactly the
coefficients (33). The constant coefficient is r²e, giving (34).

Conversely, (33) solves all five even-coefficient equations, so
P(x)+P(-x)=2r²e=C≠0. The stated splitting conditions give eleven rational
entries with exactly the required multiplicities. A root common to P(x)
and -P(-x) would force C=0, so the two supports are disjoint. The
constant-difference identity gives equal power sums through order ten. ∎

This is a complete rational splitting formulation of the eight-point
symmetric boundary, not a construction. Scaling u,v,w scales the entire
pair; since u≠0 one may set u=1, leaving two rational parameters. No
claim that the resulting quintic has rational roots is made.

### 4.9. Real splitting does not eliminate the eight-point boundary

The hypotheses on rational splitting in Theorem 12 cannot be replaced
by a real-root argument. Take (u,v,w)=(1,5,9). Formula (33) gives

\[
Q(x)=\frac{913x^5+27390x^4+312646x^3+1682790x^2
+4236661x+3951360}{913}.
\tag{37}
\]

Let N=913Q. Its values at -10,-9,-8,-6,-5,-4,-3,-2 are, respectively,

\[
-182250,\ 2520,\ -45864,\ -21750,\ 22680,\ -3060,\ -58224,\ 117054.
\]

Hence Q has a real root in each of the five disjoint intervals
(-10,-9), (-9,-8), (-6,-5), (-5,-4), and (-3,-2). As its degree is
five, these are five distinct real roots. They are disjoint from the
positive double roots 1,5,9. Thus the construction gives a genuine
symmetric real PTE pair of the prescribed multiplicity pattern.

It does not give a rational pair. Modulo 5,

\[
N(x)\equiv3x(x-1)(x+1)(x^2-2),
\tag{38}
\]

and x²-2 is irreducible over F₅. Since 913 is a unit modulo 5, complete
rational splitting of N would reduce to complete splitting modulo 5:
by the rational-root theorem all root denominators divide 913. This
contradicts (38). The example therefore separates the genuine rational
splitting obstruction from the already satisfiable real-root conditions.


### 4.10. An exact completely split Pell criterion

The genus-two torsion condition is only necessary. Requiring one more
polynomial to split gives an exact equivalence with the smallest
remaining general support case.

**Theorem 13.** A nontrivial rational ideal size-eleven pair with total
support size fourteen exists if and only if there are polynomials
D,T,S∈Q[x] and k∈Q\{0} such that

\[
S^2-DT^2=k^2,\qquad
\deg D=6,\quad\deg T=8,\quad\deg S=11,
\tag{39}
\]

where D is monic, squarefree and completely split over Q, and T also
splits completely over Q. No prior squarefreeness or coprimality
assumption on T is necessary; those properties follow from (39) and
complete rational splitting.

**Proof.** Suppose first that the PTE pair exists. Its monic root
polynomials P,Q differ by C≠0. By Theorem 5, the two disjoint supports
have eight double roots and six simple roots in total: if these numbers
are d and s, then d+s=14 and 2d+s=22. Let U be the monic product of the
double roots and D the monic product of the simple roots. Both split
over Q, D is squarefree, and PQ=DU². Set S=P+Q, T=2U, and k=C.
Then S²-DT²=(P-Q)²=C², with the required degrees.

Conversely, let s be the leading coefficient of S and put

\[
P=\frac{S+k}{s},\qquad Q=\frac{S-k}{s}.
\tag{40}
\]

These are monic degree-eleven rational polynomials with difference
2k/s≠0, and

\[
PQ=D\left(\frac{T}{s}\right)^2.
\]

The right side splits completely over Q. Hence P and Q separately
split over Q, and their roots give a nontrivial ideal size-eleven pair
by the constant-difference criterion. They have disjoint supports.

Theorem 5 now forces T to be squarefree and coprime to D. Indeed, a
repeated root of T would have multiplicity at least four in PQ, while
a common root of D and T would have multiplicity at least three.
Disjointness would put the entire multiplicity in P or Q, contradicting
Theorem 5. There are therefore exactly six simple roots from D and
eight double roots from T across the two supports, giving total support
fourteen. ∎

Clearing denominators in (40) gives an integer PTE pair after scaling
its rational roots. The theorem does not construct D,T,S, nor does it
identify complete splitting with the existence of a rational torsion
point alone. The exact task left is simultaneous complete rational
splitting of D and T in (39). In the symmetric eight-value case,
Theorem 12 gives the separate exact quintic criterion, including its
necessary distinctness and disjointness conditions.


## 5. A distinct-entry construction route: shift chains

These results concern lifting a symmetric ideal size-ten pair. Theorem 14
in Section 5.32 proves normalized finiteness for all eighteen chain types,
and Section 5.31 excludes the previously named sixth mode. No size-eleven
witness or general nonexistence theorem is established. The lifting
mechanism and Faltings's theorem are standard; the explicit all-type
reductions and their exceptional-branch analysis are given here.

### 5.1 Exact lifting criterion and eighteen chain types

Let A,B be disjoint rational sets of size ten with matching powers 1
through 9. Their monic polynomials satisfy P-Q=C != 0. For a nonzero
rational shift t put F=P(x)Q(x-t), G=Q(x)P(x-t). Then

 F-G=C(Q(x-t)-Q(x))

has degree nine and leading coefficient -10tC. The monic gcd H of F,G
has degree

 c=|A intersect (A+t)|+|B intersect (B+t)|.

Indeed, every common root lies in exactly one of these two intersections
and has multiplicity one in both products. Cross-side overlaps cannot
be common roots because A and B are disjoint. Thus c<=9. After cancelling
H, eleven entries remain on each side exactly when c=9; the resulting
monic polynomials then differ by the nonzero constant -10tC. Their roots
are distinct exactly when A intersect (B+t) and B intersect (A+t) are
both empty. This is a necessary and sufficient criterion for this lift.

Assume now A=-A and B=-B, and normalize t=1 by scaling. Join input
vertices whose difference is one. For a valid lift there are twenty
vertices, nine edges, and no edge between colors A and B. Each component
is a monochromatic consecutive arithmetic chain; the graph is a union
of paths, so it has eleven components. Reflection permutes them.
A reflection-fixed chain has center zero. It cannot have odd length,
since that would put zero in a symmetric set of even cardinality.
Two even centered chains cannot coexist: both contain -1/2 and 1/2.
The odd number of components therefore forces exactly one centered
chain of length 2m, 1<=m<=5. Exchange colors so this chain belongs to A.
The remaining components form five reflected pairs.

If aa and bb are the multisets of lengths of the reflected chain pairs
in A and B, respectively, their exact constraints are

 len(aa)+len(bb)=5, sum(aa)=5-m, sum(bb)=5.

The positive integer partitions give eighteen types. This is a finite
type classification, with no height bound on the rational centers.
For a length L define O_L={j-(L-1)/2:0<=j<L}. A reflected pair of
center +/-z consists of z+O_L and -z+O_L. The centered chain is O_(2m).
Matching the four even moments 2,4,6,8 is sufficient for matching all
moments through 9. Conversely, rational centers solving these equations
produce a valid lift if all twenty vertices are distinct and no
cross-color unit difference occurs: the prescribed chains already give
nine edges, and the degree-nine gcd bound forbids additional edges.
Scaling clears denominators of a successful rational output.

### 5.2 Four modes reduce to conics and cubic splitting

Consider aa=(5-m), bb=(1,1,1,2), for m=1,2,3,4. Let +/-u be the
centers of the noncentral A chain pair, +/-v the centers of the length-two
B pair, and +/-r_i the three singleton B pairs. Set a=u^2, b=v^2.
For k=1,2,3,4 define

 p_k = (1/2) sum_{d in O_(2m)} d^(2k)
       + sum_{d in O_(5-m)} (u+d)^(2k)
       - sum_{d in O_2} (v+d)^(2k).

Symmetry of each offset set makes these polynomials in a,b. The required
relations are sum_i (r_i^2)^k=p_k. Set

 e1=p1, e2=(p1^2-p2)/2, e3=(p3-e1*p2+e2*p1)/3,
 K(X)=X^3-e1 X^2+e2 X-e3.

The first three moment conditions say the r_i^2 are the roots of K.
Newton's identities give the fourth condition F=p4-e1*p3+e2*p2-e3*p1=0.
Direct polynomial expansion yields, in order m=1,2,3,4,

 F1=-4(a-b)(a-5b-4)(a^2-2ab-2a+b^2-2b+1),
 F2=(b-1)(4a-10b-9)(16a^2-32ab-8a+16b^2-8b+1)/8,
 F3=-4(a-b)(b-4)(b-1)(3a-5b-3),
 F4=4(b-9)(b-4)(b-1)(4a-5b-1).

All factors except the indicated linear relation in a force input
collisions. The a=b factor makes the length-two B chains overlap the
length-four or length-two A chains. The quadratic factors in F1 and F2
are, respectively,

 [u^2-(v+1)^2][u^2-(v-1)^2],
 16[u^2-(v+1/2)^2][u^2-(v-1/2)^2];

they also force A/B overlaps. Finally b=j^2, 1<=j<m, puts a length-two
B chain into the centered A chain. Hence any valid input in these modes
satisfies, respectively,

 a=5b+4, a=(10b+9)/4, a=(5b+3)/3, a=(5b+1)/4.

These conics have rational parametrizations: write u^2=alpha*v^2+s^2,
with (alpha,s)=(5,2),(5/2,3/2),(5/3,1),(5/4,1/2), and put

 v=2st/(alpha-t^2), u=s+tv.

The exceptional points with v=0 are invalid because the reflected B chains
coincide. On each conic, completely splitting K into three distinct
rational-square roots, together with the noncollision conditions in
5.1, suffices for an eleven-entry lift. Square discriminant alone does not.

### 5.3 Explicit genus-two conditions and normalized finiteness

After substituting the four respective conic relations in K, its
discriminants are

 D1=-4h1, D2=-b*h2/4, D3=-4b*h3/243, D4=-h4/64,

where

 h1=800b^6-2625b^5+12050b^4-529b^3-54448b^2-53744b-14400,
 h2=675b^5-7950b^4-26125b^3-35704b^2-23248b-6144,
 h3=900b^5-7100b^4-21675b^3+71658b^2+3321b-93312,
 h4=50b^6-1875b^5+20600b^4-50624b^3+102912b^2-2392064b-3686400.

A valid input has b=v^2 with v!=0. If K splits over Q, its discriminant
is d^2 for rational d. It therefore gives a rational point on, respectively,

 C1: y^2=-h1(b), with y=d/2;
 C2: y^2=-h2(b), with y=2d/v;
 C3: y^2=-3h3(b), with y=27d/(2v);
 C4: y^2=-h4(b), with y=8d.

The right sides are squarefree polynomials of degrees 6,5,5,6. Thus the
smooth projective completions all have genus two. The exact certificate
`scripts/check_shift_chains.py` checks the discriminant identities and
constant gcds with the derivatives, as well as the four factorizations
using two expansions of the chain moments.

**Finiteness theorem.** Each of these four chain modes has only finitely
many valid rational configurations after the shift is fixed to one and
the symmetry center to zero.

**Proof.** Faltings's theorem gives finitely many rational points on each
C_m. Each b gives finitely many u,v with the required squares. K then
fixes the three squared singleton centers, leaving finitely many sign
and permutation choices. Therefore there are finitely many valid seeds
and outputs in these normalized modes. Arbitrary scaling and translation
are excluded from the normalization. This proves no effective height
bound and no nonexistence assertion.

The curves have the rational control points (b,y)=(1,336),(4,1920),
(9,3240),(16,384). They produce only the trivial A=B consisting of the
ten half-integers +/-1/2,+/-3/2,+/-5/2,+/-7/2,+/-9/2; the certificate
explicitly verifies this rejection. Their rational points have not been
completely determined. Any further point must still pass the rational-
square, complete-splitting, and noncollision conditions. These four
reductions alone do not treat the other fourteen chain types; the complete
normalized finiteness theorem for all eighteen appears in Section 5.32.
Constructions outside this shift-lift route remain outside its scope.
The finiteness input is the established theorem of G. Faltings,
*Endlichkeitssaetze fuer abelsche Varietaeten ueber Zahlkoerpern*,
Inventiones Mathematicae 73 (1983), 349-366,
https://doi.org/10.1007/BF01388432.

### 5.4 A fifth mode: a smooth plane quartic

This gives a fifth normalized finiteness result, not a distinct size-eleven construction or an exclusion of the mode.

Consider m=1, aa=(1,3), bb=(1,1,3). The centered A chain is {-1/2,1/2}. Let the length-three A and B chain centers be +/-u and +/-v, and the A singleton centers +/-w. Put a=u^2,b=v^2,c=w^2. The two B singleton-center squares must have power sums
 p_k=4^(-k)+sum_{d=-1}^1 (u+d)^(2k)+c^k-sum_{d=-1}^1(v+d)^(2k), k=1,2,3,4.
Set e1=p1,e2=(p1^2-p2)/2. Their quadratic is X^2-e1 X+e2. The remaining two Newton conditions are
 F=p3-e1*p2+e2*p1=0, G=p4-e1*p3+e2*p2=0.
Exact expansion shows
 F=3(a-b)[4a^2-32ab+12ac-29a+40b^2-24bc+106b-21c+34]/4.
The factor a=b forces overlap of the two length-three pairs and is invalid.
Write D=4a-8b-7. If D=0, F/(a-b)=-3(4b-9)(4b-1)/8, so b=1/4 or9/4. Either v=+/-1/2 or+/-3/2 puts a B length-three point at +/-1/2, an A vertex. Thus D!=0 for any valid input, and
 c=-(4a^2-32ab-29a+40b^2+106b+34)/(3D).
Substitution in the other Newton condition gives
 G/(a-b)=2[(a-b-1)^2-4b]H/(3D),
 H=8a^2-40ab-58a+80b^2+50b+95.
The bracketed factor is [u^2-(v+1)^2][u^2-(v-1)^2], again forcing overlap of the length-three chains. Therefore H=0 is necessary.

Consequently (u,v) gives a rational point on the projective plane quartic
 T(U,V,W)=8U^4-40U^2V^2-58U^2W^2+80V^4+50V^2W^2+95W^4=0,
with W!=0, u=U/W,v=V/W.
This quartic is smooth even over the algebraic closure. Write T as a quadratic in U^2,V^2,W^2 with symmetric derivative matrix
 M=[[16,-40,-58],[-40,160,50],[-58,50,190]].
At a singular projective point, the nonzero coordinates' squares would lie in the kernel of the corresponding principal submatrix. All nonempty principal minors are nonzero: diagonal entries 16,160,190; two-by-two minors 960,-324,27900; determinant -163840. No singular point exists. A smooth plane quartic is geometrically integral and has genus 3. Faltings therefore gives finitely many rational points, hence finitely many normalized valid configurations in this fifth mode. Indeed a,b determine c, and the quadratic fixes the other two squared centers, leaving finitely many sign choices.

A rational point alone is not sufficient: a,b,c must be rational squares; the quadratic must split into rational squares; and all twenty input vertices must be distinct with no cross-color unit edge. The visible points (u,v)=(3/2,1/2),(5/2,3/2) belong to the excluded D=0 branch. No complete rational-point determination is claimed.

For effective arithmetic, one elliptic quotient is
 y^2=-240v^4+760v^2+81, y=8u^2-20v^2-29.
This may assist a future descent or elliptic Chabauty computation, but neither ranks nor a complete point set have been certified here. The necessary quartic has genuine rational points, so a simple local-solubility obstruction cannot exclude the whole projective curve.

The exact SymPy certificate `scripts/check_fifth_mode.py` checks the
elimination identities, exceptional branch, principal minors, quotient
identity, and excluded control points.

### 5.5 A sixth mode and a nonexceptional rational control

For m=2, aa=(1,2), bb=(1,1,3), let the centered A chain be
{-3/2,-1/2,1/2,3/2}. Write +/-u for the A length-two centers,
+/-w for its singleton centers, and +/-v for the B length-three centers.
Put a=u^2,b=v^2,c=w^2. The two B singleton-square roots must have moments

 p_k=(1/2)sum_{d in O_4}d^(2k)
      +sum_{d in O_2}(u+d)^(2k)+c^k
      -sum_{d in O_3}(v+d)^(2k), k=1,2,3,4.

As in 5.4 set e1=p1,e2=(p1^2-p2)/2,
F=p3-e1*p2+e2*p1 and G=p4-e1*p3+e2*p2. Define

 D=8a^2-48ab+4a+48b^2+24b-9,
 N=48a^2b-8a^2-192ab^2-168ab+50a+160b^3+352b^2+10b-27.

Exact expansion gives F=-3(N-2Dc)/16. The resultant in a is

 Res_a(N,2D)=-128(4b-25)(4b-9)^2(4b-1)^3.

Thus D=0 and F=0 force b=1/4,9/4,25/4. In each case the B
length-three pair has a point in the central A chain. Valid inputs have
D!=0 and c=N/(2D). Substitution in G yields

 G=(4b-9)(4b-1)JH/(64D),
 J=16a^2-32ab-8a+16b^2-8b+1,
 H=24a^2-80ab-100a+80b^2+40b+81.

The factors 4b-9 and 4b-1 force the same central-chain collisions.
Also J=16[u^2-(v+1/2)^2][u^2-(v-1/2)^2], so J=0 makes an
A length-two chain overlap a B length-three chain. Hence H=0 is necessary.
This gives the plane quartic

 24U^4-80U^2V^2-100U^2W^2+80V^4+40V^2W^2+81W^4=0,

with u=U/W,v=V/W and W!=0. Its symmetric derivative matrix in
U^2,V^2,W^2 is

 [[48,-80,-100],[-80,160,40],[-100,40,162]].

The nonempty principal minors are 48,160,162,1280,-2224,24320,-829440.
The support/kernel argument of 5.4 proves geometric smoothness, hence
genus three. Faltings and the finite choices for c, square roots, and
the two B singleton roots give normalized finiteness in this sixth mode.
This is not an effective point classification or nonexistence proof.

There is a nonexceptional rational point

 (u,v)=(64/61,1/122), D=70543200/13845841,
 c=N/(2D)=26785/14884=5*11*487/(122^2).

It lies on the quartic and has D!=0, but c is not a rational square:
its 5-adic valuation is one. Consequently it cannot give the required
rational singleton center w and is not a PTE witness. This control
shows that neither rationality of the quartic point nor exclusion of
the exceptional denominator alone solves the remaining splitting problem.

The exact certificate `scripts/check_sixth_mode.py` verifies the moment
identities, exceptional resultant, collision factor, smoothness minors,
and this nonexceptional rejection. No finite search failure is used.
The remaining twelve types are treated for normalized finiteness in
Section 5.32. No eleven-distinct-integer solution or general impossibility
theorem is established.

### 5.6 Effective rank computation for the sixth-mode Jacobian

The sixth-mode curve C has Jacobian rank exactly seven over Q. This
certifies a limitation of ordinary Chabauty-Coleman here, not a PTE
nonexistence theorem or a determination of C(Q).

Write its affine equation as

 T(u,v)=24u^4-80u^2v^2-100u^2+80v^4+40v^2+81=0.

The three sign involutions give these genus-one quotients:

 E_v: y^2=-80v^4+760v^2+139, y=12u^2-20v^2-25;
 E_u: z^2=-5(u^2-1)(u^2-19), z=10v^2-5u^2+5/2;
 E_r: h^2=139r^4+1120r^2-1520,
      r=u/v, h=(81/v^2-50r^2+20)/2.

They are degree-two quotients by u->-u, v->-v, and simultaneous
negation, respectively. The quartic polynomials are squarefree. The
point (u,v)=(1,1/2) gives rational points on all three smooth projective
models, so they are elliptic curves after choosing origins. The rational
maps extend to morphisms of the smooth projective curves.

The Jacobian of C is Q-isogenous to E_v x E_u x E_r. Indeed, the regular
differential basis du/T_v, u du/T_v, v du/T_v has characters (-1,-1),
(+1,-1), (-1,+1) under the two independent sign involutions. Pullbacks
from the three elliptic quotients span these three distinct lines.
The induced homomorphism from J(C) to the product therefore has
invertible differential and is an isogeny.

For explicit rank computations, a quartic q^2=A t^4+B t^2+C maps to
Y^2=X(X^2+B X+AC) by X=A t^2,Y=A t q. This is a nonconstant
morphism between genus-one curves, and hence gives an isogeny after
choosing origins and translating the target origin. A standard
2-isogeny then gives Y^2=X(X^2-2B X+B^2-4AC). Scaling X=s^2 x,
Y=s^3 y with s=4,1,2 for E_v,E_u,E_r, respectively, gives models

 L_v: y^2=x(x^2-95x+2430),
 L_u: y^2=x(x^2-200x+8100),
 L_r: y^2=x(x^2-560x+131220).

These isogenies preserve Mordell-Weil rank.

Here is a complete elementary 2-isogeny descent certificate, rather
than a rank estimate from searching for points. For
L:y^2=x(x^2+A x+B), its 2-isogenous partner is
L':y^2=x(x^2-2A x+A^2-4B). The classical descent homomorphism alpha
sends a point with x!=0 to its rational squareclass x, infinity to 1,
and (0,0) to the squareclass of B. Its image is L(Q) modulo the image
of the dual isogeny. If alpha and alpha' have image sizes m,m', then

 2^rank(L(Q))=m*m'/4.

For completeness the factor four follows from composing the two
isogenies to [2]: L(Q)/2L(Q) has order 2^rank*|L[2](Q)|, while the
kernel correction in the two quotient indices is 4/|L[2](Q)|.
The squareclass descent and its quartic criterion are the classical
2-isogeny descent (see J. H. Silverman and J. Tate, *Rational Points on
Elliptic Curves*, Springer).

Every image class has a signed squarefree representative d dividing B.
Indeed, at any prime not dividing B, the equation forces the valuation
of x to be even; negative valuations are also even. For x=d(U/V)^2
with coprime integers U,V, the corresponding necessary homogeneous
quartic equation is

 Z^2=d U^4+A U^2 V^2+(B/d)V^4.

A rational Z here is automatically an integer. If d and B/d are both
negative and A<=0, there is no real point. Otherwise check this equation
modulo 256. Since U,V are coprime, either V is odd or U is odd and V is
even. Scaling the odd coordinate to one leaves exactly the two patches
(U,V)=(t,1) for all t modulo 256, and (1,t) for even t modulo 256;
Z is scaled by the square of the same unit. Exhaustively checking these
patches gives a rigorous upper bound on the descent image.

The resulting bounds are attained, because the certificate supplies
an exact rational point for every surviving class. The six exact images,
with each partner listed immediately after its curve, are:

 (-95,2430):       {1,5,6,30};
 (190,-695):       {1,-5,139,-695};
 (-200,8100):      {1,5,6,30};
 (400,7600):       {1,-5,19,-95};
 (-560,131220):    {1,5,6,30};
 (1120,-211280):   {1,-5,19,139,-95,-695,2641,-13205}.

For example, L_v has points (24,132) and (45,90), which supply classes
6 and 5, while infinity and (0,0) supply 1 and 30. Its partner has
(-125,1050) and (139,2502), supplying -5 and 139. All witnesses for the
other models, together with exact curve-equation checks, are listed in
`scripts/check_sixth_exact_ranks.py`. No point-search completeness is
assumed: the modular exclusions give the upper bounds and the explicit
points attain them. The image-size products are 16,16,32, so the ranks
are exactly 2,2,3. Consequently

 rank J(C)(Q)=2+2+3=7.

The dependency-free certificate verifies the complete squareclass
lists, real/modulo-256 exclusions, every witness, the partner models,
and the normalization from the quartic coefficients. The companion
`scripts/check_sixth_rank_maps.py` checks the quotient identities and
squarefreeness using SymPy.

Since genus(C)=3, the ordinary Chabauty-Coleman rank<genus hypothesis
fails. This does not rule out covering methods, other refinements of
Chabauty, or another proof. No complete C(Q), effective height bound,
valid size-eleven pair, or general impossibility result follows.

### 5.7 Incorporating the missing center-square condition

Rank seven in 5.6 concerns the necessary quartic, not all PTE lifting
conditions. A more selective route is to impose the recovered A singleton
center at the outset. It has an explicit rank-one elliptic projection.
This does not yet classify its simultaneous rational lifts.

Parametrize the conic H(a,b)=0 of 5.5 by

 d=3t^2-10t+10, p=6t^2-5t+5,
 q=3t^2+26t+10, r=3t^2-10t+100,
 a=2p/d, b=q/(4d), c=r/(4d).

The parameter comes from a=1+t(b-1/4), and t ranges over P^1(Q).
The quadratic d has no rational root. The identity recovering c is

 4N-Dold(-6a+20b+21)=(2a-4b+1)H,

where Dold and N denote D,N in 5.5. Thus, on H=0 away from Dold=0,
c=(-6a+20b+21)/8. Under this parametrization

 Dold=288t^2(t-4)(3t-4)/d^2.

The exceptional rational parameters are t=0,4,4/3, all corresponding
to the already excluded collisions. The parameter at infinity gives
(a,b,c)=(4,1/4,1/4); it must be retained during parametrization, but its
input has a cross-color collision and is not a PTE witness.

The three center conditions are the simultaneous rational squares

 u^2=2p/d, v^2=q/(4d), w^2=r/(4d).

Their smooth projective normalization is a geometrically connected
multiquadratic cover of P^1 of degree eight and genus nine. Indeed the
four quadratics p,q,r,d are squarefree and pairwise coprime. Their eight
roots are distinct. The three squareclasses are independent in Qbar(t)^*/Qbar(t)^{*2}, where Qbar is an
algebraic closure of Q, since each numerator has roots absent from
the others. At the roots of p,q,r the inertia labels are respectively
(1,0,0),(0,1,0),(0,0,1), and at roots of d the label is (1,1,1).
Each branch point has inertia two; infinity is unramified. Riemann-Hurwitz
gives 2g-2=8(-2)+8*4=16, hence g=9. This is a necessary square-center
cover, not the full PTE condition: the two B singleton-square roots
must still be rational squares, and collisions must still be excluded.

The w-square projection has a particularly simple form. Put z=3t-5.
Then d=(z^2+5)/3 and r=(z^2+275)/3, and with

 s=2w(z^2+5)

its equation is

 s^2=(z^2+5)(z^2+275)=z^4+280z^2+1375.

This smooth genus-one quartic has a birational elliptic model

 E: Y^2=X(X^2-560X+72900),
 X=2s+2z^2+280, Y=2zX;
 z=Y/(2X), s=(X-2z^2-280)/2.

The displayed inverse applies where X!=0; the maps extend on the smooth
projective models, including the rational points at infinity. Thus this
is not merely a rank comparison with an unrelated elliptic curve.

The rank of E(Q) is exactly one. Its 2-isogenous partner has coefficients
(A,B)=(1120,22000). The same complete descent as in 5.6, now using the
real condition and residues modulo 256 and 3, gives exact images

 E: {1,5}; E': {1,-5,-11,55}.

Infinity and (0,0) give classes 1 and the squareclass of B. The point
(180,900) on E supplies 5. On E', (-1100,0) supplies -11 and
(-500,12000) supplies -5. Thus every surviving class is attained and
2^rank=2*4/4=2. These are exact upper and lower bounds, not a point-search
inference. For controls the quartic has points
(z,s)=(5/2,225/4) and (19/4,1449/16).

The exact certificate `scripts/check_sixth_center_lift.py` verifies the
parametrization, exceptional parameters, branch polynomial checks,
birational identities, witnesses and complete descent images.

This supplies an alternative to working directly with the rank-seven
Jacobian: an explicit simultaneous-lifting problem over a rank-one
elliptic curve. The reduction alone does not classify the simultaneous
lifts, supply their height bound, or exclude or construct a size-eleven
pair. Section 5.8 certifies the full elliptic generator. Ordinary
Chabauty on the rank-seven curve has not thereby been restored.

### 5.8 Full Mordell-Weil generator for the center-square projection

Let

 E: y^2=x(x^2-560x+72900),   P=(180,900),   T=(0,0).

Assume the already certified two-isogeny descent giving rank E(Q)=1.
Then

 E(Q) = Z P direct-sum (Z/2Z) T.

This assertion is a saturation theorem, not an inference from a bounded search.
The finite search below has a proved height cutoff covering every possible
proper generator of P.

#### Torsion and parity of a hypothetical index

The discriminant of this integral model is

 Delta=16*72900^2*22000=2^12*3^12*5^7*11.

Good reductions at 7 and 13 have respectively 8 and 12 points. The standard
injectivity of rational torsion under good reduction at primes greater than
2 implies that its order divides 4. Alternatively, prime-to-p injectivity
at these two primes gives the same conclusion. The only nonzero rational
2-torsion point is T: the quadratic factor has discriminant 22000, which is
not a rational square. A rational point R with 2R=T would, from duplication,
have x(R)^2=72900, hence x(R)=270 or -270. At these two values the right-hand
side of the curve equation is respectively 270^2*(-20) and 270^2*(-1100),
both negative. Thus E(Q)_tors={O,T}. In particular P has infinite order.

For this equation the standard 2-isogeny descent homomorphism is
alpha(O)=1, alpha(T)=72900 modulo squares, and alpha(R)=x(R) modulo squares
otherwise. Thus alpha(T)=1 and alpha(P)=5. If Q is any generator of the
free quotient and P=nQ+epsilon*T, then n is odd. If P were not a generator,
we could therefore choose |n|>=3.

#### An elementary global height bound

Use the scaled coordinate u=x/270. For u=a/b in lowest terms, including
u(O)=1/0, put H=max(|a|,|b|) and h(R)=log H. Its duplication map is

 F=27(a^2-b^2)^2,
 G=4ab(27a^2-56ab+27b^2),
 u(2R)=[F:G].

The following two bounds hold for every coprime pair (a,b):

 gcd(F,G) divides 11880=2^3*3^3*5*11,
 (121/48)*H^4 <= max(|F|,|G|) <= 440*H^4.

Here gcd means the positive gcd of |F|,|G|, also when one is zero.

For the gcd assertion, at 3 either both a,b are units, when G is a unit,
or exactly one is divisible by 3, when v_3(F)=3. At 2, unless both a,b
are odd F is odd; if they are both odd, the parenthesis in G is 2 modulo
4, so v_2(G)=3. For an odd prime p other than 3 dividing F and G, both
a,b are units and a=+/-b modulo p. The plus sign gives the parenthesis
in G congruent to -2b^2 and is impossible. The minus sign forces p=5 or
11. Writing that parenthesis as 27(a+b)^2-110ab then gives valuation
exactly one at either of these primes.

The upper real bound follows by the triangle inequality. For the lower
bound, exchange a,b if necessary (both homogeneous polynomials are
symmetric) and put v=b/a, with |v|<=1. If |v|<=5/6, then
F/H^4=27(1-v^2)^2 >=121/48. If -1<=v<=-5/6, then
|G|/H^4=4|v|(27+56|v|+27v^2)>=90. If 5/6<=v<=1, the polynomial
-27v^2+56v-27 is increasing and at least 11/12; hence
|G|/H^4 >=4*(5/6)*(11/12)=55/18>121/48.

Consequently

 4h(R)-log(51840/11) <= h(2R) <=4h(R)+log 440.

Define the doubled canonical height

 h_c(R)=lim_{k->infinity} 4^(-k) h(2^k R).

Changing x to x/270 changes its naive logarithmic height by at most log270,
so this is exactly twice the usual canonical height. Telescoping the last
inequality proves

 h(R)-log(51840/11)/3 <= h_c(R) <= h(R)+log440/3.

The canonical-height theorem gives h_c(nR+torsion)=n^2*h_c(R).
Since u(P)=2/3, h(P)=log3, and therefore

 h_c(P) <= log3+log440/3.

If P=nQ+epsilon*T with |n|>=3, the nonnegativity and quadraticity of
canonical height imply

 h(Q) <= (log3+log440/3)/9 + log(51840/11)/3 < log24.

The last strict inequality is certified entirely with integers by

 27*440*51840^9 < 24^27*11^9.

#### Exhaustive classification below the proved cutoff

The exact script checks all reduced a/b with |a|<=24, 1<=b<=24 and
x=270a/b. The only x-coordinates for which the curve's right-hand side
is a rational square are

 u=0       : T;
 u=2/3     : +/-P;
 u=3/2     : +/-(P+T);
 u=15/8    : +/-2P;
 u=8/15    : +/-(2P+T).

The point O is handled separately. The corresponding nonnegative y values
are 0, 900, 2025, 38475/8, 1368. Thus every possible Q under the hypothetical
proper-index bound is torsion or has free class +/-[P] or +/-2[P]. None can
satisfy [P]=n[Q] with |n|>=3 in the torsion-free quotient. The hypothetical
proper index is impossible, proving the displayed group decomposition.

#### Scope and standard inputs

The standard inputs beyond elementary arithmetic are the rank-one descent
already proved in Section 5.7, the 2-isogeny descent homomorphism,
good-reduction torsion injectivity, and canonical-height quadraticity and
torsion invariance. References for the latter two are Silverman,
The Arithmetic of Elliptic Curves, Chapters VII section 3 and VIII section 9.
No numerical height approximation, conjectural analytic rank, unverified
software saturation routine, or unproved search cutoff is used here.

This classifies rational points of the elliptic projection as nP+epsilon*T.
It does not classify which multiples meet the remaining simultaneous square
conditions, or resolve the original distinct size-eleven PTE problem.

The exact certificate is `scripts/check_rank_one_generator.py`. Combined
with the birational maps in 5.7, this reduces the remaining necessary
simultaneous-square lifting to points nP+epsilon*T with n an integer and
epsilon in {0,1}. The points at infinity and excluded collision cases
remain subject to the handling in 5.7. No bounded test of multiples is
used to infer completeness of those lifting conditions.

### 5.9 Uniform index restrictions and a scaled elliptic height bound

Write an elliptic candidate uniquely as R=nP+epsilon*T using 5.8,
where n is an integer and epsilon is 0 or 1. A necessary condition for
simultaneous rational center squares is

 epsilon=0: n is 0,1,120 modulo 240;
 epsilon=1: n is 0,120,239 modulo 240.

This excludes entire infinite residue classes, not just tested multiples.
It is not a complete lifting criterion.

**Proof of the congruences.** Use the good reduction primes 13,37,59,113.
The respective orders of P modulo those primes are 6,40,16,120. For a
finite reduced point (X,Y) with X!=0, form

 z=Y/(2X), t=(z+5)/3, d=3t^2-10t+10.

If d!=0, require that 2(6t^2-5t+5)/d and
(3t^2+26t+10)/(4d) are squares in the residue field, allowing zero.
If the point is infinity, X=0, or d=0, retain it conservatively; do not
exclude a rational candidate because an affine formula is undefined.
For every retained affine computation, X and d are units and reduction
of the rational center squares is valid. Good reduction is a group
homomorphism, so these conditions depend only on n modulo the displayed
order and epsilon. The exact finite-field certificate computes every
point in those four cyclic groups and intersects their conditions modulo
the least common multiple 240, obtaining exactly the two displayed sets.
This finite calculation proves a uniform necessary condition for all
integers n, with no unproved search cutoff.

Among |n|<120 the surviving points are O,T,P,-P+T. The first two map
to t=infinity, with (a,b,c)=(4,1/4,1/4); the B length-three chain
meets the centered A chain. The latter two map to t=5/2, with
(a,b,c)=(16,25/4,25/4); the B point v-1=3/2 again belongs to the
central A chain. Thus every valid sixth-mode PTE lift must have |n|>=120.

This also yields an explicit lower bound on the **scaled elliptic
coordinate height**, not directly on the original integer PTE entries.
With h(R)=log H(x(R)/270) and h_c as in 5.8,

 h_c(P)=h_c(4P)/16 >= [log77763-log(51840/11)/3]/16,

because the exact scaled-coordinate height of 4P is 77763. For every
remaining valid candidate, torsion invariance and quadraticity give

 h(R) >= h_c(R)-log440/3
      >= 900log77763-300log(51840/11)-log440/3.

Consequently H(x(R)/270)>10^3298. The last strict inequality is
verified without floating-point logarithms by

 77763^2700 * 11^900 > 10^9894 * 51840^900 * 440.

The dependency-free certificate `scripts/check_lift_congruences.py`
checks the orders, complete residue tests, collision controls, exact 4P
height, and integer inequality. This is a rigorous bounded-index exclusion
and a necessary height bound for one construction mode, not an upper
bound or a proof that the remaining progressions contain no valid lift.
The simultaneous-square classification and the original size-eleven
problem remain unresolved.

### 5.10 Two squareclass covers and explicit genus-two quotients

This is a finite covering reduction, not a determination of the rational
points on the covers. It supplies explicit genus-two curves on which a
certified rank-and-Chabauty calculation could be attempted. No favorable
Jacobian rank or complete rational-point list is assumed.

Put z=3t-5 and write z=A/B with coprime integers A,B, B>0. Define

 Q0=A^2+5B^2,
 Q1=4A^2+30AB+80B^2,
 Q2=A^2+36AB+185B^2,
 Q3=A^2+275B^2.

The center conditions become a=Q1/Q0, 4b=Q2/Q0, 4c=Q3/Q0.
Thus they hold exactly when there are one positive squarefree integer
Delta and integers C0,C1,C2,C3 satisfying

 Qi=Delta*Ci^2  (i=0,1,2,3).

Indeed, equality of the rational squareclasses of the four positive
integers is equivalent to the three square-ratio conditions. The
quadratics have no rational zero, so zero values need not be separated.

Only Delta=1 and Delta=5 are possible. If a prime divides Delta, it
divides all Qi, but it cannot divide B, since it divides Q0 and gcd(A,B)=1.
The difference Q3-Q0=270B^2 therefore restricts its prime factors to
2,3,5. If 2 divides Delta, then A,B are both odd; but Q0=6 modulo 8
has 2-adic valuation one and Q3=4 modulo 8 has valuation two, a
contradiction to equal squareclasses. The remaining mixed-parity
primitive cases have Q0=1 modulo 4. Its squarefree part must therefore
be 1 modulo 4, excluding 3 and 15 and leaving 1 or 5. In fact A must
be odd and B even: A even and B odd would make Q3 equal 3 or 7 modulo
8, impossible for Delta times an odd square with Delta=1 or 5.

This squareclass restriction is also implicit in the elliptic descent.
For a lift, write Q0=Delta*C0^2, Q3=Delta*C3^2 and choose signs such
that the quartic ordinate is s=Delta*C0*C3/B^2. The birational map gives

 X=2s+2z^2+280=Delta*(C0+C3)^2/B^2.

The factor C0+C3 is nonzero for finite z because Q3-Q0=270B^2.
Consequently Delta is exactly the two-isogeny class alpha(R). Since
alpha(P)=5 and alpha(T)=1, a point R=nP+epsilon*T has Delta=1 when
n is even and Delta=5 when n is odd. This is not an additional
restriction beyond the already certified elliptic group and descent.
The points at t=infinity have Delta=1 and are retained separately.

For i=0,1,2,3 put qi(z)=Qi(z,1). For each omitted index j and each
Delta in {1,5}, define the smooth projective curve

 G(j,Delta): Delta*y^2 = product of qi(z) over i != j.

Every simultaneous center lift with finite z maps to all four curves
with its own Delta, by

 y=Delta*product of Ci over i != j divided by B^3.

These are eight explicit genus-two curves. In the order j=0,1,2,3,
the sextic right-hand sides are

 4z^6+174z^5+3000z^4+56280z^3+537300z^2+2318250z+4070000;
 z^6+36z^5+465z^4+10080z^3+53175z^2+49500z+254375;
 4z^6+30z^5+1200z^4+8400z^3+27900z^2+41250z+110000;
 4z^6+174z^5+1920z^4+9300z^3+24300z^2+42150z+74000.

The four quadratic discriminants are -20,-380,556,-1100. Their six
pairwise resultants, in lexicographic order 01,02,03,12,13,23, are
8100,38880,72900,131220,1287900,364500. They are all nonzero, so each
sextic is squarefree, and a degree-six hyperelliptic equation has genus
two. At infinity each Delta=1 curve has two rational points. The known
finite center lift z=5/2 belongs to the Delta=5 curves; its positive
y-values for j=0,1,2,3 are 3375/2,3375/8,675/2,675/2.

One may select any one j for Delta=1 and any one j for Delta=5. A
complete rational-point classification on those two selected curves
would give a finite list of possible z, which must then be checked
against all center squares, the two B singleton-square conditions,
and the collision exclusions. The map is only a necessary condition:
a point on one genus-two quotient does not by itself give a center
lift, much less a PTE witness.

The full equations Qi=Delta*Ci^2 form two twists of a genus-seventeen
cover; the value of the reduction is their genus-two quotients, not
a claim that the full cover has smaller genus. The degree-sixteen
multiquadratic cover of the z-line has eight simple branch points,
so Riemann-Hurwitz gives 2g-2=-32+8*8=32. Its map to the genus-nine
center cover is an unramified degree-two map.

The certificate `scripts/check_sixth_genus_two_covers.py` verifies the
coordinate identities, the elementary local squareclass exclusions,
the sextic expansions, discriminants, resultants, the elliptic
squareclass identity, and the known rational controls. It supplies no
Jacobian rank bound, Mordell-Weil basis, Chabauty computation, rational
point completeness, or effective upper height bound. Those remain the
precise blockers to using this reduction as a complete classification.

### 5.11 The genus-nine center cover also has rank above genus

Let C9 be the genus-nine simultaneous-center cover of 5.7. Its Jacobian
has rank at least twelve over Q. This is an obstruction to ordinary
Chabauty-Coleman on C9, not a rational-point classification, an exclusion
of other methods, or a PTE nonexistence result.

For z=3t-5 use

 f0=z^2+5, f1=4z^2+30z+80,
 f2=z^2+36z+185, f3=z^2+275.

The seven nontrivial characters of the group (Z/2Z)^3 acting on C9
give its quotient curves. Six are the genus-one curves y^2=fi*fj
for 0<=i<j<=3; the seventh is the genus-three curve
H3:y^2=f0*f1*f2*f3. Pullback differentials from these quotients belong
to seven different character spaces, of dimensions 1,1,1,1,1,1,3.
Their direct sum is the nine-dimensional differential space of C9.
The resulting homomorphism of Jacobians therefore has invertible
differential and is an isogeny over Q:

 J(C9) ~ E_a x E_b x E_c x E_ab x E_ac x E_bc x J(H3).

The ranks of E_a,E_b,E_ab and E_c are already certified as 2,2,3,1.
The remaining two elliptic quotients have the models

 E_ac: y^2=x(x^2-2360x+1287900),
 E_bc: y^2=x(x^2-920x+364500).

For example, the binary quartics f1*f3 and f2*f3 have invariants
I=12ae-3bd+c^2 and J=72ace+9bcd-27ad^2-27eb^2-2c^3
satisfying I=A^2-3B, J=-2A^3+9AB for these two pairs (A,B).
The standard quartic Jacobian y^2=x^3-27Ix-27J is obtained from the
displayed cubic after translating by A/3 and scaling by 3. Both
quartics have rational points at infinity (leading coefficients 4 and 1),
so they are isomorphic to their Jacobians after choosing an origin.

A complete two-isogeny descent gives the following exact images, with
each partner immediately after its curve:

 (-2360,1287900): {1,5,159,795};
 (4720,418000): {1,-11,-95,1045};
 (-920,364500): {1,5};
 (1840,-611600): {1,-5,-11,139,55,-695,-1529,7645}.

For the upper bounds, enumerate the signed squarefree divisors d of B
and test the primitive quartic
Z^2=dU^4+A U^2 V^2+(B/d)V^4. The real obstruction and congruences
modulo 256, 3 and 125 leave exactly the displayed classes. Primitive
local solutions are checked on the patches (U,V)=(t,1) and (1,t), with
t divisible by the relevant prime in the second patch.

Infinity and (0,0) witness classes 1 and squareclass(B). The other
classes are attained by these exact rational points:

 E_ac: (405,14175), (3180,111300);
 E_ac': (-4400,66000), (-95,1425);
 E_bc: no additional point is needed;
 E_bc': (-1620,39600), (-1100,39600), (-695,31275),
        (556,20016), (880,39600), (30580,5504400).

The descent formula 2^rank=|alpha|*|alpha'|/4 gives rank 2 for each
new quotient. Consequently

 rank J(C9)(Q)=12+rank J(H3)(Q) >=12>9.

The dependency-free checker `scripts/check_center_cover_rank_obstruction.py`
verifies all finite local tests, class witnesses, partner coefficients,
and quartic invariants. It does not compute the rank of J(H3), rule out
quadratic/elliptic Chabauty on other covers, or solve any rational-point
classification. In particular no favorable rank assertion is made for
the genus-two targets in 5.10.

### 5.12 Genus-two Richelot targets: certified torsion and positive rank

This is a structural and arithmetic improvement of the genus-two reduction,
not a rank upper bound or rational-point classification. It singles out the
pair with omitted index j=2 for a more accessible descent, certifies that all
eight Jacobians have positive rank, and rules out a decomposition into elliptic
curves over Q. No favorable Chabauty rank condition is asserted.

The accompanying dependency-free certificate is
`scripts/check_genus_two_richelot_certificate.py`. Run it with Python 3. It verifies
polynomial identities, good reduction at every prime used, finite-field curve
and Jacobian counts, rational control points, and the integer-factor tests
below. Standard theorems invoked in the deductions are explicitly identified.

#### The exact Richelot isogeny for j=2

Write

- q0 = z^2+5;
- q1 = 4z^2+30z+80;
- q3 = z^2+275;
- C_delta: delta*y^2 = q0*q1*q3, for delta=1 or 5.

For polynomials A,B set [A,B]=A'B-AB'. Use the determinant with columns in
**ascending** coefficient order, namely constant, linear, quadratic. Then

    Delta = det [[5,0,1],[80,30,4],[275,0,1]] = -8100,
    L1 = [q1,q3] = -30*(z^2-68z-275),
    L2 = [q3,q0] = -540*z,
    L3 = [q0,q1] = 30*(z-1)*(z+5).

Consequently, on putting

    F(z)=z*(z-1)*(z+5)*(z^2-68z-275),

we obtain

    L1*L2*L3 = 486000*F,
    (L1*L2*L3)/Delta = -60*F.

The Richelot theorem supplies a Q-defined degree-four (2,2)-isogeny from
J(C_1) to the Jacobian of v^2=-60*F. Dividing v by 2 gives the model

    D_1: Y^2=-15*F(z).

For the second twist, first replace delta*y^2=f by the Q-isomorphic equation
w^2=delta*f using w=delta*y. Apply Richelot to (delta*q0,q1,q3): the determinant
is multiplied by delta, while the bracket product is multiplied by delta^2.
The new right-hand side is therefore -60*delta*F. For delta=5, dividing v by
10 gives

    D_5: Y^2=-3*F(z).

Thus J(C_delta) and J(D_delta) are Q-isogenous and have the same Mordell-Weil
rank. The negative signs in these two equations are essential.

For a reference fixing the convention, see Proposition 7 in
[Efficient computation of (2^n,2^n)-isogenies](https://link.springer.com/article/10.1007/s10623-024-01366-1),
which uses the ascending coefficient determinant and h_i=L_i/Delta.
Its equation Y^2=h1*h2*h3 is equivalent to v^2=(L1*L2*L3)/Delta by v=Delta*Y.

This is an isogeny of **Jacobians**, not a morphism C_delta to D_delta. A list
of D_delta(Q), even if complete, does not by itself classify C_delta(Q).
Rank bounds transfer, but original-curve point determination still requires
an appropriate Mordell-Weil subgroup and a justified point-completeness
method on C_delta, or a fully worked-out use of the correspondence.

#### Exact finite-field counts

For a good odd prime p, let N1=#C(F_p) and N2=#C(F_(p^2)), and set

    a1 = p+1-N1,
    a2 = (N2-p^2-1+a1^2)/2,
    P_p(T) = T^4-a1*T^3+a2*T^2-p*a1*T+p^2.

Then #J(F_p)=P_p(1). Direct enumeration gives the following values for both
C_delta and D_delta. Equality is also a useful independent check on the sign
and twist of the Richelot formula, though finitely many matching counts alone
would not prove an isogeny.

| delta | p | N1 | N2 | a1 | a2 | #J(F_p) |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 7 | 8 | 62 | 0 | 6 | 56 |
| 1 | 13 | 14 | 182 | 0 | 6 | 176 |
| 1 | 17 | 22 | 334 | -4 | 30 | 392 |
| 1 | 23 | 26 | 586 | -2 | 30 | 608 |
| 5 | 7 | 8 | 62 | 0 | 6 | 56 |
| 5 | 13 | 14 | 182 | 0 | 6 | 176 |
| 5 | 17 | 14 | 334 | 4 | 30 | 248 |
| 5 | 23 | 22 | 586 | 2 | 30 | 512 |

The certificate computes in F_(p^2)=F_p[s]/(s^2-d) for a verified nonsquare d.
The quadratic character of u+v*s is the Legendre symbol of its norm
u^2-d*v^2. It explicitly includes the correct points at infinity: one for the
odd-degree models, and 1+chi(leading coefficient) for the sextics. It checks
that the polynomial and its derivative are coprime modulo every prime used.

#### Exact torsion of both dual Jacobians

The branch points of D_delta are 0, 1, -5, infinity, and the two roots of
z^2-68z-275. The last quadratic has discriminant 36*159, a nonsquare, so the
first four are rational Weierstrass points and the last pair is conjugate.
The three classes

    [(0,0)-infinity], [(1,0)-infinity], [(-5,0)-infinity]

are independent rational 2-torsion. One way to see independence is the
standard description J[2] as even subsets of six branch points modulo
complement. No nonempty subset of these three classes gives the empty or
full branch subset. Hence (Z/2Z)^3 is contained in the torsion group.

For the upper bound, use the standard good-reduction theorem that rational
ell-primary torsion injects into J(F_p) whenever ell differs from the good
reduction prime p. Counts at p=7 and p=13 suffice here:

    #J(F_7)=56=2^3*7,
    #J(F_13)=176=2^4*11.

For ell other than 7 and 13, the ell-primary order divides both numbers.
The 7-primary part is excluded by the count 176 at p=13; the 13-primary part
is excluded by 56 at p=7. Thus the whole torsion group has order dividing 8.
Combining lower and upper bounds proves, for both twists,

    J(D_delta)(Q)_tors = (Z/2Z)^3.

No assumption about p-primary injectivity at p was needed.

#### Certified infinite-order classes

There are exact rational points

    P_1=(-15/4, 1425/32) on D_1,
    P_5=(-11, 1188) on D_5.

The certificate checks both equations using rational arithmetic. Each point
has nonzero ordinate, so it is not a Weierstrass point.

For an odd-degree genus-two model, infinity is a Weierstrass point and
L(2*infinity) has basis 1,z. If 2*[P-infinity] were zero, a nonconstant
function would have divisor 2*P-2*infinity. It would lie in L(2*infinity),
but a function linear in z has a double zero only at a Weierstrass point.
This contradicts the nonzero ordinate of P. Thus [P-infinity] is not
2-torsion; as all torsion has exponent two, it has infinite order.
Consequently

    rank J(C_1)(Q)=rank J(D_1)(Q) >= 1,
    rank J(C_5)(Q)=rank J(D_5)(Q) >= 1.

These are lower bounds only. The number of small points found in a search
has no role in the proof and supplies no rank upper bound.

#### Additional certified facts about all eight original targets

Let G(j,delta) be the original genus-two curve with factor q_j omitted.
The following good-prime Jacobian orders bound their rational torsion.

| j | primes | orders for delta=1 | orders for delta=5 | torsion order bound |
|---|---|---|---|---:|
| 0 | 7,13,17 | 52,280,292 | 52,112,292 | 4 |
| 1 | 7,13,53 | 56,216,3156 | 56,160,2508 | 4 |
| 2 | 7,13,17 | 56,176,392 | 56,176,248 | 8 |
| 3 | 7,11,47 | 48,112,1980 | 48,112,2556 | 4 |

The code bounds each ell-primary part using only primes p different from
ell, so the displayed bounds do not silently assume injectivity on
p-primary torsion at p. For j=0,1,3 the quadratic-factor branch classes give
an independent (Z/2Z)^2 subgroup. Hence the torsion is **exactly** (Z/2Z)^2
for those six Jacobians. For the original j=2 models, this calculation gives
only a torsion order dividing eight, not the exact torsion structure.

All eight have rank at least one. For j=2 this was proved above. For the
other indices, take a known rational non-Weierstrass point P: one of the two
rational infinities for delta=1, or the point with z=5/2 and positive ordinate
3375/2, 3375/8, 675/2, 675/2 in index order for delta=5. The certificate
checks the latter coordinates. Let iota be the hyperelliptic involution.
The class [P-iota(P)] cannot be zero or have order two. Otherwise 2P and
2*iota(P) would be distinct linearly equivalent effective degree-two
divisors, giving a degree-two map other than the unique hyperelliptic pencil;
that pencil has a doubled fiber only at a Weierstrass point. Since rational
torsion has exponent two for these indices, the class has infinite order.
Thus none of the eight targets provides a rank-zero shortcut.

Finally, each original Jacobian is **Q-simple**. At the indicated good
prime its Frobenius polynomial is irreducible over Q:

- j=0, p=7: T^4+2*T^2+49;
- j=1 and j=2, p=7: T^4+6*T^2+49;
- j=3, p=11: T^4-10*T^2+121.

These formulas hold for both twists because a1=0 at the selected primes.
The certificate tests irreducibility by Gauss's lemma: it excludes every
possible monic integer linear and quadratic factor, using the divisors of
the constant coefficient and coefficient comparisons. An elliptic factor
defined over Q would yield a proper Frobenius factor after good reduction,
contradicting these polynomials. Hence there is no Q-isogeny to a product
of elliptic curves and no Q-defined elliptic quotient of the original
curves. This does **not** prove absolute simplicity over Qbar, and no such
claim is made. In particular, the even Frobenius polynomials used here can
split after extension of the finite field.

#### Remaining rational-point certification

Section 5.13 below performs full two-descent on the selected Richelot duals
and proves both ranks exactly two. This transfers to the original j=2 pair,
so ordinary genus-two Chabauty's rank-below-genus hypothesis fails there.

A complete rational-point certificate on the original curves is still
needed, using sufficient Mordell-Weil data or another justified effective
method. A dual point alone does not provide original-curve points. Every
remaining z must also pass the other center-square, singleton-square and
collision conditions. The original distinct-eleven-integer problem remains
unresolved.

### 5.13 Exact rank two for both Richelot targets

This supplies both a complete **two-Selmer computation** and independent
rational divisor classes proving **both ranks are exactly two**. It does
not determine all rational points.

For

    D_k: Y² = k z(z−1)(z+5)(z²−68z−275),    k = −15 or −3,

let J_k be its Jacobian. The result is

    dim_F2 Sel_2(J_−15/Q) = dim_F2 Sel_2(J_−3/Q) = 5,
    rank J_−15(Q) = rank J_−3(Q) = 2.

The rational two-torsion dimension is three. These exact ranks transfer to the
original j=2 Jacobians by the already established Richelot isogenies.
In particular, the rank is equal to the genus, so this pair does **not**
satisfy the rank condition for ordinary Chabauty–Coleman.

The self-contained certificate is `scripts/check_richelot_two_selmer.py`. Run it with
Python 3. It uses exact integers and rational numbers, no CAS, no numerical
approximations, no external data, and no random search. Its witnesses are
explicit local divisors. Reaching the known local quotient dimensions proves
that they span the **entire local Jacobian image**, not just a collection of
curve-point images.

#### The odd-degree full two-descent model

Set x=kz and y=k²Y. This gives the monic odd-degree model

    C_k: y² = f_k(x)
    f_k(x) = x(x−k)(x+5k)(x²−68kx−275k²).

Put s=√159 and K=Q(s). Its étale algebra is

    L = Q × Q × Q × K,

with the image of x equal to

    T = (0, k, −5k, k(34+3s)).

The usual odd-degree two-descent theorem identifies H¹(Q,J_k[2]) with

    H = ker(N: L*/L*² → Q*/Q*²).

For a non-Weierstrass point P=(x,y), the Kummer map on [P−∞] is

    δ(P) = (x, x−k, x+5k, x−34k−3ks).

At a rational branch point e in {0,k,−5k}, replace the zero coordinate by
f_k′(e), leaving the other coordinates unchanged. These formulas hold over
all completions. The Selmer group is the subgroup of H whose localization
belongs to δ_v(J_k(Q_v)) at every place v.

This is full two-descent, not the fake descent for an even-degree model.
References for the theorem are E. F. Schaefer, *2-Descent on the Jacobians of
Hyperelliptic Curves*, Journal of Number Theory (1995),
[DOI 10.1006/jnth.1995.1044](https://doi.org/10.1006/jnth.1995.1044), and the
odd-degree discussion in
[Poonen–Schaefer, Explicit descent for Jacobians of cyclic covers of the projective line](https://math.mit.edu/~poonen/papers/descent.pdf).

#### The exact global square-class space

Take

    S = {2,3,5,11,19,53}.

The checker computes, by an exact Sylvester determinant,

    disc(f_1) = 2^8 · 3^13 · 5^8 · 11² · 19² · 53,
    disc(f_k) = k^20 disc(f_1).

The second identity follows directly by scaling all five roots by k. Thus
outside S the monic curve has good odd reduction, and every Selmer class has
even valuation in every factor of L. Write Q(S,2) and K(S,2) for the resulting
square-class groups. A Q(S,2) basis is

    −1, 2, 3, 5, 11, 19, 53.

A K(S,2) basis, in the order used by the certificate, is

    −1, 2, 3, 5, 11, 19, s, 13+s, 12+s, 25+2s, 11+s.

Here is a complete justification of the dimension and spanning assertion.
The ring of integers is Z[s] and its discriminant is 636. Minkowski's bound
is √636/2 < 13. The primes 2 and 3 are ramified, 5 and 11 split, and 7 is
inert. Thus the ideal class group is generated by prime ideals over
2,3,5,11. The norm identities

    N(13+s)=10,   N(12+s)=−15,   N(25+2s)=−11

show respectively that each prime-ideal class over 5 is the class of the
prime over 2, each class over 3 is the same class, and the primes over 11
are principal. The prime over 2 has square (2), so the class group is
cyclic of order at most two. Inverting S therefore makes the class group
trivial. This argument needs no class-number software and no fundamental
unit calculation.

There are nine primes of K above S: one each above 2,3,53, and two each
above 5,11,19. The S-unit theorem gives

    dim_F2 K(S,2) = (real unit rank 1) + 9 + (torsion rank 1) = 11.

All eleven listed elements have norm supported on S; the additional norm
identity N(11+s)=−38 is checked. The certificate computes their simultaneous
local square classes at R and the six primes in S and verifies rank eleven.
A square in K would be square at every completion, so this proves global
independence. The eleven elements are consequently a basis by the dimension
calculation, without an assumption about a fundamental unit.

The ambient product Q(S,2)^3 × K(S,2) has dimension 3·7+11=32. Its norm map
to Q(S,2) is surjective, already on its first rational factor. Its norm
kernel therefore has dimension 25. The checker constructs that kernel by
seven exact F2 linear equations.

#### Exact local square-class calculations

At an odd rational prime, Q_p*/Q_p*² is represented by valuation parity and
the Legendre symbol of the unit part. At 2, use valuation parity and the
odd unit modulo 8. At the real place, use sign.

For K, the primes 5,11,19 split. The checker uses the two Hensel lifts of
s²=159 and applies the rational square-class calculation to both components.
Precision is increased until the valuations and unit square classes are
stable; only exact modular arithmetic is used.

At 3 and 53 the extension is ramified, with uniformizer s. A class is
specified by the parity of its uniformizer valuation and the quadratic
character of its normalized residue in F_p.

At 2 the extension is ramified, with uniformizer π=1+s and residue field F2.
For α=a+bs, its π-valuation is v_2(N(α)). After dividing by that power of π,
reduce the unit coefficient pair modulo 8. There are 32 unit pairs, and their
squares are exactly

    (1,0), (5,4), (7,0), (3,4).

These are enumerated in the checker. Membership in these four residues is
also sufficient for a unit to be a square: a remaining factor in 1+8O_K
has valuation at least 6, strictly greater than 2v_π(2)=4, so Hensel's lemma
applies. Thus the unit quotient has dimension three and the full local
square-class group dimension four. The program verifies its multiplication
law exhaustively. It additionally checks the local homomorphism law on every
pair of the eleven global K generators at every relevant place.

#### Local images, with a completeness certificate

For odd p,

    dim J_k(Q_p)/2J_k(Q_p) = dim J_k[2](Q_p).

At 2 one adds g=2. This follows from the structure of the p-adic Lie group:
its sufficiently small open subgroup is Z_p^g, and the finite quotient
contributes the rational two-torsion. The quadratic factor of f_k remains
irreducible at 2,3,53 and splits at 5,11,19. Thus the local dimensions are

| Place | R | 2 | 3 | 5 | 11 | 19 | 53 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Dimension | 2 | 5 | 3 | 4 | 4 | 4 | 3 |

At R all six branch points, including infinity, are real; the real Jacobian
has 2² connected components, giving dimension two.

The following rational x-coordinates specify the additional Q_p points used
as local witnesses. A coordinate need not give a Q-rational point: the
checker verifies f_k(x) is a nonzero square in the indicated local field,
which is exactly what is required. All listed torsion witnesses are the
classes at x=0,k,−5k. A subset suffices in the rows indicated.

| k | Place | Branch-point witnesses | Additional point x values |
|---|---|---|---|
| −15 | R | 0,−15 | none |
| −15 | 2 | 0,−15,75 | −100,−81 |
| −15 | 3 | 0,−15,75 | none |
| −15 | 5 | 0,−15,75 | −100 |
| −15 | 11 | 0,−15,75 | −95 |
| −15 | 19 | 0,−15,75 | −99 |
| −15 | 53 | 0,−15 | −94 |
| −3 | R | 0,−3 | none |
| −3 | 2 | 0,−3,15 | −100,−95 |
| −3 | 3 | 0,−3 | −75 |
| −3 | 5 | 0,−3,15 | −97 |
| −3 | 11 | 0,−3,15 | −100 |
| −3 | 19 | 0,−3,15 | −93 |
| −3 | 53 | 0,−3 | −100 |

Every witness belongs to the genuine local Kummer image. The checker verifies
that their square classes are independent and reach the known dimension in
each row. They therefore span the **whole local Jacobian image**. No claim
that all local divisor classes come from single curve points is needed, and
no omitted degree-two divisor could enlarge these already full images.

#### Linear algebra and the resulting Selmer groups

For each global generator the checker computes its localization, reduces it
modulo the certified local image, and adds the resulting homogeneous linear
conditions. The dimensions of the surviving global space are:

| Conditions imposed | Norm | + R | + 2 | + 3 | + 5 | + 11 | + 19 | + 53 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| k=−15 | 25 | 23 | 18 | 15 | 11 | 8 | 5 | 5 |
| k=−3 | 25 | 23 | 18 | 15 | 11 | 8 | 6 | 5 |

Outside S, good odd reduction identifies the local Kummer image with the
unramified subgroup of H¹(Q_p,J_k[2]). All elements in the global space
already have even valuation outside S. Hence there are no further local
conditions: these are the exact two-Selmer groups, not just supersets.

The program gives an explicit basis by 32-bit exponent vectors. The first
seven bits refer to the Q(S,2) basis in the first rational coordinate, the
next seven to the second, the next seven to the third, and the final eleven
to the K(S,2) basis above. In this convention the five basis vectors are:

    k=−15:
    1107444881, 2160185348, 3523219480, 2445333760, 2445366784.

    k=−3:
    1673724049, 1929412628, 3523219480, 2969669888, 2428719616.

They satisfy all 27 independent equations and are independent. Consequently
both two-Selmer groups have dimension five. Since J_k(Q)[2] has dimension
three, the Mordell–Weil rank is at most 5−3=2.

#### An exact rank lower bound of two

The checker first verifies the known rational points

    (z,Y)=(−15/4,1425/32) for k=−15,
    (z,Y)=(−11,1188) for k=−3,

and checks that each of their global Kummer classes, together with the three
rational branch classes, has dimension four after localization. A further
rational degree-two divisor supplies the missing independent class.

On D_−15 take the divisor cut out by

    u(z)=z²−68z+361,    Y=v(z)=906z−6954.

On D_−3 take the divisor cut out by

    u(z)=z²−20z+361,    Y=v(z)=342z−4446.

Both divisors are rational, although their individual points are conjugate.
The discriminants of u are 3180 and −1044, respectively. The exact identities,
with F=z(z−1)(z+5)(z²−68z−275), are

    −15F−(906z−6954)²
      = (z²−68z+361)(−15z³−60z²+9615z−133956),

    −3F−(342z−4446)²
      = (z²−20z+361)(−3z³+132z²+5379z−54756).

Thus u divides kF−v² and defines an effective rational degree-two divisor E
on D_k. Its Jacobian class is [E−2∞]. These are genuine rational Jacobian
classes, not merely locally soluble covers or points found over an unspecified
extension.

For u(z)=z²+Az+B, the transformed monic polynomial in x is

    u_x(x)=x²+kAx+k²B.

The Kummer image of [E−2∞] is u_x(T), evaluated in the four factors of L.
For the two divisors above the images are respectively

    k=−15: (81225, 66150, 163350, 143100),
    k=−3:  (3249, 3078, 4374, 20412+1296s).

The checker verifies the polynomial divisibilities, evaluates these four
coordinates directly, and checks that adjoining this class to the three
branch classes and the known point class increases their simultaneous local
Kummer rank from four to **five**, for each twist. Global square-class
independence follows from local independence. Consequently

    dim J_k(Q)/2J_k(Q) ≥ 5,
    rank J_k(Q) ≥ 5 − dim J_k(Q)[2] = 2.

Combined with the two-Selmer upper bound, this proves both ranks equal two.
The same exact sequence also shows Sha(J_k/Q)[2]=0 for these two dual
Jacobians. No assumption that the entire Tate–Shafarevich group is finite is
needed for that assertion.

#### What remains open

The two selected original j=2 Jacobians also have exact rank two, by
isogeny. Ordinary genus-two Chabauty–Coleman therefore cannot provide the
proposed rank-one route on this pair. This is a proved obstruction to that
specific route, not to all rational-point methods.

The calculation supplies neither a fully saturated Mordell–Weil basis nor a
complete rational-point set on the original curves. A subgroup generated by
the displayed classes and the rational two-torsion surjects onto J_k(Q)/2J_k(Q)
and hence has finite odd index in J_k(Q), and Section 5.16 below proves that 13 does not divide this index.
Saturation at the other odd primes is not claimed. Moreover, the original curves still require their own
explicit Mordell–Weil data or use of the isogeny correspondence. The original
distinct-eleven-integer PTE classification remains open.

### 5.14 A genus-three quadratic-Chabauty route for the remaining center lifts

This gives a genuinely applicable finiteness criterion, not a rational-point
classification. The genus-two rank-two targets themselves have rational
Néron–Severi rank one. However, two explicit genus-three covers contain every
remaining PTE center lift in the sixth chain mode and have Jacobian rank three and rational
Néron–Severi rank two. They therefore satisfy the usual quadratic-Chabauty
finiteness inequality.

The exact Python 3 checker is `scripts/check_genus_three_qc_targets.py`. It uses the
polynomial and finite-field routines in `scripts/check_genus_two_richelot_certificate.py`,
with no third-party packages. The previously certified genus-two rank-two
result in Section 5.13 remains an input.

#### The genus-two rational Néron–Severi rank is exactly one

Write

    q0=z²+5, q1=4z²+30z+80, q3=z²+275,
    Cδ: δ y²=q0 q1 q3, δ∈{1,5}.

For C1 the Frobenius polynomials at the good primes 17 and 23 are

    P17(T)=T⁴+4T³+30T²+68T+289,
    P23(T)=T⁴+2T³+30T²+46T+529.

For C5 the odd-degree coefficients change sign. The checker independently
counts points over Fp and Fp², verifies these polynomials are irreducible over
Q, and checks their middle coefficient is prime to p. Thus both reductions
are simple ordinary abelian surfaces. By the finite-field endomorphism
theorem, their endomorphism algebras over Fp are the commutative quartic
fields Kp=Q(πp).

Specialization of rational endomorphisms embeds End_Q(J(Cδ))⊗Q into each Kp
and respects the Rosati involution induced by the principal polarization.
On Kp the involution sends πp to p/πp. Its fixed field is generated by
Up=πp+p/πp, satisfying

    Up²−a1 Up+(a2−2p)=0.

The discriminants at 17 and 23 are respectively 32 and 68. Consequently the
two fixed fields are Q(√2) and Q(√17), for either twist.

Because the rational endomorphism algebra embeds in a field, it is itself
commutative. Its Rosati-fixed part is a subfield, contains Q, and embeds in
both of these nonisomorphic quadratic fields. It must therefore be Q.
The standard identification of NS(J)⊗Q with the Rosati-symmetric rational
endomorphisms now gives

    rho_Q(J(C1))=rho_Q(J(C5))=1.

This argument concerns Néron–Severi classes defined over Q, not geometric
Néron–Severi rank over Qbar. It also proves Q-simplicity, since a rational
elliptic factor would contradict irreducibility of either Frobenius
polynomial. For genus two and rank two, the usual sufficient inequality
r<g+rho−1 reads 2<2 and does not apply. This is not a statement that all
nonabelian methods fail.

#### Two covers contain every relevant original center lift

Let Xδ be the smooth projective normalization of

    δ a²=q1(z),   b²=q0(z)q3(z),   δ∈{1,5}.

There are genuine rational maps

    Xδ → Cδ,  (z,a,b) ↦ (z,y=ab),
    Xδ → E0,  (z,a,b) ↦ (z,b),

where E0 is the genus-one quartic b²=(z²+5)(z²+275).
The first map has degree two and is unramified. For example, over Qbar the
cover of the z-line has degree four, with six distinct branch points of
inertia order two: the two roots of q1 and the four roots of q0 q3.
Riemann–Hurwitz gives 2g(Xδ)−2=−8+12=4, hence g(Xδ)=3. The quotient by
(a,b)↦(−a,−b) is Cδ; its degree-two map has no ramification by the same formula.

The established sixth-mode center-square reduction says that a valid finite
center lift has qi(z)=δ ci² for i=0,1,2,3, with δ=1 or 5 and rational ci. It therefore
lifts to Xδ by taking a=c1 and b=δ c0 c3. Thus these two covers suffice as
necessary targets for the PTE problem. They are not asserted to cover all
of Cδ(Q), and a point on Xδ does not itself satisfy the remaining q0
squareclass, q2 square, singleton-root or collision requirements.

The three involution quotients of Xδ have genera 0, 1 and 2: the conic
δa²=q1, E0, and Cδ. Pullback of divisor classes gives a Q-isogeny

    J(Xδ) ~_Q J(Cδ) × E0.

The elliptic quartic E0 is birational to

    E: Y²=X(X²−560X+72900),
    X=2b+2z²+280, Y=2zX.

Its rank is exactly one, as already proved by two-isogeny descent. The
checker repeats that argument: the exact descent images on E and its
2-isogenous partner Y²=X(X²+1120X+22000) are {1,5} and {1,−5,−11,55}.
The displayed classes are all witnessed by rational points. The point
(180,900) on E has doubled x-coordinate 2025/4, independently confirming
it is nontorsion by the integral-coordinate torsion theorem.

Therefore

    rank J(Xδ)(Q)=2+1=3.

Q-simplicity of J(Cδ), and the difference in dimensions, imply
Hom_Q(J(Cδ),E0)=0. Néron–Severi rank is invariant under Q-isogeny, so

    rho_Q(J(Xδ))=rho_Q(J(Cδ))+rho_Q(E0)=1+1=2.

#### Explicit hyperelliptic equations and maps

For each δ define Nδ(t), Dδ(t), Mδ(t) by

    δ=1:
    N1=80−t²,
    D1=4t−30,
    M1=2t²−30t+160;

    δ=5:
    N5=25t²−120t+80,
    D5=10t²−8,
    M5=−60t²+100t−48.

They satisfy the polynomial identity

    δ Mδ²=4Nδ²+30NδDδ+80Dδ².

Consequently Xδ has the hyperelliptic model

    W²=(Nδ²+5Dδ²)(Nδ²+275Dδ²),
    z=Nδ/Dδ, a=Mδ/Dδ, b=W/Dδ²,
    y=Mδ W/Dδ³.

The conic parametrizations are birational. For δ=1 their inverse is
t=a−2z. For δ=5 it is t=(a−6)/(z−5/2), with the base point and tangent
interpreted on the smooth projective conic.

For δ=1 the two quartic factors are

    t⁴−80t²−1200t+10900,
    t⁴+4240t²−66000t+253900.

For δ=5 they are

    1125t⁴−6000t³+17600t²−19200t+6720,
    28125t⁴−6000t³−25600t²−19200t+24000.

The checker verifies degree eight and squarefreeness modulo 7 and 13,
which in particular proves squarefreeness over Q and genus three. It also
checks the quotient identities and rational base points:

- X1 has its two rational hyperelliptic points at t=∞. The other two
  points lying above z=∞ have t=15/2 and W=±9025/16.
- On X5, t=5/6 gives z=5/2, a=6 and W=±9025/144, hence b=±225/4.
  Its two points at t=∞ have z=5/2, a=−6 and b=±225/4.

The rational formulas extend to morphisms of smooth projective curves.
A future point-completeness computation must include these points and
must not silently drop denominator-zero charts.

#### The quadratic-Chabauty criterion really does apply

For both Xδ,

    r=3, g=3, rho_Q=2, so r=3<g+rho_Q−1=4.

The Mordell–Weil-restricted quadratic-Chabauty finiteness lemma therefore
applies to these targets.
This is the standard criterion of Balakrishnan–Dogra, restated with its
rational Néron–Severi convention in
[Unlikely intersections and the Chabauty–Kim method over number fields](https://doi.org/10.1007/s00208-023-02638-2).
Here the restricted (modified) depth-two Selmer set uses the image of
J(Q) tensor Qp at depth one, as specified in Dogra’s Remark 1.5, which refers to
Balakrishnan–Dogra’s Definition 2.2(c). The analogous statement for the unrestricted Selmer set would
require additional control of the p-primary Tate–Shafarevich group. No
finiteness of Sha[13^infinity] is assumed or proved; Sha[2]=0 for a dual
Jacobian does not establish it.

This numerical finiteness statement is distinct from the extra hypotheses
and arithmetic input needed by a particular explicit p-adic-height algorithm.

There is also a concrete rational symmetric correspondence. If f:Xδ→E0
is the degree-two elliptic quotient, put

    Z=3 f* f_*−2 id  on J(Xδ),

where f* denotes pullback. The operator f* f_* acts by 2 on the elliptic
factor and by 0 on the genus-two factor. Thus Z acts by 4 and −2,
respectively; its trace on holomorphic differentials is 4−2−2=0.
It is nonzero and Rosati symmetric, and supplies a rational trace-zero
Néron–Severi class of the type used by quadratic Chabauty. No associated
height function has been calculated here.

A concrete good ordinary prime is p=13. The elliptic factor has 12 points
and trace 2, while each genus-two factor has Frobenius polynomial
T⁴+6T²+169. The genus-three Frobenius polynomial is therefore

    (T²−2T+13)(T⁴+6T²+169)
    =T⁶−2T⁵+19T⁴−12T³+247T²−338T+2197.

Its middle coefficient is prime to 13, and the explicit models have good
reduction there. The checker verifies both facts.

#### What this does not yet certify

This is a new applicable route, not a completed quadratic-Chabauty
calculation. In particular, it supplies no finite p-adic candidate list,
height bound, or complete rational-point set. An explicit implementation
still needs sufficient Mordell–Weil data on Xδ, the relevant index and
saturation control (the logarithm rank and 13-saturation are now certified
in 5.16), local and global height calculations, rigorous p-adic root
isolation, and an appropriate final Mordell–Weil sieve or equivalent proof
of rational-point completeness. The surviving z-values must then satisfy
all remaining PTE lifting and collision conditions.

The original distinct-eleven-integer PTE classification remains open.

### 5.15 Exact de Rham input for the two genus-three quadratic-Chabauty targets

This supplies an explicit rational de Rham basis, cup-product matrix,
holomorphic-first symplectic basis, and genuine trace-zero correspondence.
It does not compute a Frobenius matrix, a p-adic logarithm or height, a
quadratic-Chabauty root, or a rational-point classification.

The deterministic checker is `scripts/check_genus_three_qc_de_rham.py`, requiring
Python 3 and SymPy. Run it with

    python3 scripts/check_genus_three_qc_de_rham.py

The existing genus-three certificate supplies the smooth curves, quotient
maps, isogeny decomposition, and Frobenius polynomials used below. The
checker verifies the additional identities exactly, not numerically.

#### Curves and an explicit quotient involution

Write

    E0: b^2=(z^2+5)(z^2+275),
    C_delta: delta*y^2=(z^2+5)(4z^2+30z+80)(z^2+275).

For delta=1 and 5, use the degree-eight model

    X_delta: W^2=(N^2+5D^2)(N^2+275D^2),
    z=N/D, a=M/D, b=W/D^2, y=MW/D^3,

with

    delta=1: N=80-t^2, D=4t-30, M=2t^2-30t+160;
    delta=5: N=25t^2-120t+80, D=10t^2-8, M=-60t^2+100t-48.

The degree-two elliptic quotient has deck involution

    sigma_1(t,W)=((15t-160)/(2t-15), 9025W/(2t-15)^4),
    sigma_5(t,W)=((24-25t)/(25-30t), 9025W/(25-30t)^4).

Direct substitution verifies the curve identities and sigma^2=1. These
maps fix z,b and send a,y to -a,-y. Their rational formulas extend across
their coordinate poles on the smooth projective curves. The norm identity
for the elliptic quotient f gives f* f_*=1+sigma on the Jacobian. Hence the
previous trace-zero endomorphism is

    Z=3f* f_*-2=1+3sigma.

In the ordered basis (dt/W,t dt/W,t^2 dt/W), the rows of the following
matrices express pullbacks by sigma:

    S1 = [[45/19,   -12/19,    4/95],
          [480/19, -109/19,    6/19],
          [5120/19,-960/19,   45/19]],

    S5 = [[125/19, -300/19,  180/19],
          [120/19, -269/19,  150/19],
          [576/95, -240/19,  125/19]].

Both square to the identity and have eigenvalues 1,-1,-1.

#### Second-kind quotient differentials and the sign convention

On E0 set

    e0=dz/b,  e1=z^2 dz/b.

On C_delta set

    c0=dz/y,  c1=z dz/y,
    c2=(z^3+15z^2/4) dz/y,
    c3=(z^4+4125z^2/32) dz/y.

Every differential has no finite poles. At a finite branch point, dz/y
is regular in the uniformizer y. To check the infinities, put u=1/z and
write f=(z^2+5)(4z^2+30z+80)(z^2+275). Then

    (u^6 f(1/u)/4)^(-1/2)
       =1-15u/4-4125u^2/32+64725u^3/128+O(u^4).

Since dz/y equals -(sqrt(delta)/2)u times this series times du at one
infinity, the displayed corrections in c2,c3 cancel their residues. The
other sheet changes the signs of all differentials, so its residues also
vanish. The e forms have zero residues by the analogous even quartic
expansion. Thus all six quotient forms are of the second kind.

Fix the signed cup convention

    <alpha,beta> = sum_P Res_P(F_alpha beta),  dF_alpha=alpha,

where the sum is over the poles on the smooth projective quotient curve.
Termwise local primitives are sufficient here; their arbitrary constants
do not matter because beta has zero residue. In this convention the
elliptic cup matrix is

    CE = [[0,2],[-2,0]],

and the genus-two matrix in (c0,c1,c2,c3) is

    CC = [[0,        0,          0,          delta/4],
          [0,        0,          delta/2,   -45delta/16],
          [0,       -delta/2,     0,          525delta/2],
          [-delta/4, 45delta/16, -525delta/2, 0]].

The checker obtains every entry from independent Laurent expansions,
including <c2,c3>=525delta/2. It verifies det(CC)=delta^4/64.

#### A full de Rham basis on X_delta

From now on each e or c denotes its pullback to X_delta. Each quotient
map has degree two, so cup products of pulled-back forms are multiplied
by two. The e classes have sigma eigenvalue +1, while the c classes have
eigenvalue -1. Cup-product invariance under sigma makes the cross-pairings
zero. Therefore, in the order (e0,e1,c0,c1,c2,c3),

    CX = diag(2CE,2CC),       det(CX)=4delta^4 != 0.

This proves that the six classes are independent modulo exact
differentials. Since dim H^1_dR(X_delta)=6, they are a basis. In this
basis sigma=diag(1,1,-1,-1,-1,-1) and

    Z=diag(4,4,-2,-2,-2,-2).

For explicit rational expressions in t,W, put R=N'D-ND'. One has R=lambda M
with lambda=-2 for delta=1 and lambda=-20 for delta=5. The pullbacks are

    e0=R dt/W,
    e1=(RN^2/D^2) dt/W,
    c0=lambda D dt/W,
    c1=lambda N dt/W,
    c2=lambda (N^3/D^2+15N^2/(4D)) dt/W,
    c3=lambda (N^4/D^3+4125N^2/(32D)) dt/W.

The first, third and fourth forms, namely e0,c0,c1, are independent
polynomial multiples of dt/W of degree at most two, so they form a
holomorphic basis. The expressions with powers of D in their denominators
are meromorphic second-kind differentials, as guaranteed by pullback;
these coordinate poles must not be discarded by an implementation.

#### Holomorphic-first symplectic normalization

The ordered basis

    nu=(e0, c0, c1, e1/4,
        (2c3+45c2/4+1050c1)/delta,
        c2/delta)

has cup matrix

    J=[[0,I3],[-I3,0]].

The checker verifies B CX B^T=J for the displayed rational basis-change
matrix B. In this basis

    Z=diag(4,-2,-2,4,-2,-2).

It is self-adjoint for the cup pairing, and its trace on holomorphic
differentials is 4-2-2=0. Every displayed basis-change coefficient is
13-integral. This does not alone certify integrality or the precision
behavior of an eventual computed Frobenius matrix.

When importing these data into another package, distinguish the
endomorphism matrix from its associated alternating correspondence form.
For example, in the convention where an alternating matrix K is converted
to its endomorphism by KJ, the input is K=ZJ^(-1), not Z itself. A different
cup sign or row/column convention requires the corresponding conversion.

#### An independent polynomial projector at p=13

The previously certified Frobenius factors are

    PE(T)=T^2-2T+13,
    PC(T)=T^4+6T^2+169.

Define

    eE(T)=(2T+9)PC(T)/2704,       Zpoly(T)=6eE(T)-2.

Exact polynomial division gives eE=1 modulo PE and eE=0 modulo PC.
Consequently eE^2=eE modulo PE PC, and Zpoly has residues 4 and -2 on
the elliptic and genus-two factors. If F is a Frobenius matrix on the
compatible cohomology, eE(F) is the elliptic projector and Zpoly(F) is
the already constructed rational correspondence Z. This conclusion uses
the known quotient decomposition, not an assumption that arbitrary
Frobenius polynomials lift to rational endomorphisms.

The denominator is 2704=16*13^2. Its division can lose two 13-adic digits;
total error also depends on matrix-entry valuations, evaluation, and basis
changes. The exact diagonal correspondence above avoids this division
when supplying Z directly. No Frobenius matrix has been computed here.

For reference, direct integer discriminants of the two degree-eight
polynomials are

    delta=1: 2^52*3^33*5^30*11^2*19^12*53,
    delta=5: 2^108*3^33*5^58*11^2*19^12*53.

Both reductions modulo 13 are squarefree and have ten affine points and
two points at infinity, with no affine Weierstrass points. Thus a complete
13-adic computation must account for twelve residue discs on each curve.

#### Scope of the new certificate

These are exact, reusable cohomological inputs, not a quadratic-Chabauty
calculation. Section 5.16 below certifies logarithm rank and 13-saturation.
Still required are suitable Mordell-Weil divisor data, any other index
control needed by the final sieve, and compatible
local and global heights including all away-from-13 contributions,
certified power-series root isolation on every chart, and a final
Mordell-Weil sieve or equivalent completeness argument. The original PTE
classification remains open.

### 5.16 Certified 13-adic logarithm rank and 13-saturation

This section computes genuine p-adic arithmetic inputs, using exact finite
rings and proved error bounds. It proves that both selected genus-two
logarithm matrices have full rank at 13, that the known subgroups are
13-saturated, and that both sixth-mode genus-three covers have logarithm
rank three. It does not compute p-adic heights, quadratic-Chabauty roots,
or a complete rational-point set.

Two dependency-free certificates accompany the proof:

    python3 scripts/check_thirteen_logarithms.py
    python3 scripts/check_thirteen_saturation.py

The second uses the adjacent finite-field utilities from
check_genus_two_richelot_certificate.py. It gives an independent
finite-field proof of saturation, not a point search.

#### Models and the global classes being tested

Use the monic models from 5.13,

    J_k = Jac(y^2=f_k(x)),
    f_k=x(x-k)(x+5k)(x^2-68kx-275k^2),  k=-15,-3.

The differential basis is eta=(dx/(2y), x dx/(2y)).
The rational classes P=[(x,y)-infinity] are

    k=-15: (x,y)=(225/4,320625/32);
    k=-3:  (x,y)=(33,10692).

The second class Q=[div(u,v)-2 infinity] is specified by

    k=-15: u=x^2+1020x+81225,  v=-13590x-1564650;
    k=-3:  u=x^2+60x+3249,    v=-1026x-40014.

Their exact global identities and independence were proved in 5.13.
The calculations below concern these same classes, not newly assumed
local points.

Both monic models have good reduction at 13. Direct square counting in
F_13 and F_169 gives N1=14, N2=182, and hence

    #J_k(F_13)=(N1^2+N2)/2-13=176.

The logarithm certificate recomputes these counts rather than importing
a stored Frobenius polynomial. In particular [176]D reduces to zero for
every D in J_k(Q_13), and 176 is a 13-adic unit.

#### Exact arithmetic in an unramified extension

Let O be the integers of Q_13(s), where s^2=2. Since 2 is a nonsquare
modulo 13, this is an unramified quadratic extension. All calculations
are performed in O/13^2 O, represented as

    (Z/169Z)[s]/(s^2-2).

A denominator is inverted only when its norm is a unit modulo 169.
The checker fails rather than dividing by a nonunit.

Instead of reducing a divisor close to the Jacobian identity in a
potentially singular Mumford chart, translate it by a fixed effective
degree-two divisor R. The following coordinate pairs specify its two
points modulo 169:

    k=-15:
      (6,20), (4+8s,35+83s);

    k=-3:
      (12+12s,127+49s), (1+7s,11+166s).

Their ordinates and their x-coordinate differences are units. The curve
equations are checked modulo 169, so Hensel's lemma gives actual points
over O with these residues. R is auxiliary and need not be defined over
Q_13. The difference [R+176D]-[R] is the original Q_13-rational class
176D, so no field-trace normalization is introduced.

For D=P and Q separately, the checker performs 176 additions to R.
For coprime Mumford polynomials u1,u2, it verifies a Bezout identity
a*u1+b*u2=1, uses the Chinese remainder formula for the ordinate
polynomial, then performs the usual Cantor reduction. Every division
has a unit leading coefficient. Every reduced u is required to have
degree exactly two, with monic leading coefficient. In particular, an
apparent degree drop modulo 169 is rejected: silently dropping a point
near infinity would invalidate the computation.

These assertions keep every operation in a valid reduced Mumford chart
and make the finite-ring calculation the reduction of local Jacobian
arithmetic. The endpoint S=R+176D is checked to have the same Mumford
reduction modulo 13 as R.

#### Tiny integrals with a rigorous error bound

Write the points of R as (x_i,y_i). The endpoint polynomial u_S has a
unique root x_i' in each of their residue discs. One Newton step gives it
modulo 169:

    x_i'=x_i-u_S(x_i)/u_S'(x_i).

The derivative is a unit. The checker verifies u_S(x_i')=0,
y_i'=v_S(x_i'), the curve equation, and y_i'=y_i modulo 13.

Put h_i=x_i'-x_i, so h_i belongs to 13O. Because y_i is a unit, x is a
local parameter and both integrands have integral power-series
coefficients. For a term a_n h_i^(n+1)/(n+1), n>=1, its valuation is at least

    n+1-v_13(n+1) >= 2.

Thus the first-order formulas give rigorous congruences modulo 13^2:

    integral eta_0 = h_i/(2y_i) mod 13^2,
    integral eta_1 = x_i h_i/(2y_i) mod 13^2.

Summing over the two discs computes log(176D). This is the abelian
logarithm of S-R; it does not require a global Coleman integration path.
The resulting coordinates are checked to lie in Q_13 modulo the stated
precision.

The certified rows log(176D)/13 modulo 13 are

| k | D=P | D=Q | determinant |
|---|---|---|---:|
| -15 | (8,7) | (7,12) | 8 |
| -3 | (11,7) | (5,11) | 8 |

Dividing by 176, whose inverse modulo 13 is 2, gives

| k | log(P)/13 | log(Q)/13 | determinant |
|---|---|---|---:|
| -15 | (3,1) | (1,11) | 6 |
| -3 | (9,1) | (10,9) | 6 |

Consequently, for either curve,

    v_13(det(log(P),log(Q)))=2.

In particular the two logarithms are linearly independent over Q_13.
These are proved congruences with an error bound, not floating-point
approximations or inferences from algebraic Mordell-Weil rank.

#### Saturation from the formal logarithm

On these good odd-degree models eta is an integral basis of invariant
differentials. The formal logarithm identifies the kernel J_1(Q_13)
of reduction with 13Z_13^2. As #J(F_13)=176 is prime to 13,

    log J(Q_13)=13Z_13^2.

Indeed multiplication by 176 brings every local point into J_1, and
division by this unit preserves that lattice. Conversely J_1 is already
contained in J(Q_13).

The normalized determinant 6 shows that log(P),log(Q) span this entire
lattice. If 13R=aP+bQ for a rational point R of the Jacobian, reduction
of logarithms modulo 13^2 forces a=b=0 modulo 13. Subtracting the
resulting integral combination of P,Q leaves a point killed by 13.
There is no local 13-torsion: neither the reduction group of order 176
nor the formal group has it. Thus R belongs to the subgroup generated
by P,Q. This proves 13-saturation. Adjoining the previously known
two-primary torsion does not change the conclusion.

This settles saturation at 13 only. It does not prove saturation at
every odd prime or that the displayed generators form a full
Mordell-Weil basis.

#### An independent finite-field saturation certificate

There is a second proof not using p-adic integration. At the good primes
181 and 199, direct counts on both monic curves give

| prime ell | #C(F_ell) | #C(F_(ell^2)) | #J(F_ell) |
|---:|---:|---:|---:|
| 181 | 192 | 32970 | 34736=13*2672 |
| 199 | 200 | 39854 | 39728=13*3056 |

Let pi_ell(D)=[#J(F_ell)/13]D in the finite Jacobian. Exact Cantor
arithmetic gives the following relations; every stated nonzero
projection has order 13:

    k=-15:
      ell=181: pi(Q)=10*pi(P), pi(P) nonzero;
      ell=199: pi(P)=0, pi(Q) nonzero.

    k=-3:
      ell=181: pi(Q)=7*pi(P), pi(P) nonzero;
      ell=199: pi(Q)=3*pi(P), pi(P) nonzero.

A rational 13-divisible combination must have zero image under both
maps. The resulting coefficient rows over F_13 are respectively

    (1,10),(0,1),       and       (1,7),(1,3).

Their determinants are nonzero. Hence no nonzero combination of P,Q
modulo 13 can be globally divisible by 13. The certificate recomputes
the counts and projections, checks the Mumford curve identities after
every operation, and checks that the full group orders kill P and Q.
An independent implementation confirmed every projected divisor and
relation.

#### The elliptic logarithm and transfer to both genus-three covers

For E: y^2=x(x^2-560x+72900), use its certified generator U=(180,900)
and differential dx/(2y). Exact rational addition gives t=-x/y at 12U:

    t = 4061199830818098935573944975439863124829948
        /31468337648588340182480022250211888181370801.

Here v_13(t)=1 and t/13=12 modulo 13. The elliptic formal logarithm is
t+O(t^2), so

    log(12U)/13=12 mod 13,
    log(U)/13=1 mod 13.

Thus the elliptic logarithm is nonzero, and the product J_k x E has
logarithm rank three. Combine the Richelot isogeny J_k -> J(C_delta)
with the quotient-pullback isogeny J(C_delta) x E -> J(X_delta).
The quartic E0 and the displayed cubic E are birational, as in 5.14.
For the quotient-product map, the two norm maps compose with the
pullbacks as multiplication by two on each factor; the cross terms
vanish by the distinct involution characters. Its kernel is therefore
two-primary. The Richelot map is likewise a power-of-two isogeny.
These maps induce isomorphisms on 13-adic Lie algebras and on the
Mordell-Weil groups localized at 13.

It follows that both genus-three covers X_1 and X_5 have logarithm
rank three, and the p-adic closure of their rational Jacobian points
has finite index in the local group. The abstract image of the displayed
product subgroup is 13-saturated as well. Explicitly, for a power-of-two
isogeny phi with a companion psi satisfying phi psi=psi phi=[m],
if 13R=phi(h) and H is 13-saturated, then psi(R) belongs to H.
Choosing integers u,v with 13u+mv=1 gives
R=phi(u h+v psi(R)), proving saturation of phi(H). These conclusions do not require assuming a
nonzero p-adic regulator from the algebraic rank alone. They also do
not assert that explicit divisor coordinates for the transported
classes on X_delta have already been computed. The determinant residue
6 is for the specified dual genus-two bases only; no target-basis
logarithm determinant is asserted.

This verifies the full-logarithm-rank hypothesis needed by the
corresponding explicit quadratic-Chabauty approach, beyond the
numerical finiteness criterion proved in 5.14.

#### Remaining work

The exact de Rham inputs, logarithm rank and 13-saturation are now
certified. Still missing are the actual local/global height functions
and all permitted away-from-13 contributions (the elliptic factor is
subsequently computed in 5.17), suitably represented
transported divisors, rigorous isolation of every quadratic-Chabauty
root, and a final sieve or equivalent completeness argument with all
required index information. No complete rational-point set, valid
size-eleven integer pair, or general impossibility theorem is claimed.

### 5.17 Certified elliptic-factor 13-adic height and finite local terms

This computes an actual Coleman–Gross height on the elliptic factor of both
genus-three targets. It also computes its away-from-13 local contribution,
including the nonzero bad-prime corrections. It does not compute the
genus-two height block, the full genus-three local height attached to the
trace-zero correspondence, a quadratic-Chabauty root, or a rational-point
classification.

The reproducible checker is `scripts/check_elliptic_height13.py`; it requires only
Python 3 and SymPy. No Sage, Magma, PARI, installed package, remote execution,
or unproved numerical height estimate is used.

#### 1. Curve, normalization, and the actual result

Use the globally minimal equation

    E: Y² = X³ − 560X² + 72900X,
    P = (180,900), T = (0,0), omega = dX/(2Y).

The existing rank-one generator certificate proves E(Q)=Z P direct-sum <T>.
The computations of P itself below do not need this generator theorem;
only the assertion that the component table covers all rational points does.

The branch of log13 is normalized by log13(13)=0. The height character has
component log13 at 13 and sends a uniformizer ell to −log13(ell) at ell!=13.
Choose the complementary de Rham line

    W_E = span([(X−280)omega]).

This is a specified rational splitting, not an assertion that W_E is the
unit-root splitting. Write h_W for the symmetric Coleman–Gross pairing,
without dividing its character by 13. The certified result is

    h_W(P,P) = 59106874840 + O(13^10).

Equivalently its base-13 digits, starting at 13^0, are

    [0,7,10,12,12,6,12,5,7,5].

In particular v13(h_W(P,P))=1. Its half-height is

    lambda_W(P) = h_W(P,P)/2 = 29553437420 + O(13^10).

The formal elliptic logarithm and the rank-one quadratic coefficient are

    log_omega(P) = 87056567455 + O(13^10),
    13 h_W(P,P)/log_omega(P)^2 = 6699432669 + O(13^9).

Thus v13(log_omega(P))=1. The coefficient h_W(P,P)/log_omega(P)^2 is
known with absolute error O(13^8), not O(13^10); dividing by the square of
the logarithm loses two digits. For Q=nP+epsilon T, h_W(Q,Q)=n²h_W(P,P).

#### 2. Compatibility with the previously chosen quotient splitting

For the quartic E0: b²=(z²+5)(z²+275), use

    X=2b+2z²+280, Y=2zX,
    e0=dz/b, e1=z²dz/b.

Direct differentiation gives

    omega=e0/2,
    X omega=e1+140e0+dz,
    [e1]=[(X−280)omega].

Therefore W_E is exactly the elliptic complementary line in the existing
holomorphic-first genus-three symplectic splitting, not a new incompatible
choice. Under the signed cup convention of that certificate,
<omega,(X−280)omega>=1. The differential denoted eta in Balakrishnan–Besser
Theorem 4.1 is its negative, since their normalization is eta cup omega=1.

The origin O_E corresponds to the quartic infinity with b/z² tending to +1;
the other infinity maps to T. For reference, (z,b)=(5/2,−225/4) maps to P,
and (5/2,+225/4) maps to −P+T. A choice of base point on X_delta must still
be accounted for when using f(Q)−f(base) in a curve-level height.

#### 3. Global minimality, including the prime 2

The exact invariants are

    Delta=2^12*3^12*5^7*11, c4=1518400=2^6*5^2*13*73.

At 3,5,11 the valuation of c4 is less than 4, proving minimality; all other
odd primes have good reduction. At 2 the only possible discriminant-lowering
scale is u=2. An integral new model would arise from

    X=4X'+r, Y=8Y'+4sX'+t.

Its coefficients force s in Z2, t in 4Z2, and

    (−560+3r−s²)/4 in Z2,
    (72900−1120r+3r²−2st)/16 in Z2.

The first condition also forces r in Z2. If s is even, it forces r=0 mod4;
the numerator in the second condition is then 4 mod16. If s is odd, it
forces r=3 mod4, and that numerator is odd. Both are contradictions.
Thus the original model is minimal also at 2. The checker independently
exhausts the residue classes needed for this contradiction.

#### 4. Exact away-from-13 contribution of P

Normalize each finite local half-height lambda_ell to equal
v_ell(d(Q))*log13(ell) on the nonsingular subgroup E_0(Q_ell), where d(Q)
is the positive square-root denominator of X(Q). With this normalization
the local Neron multiplication relation is

    lambda_ell(mQ) = m² lambda_ell(Q)
                      + v_ell(psi_m(Q))*log13(ell).

Here psi_m is the usual division polynomial. These are the local functions
with the discriminant-dependent constant removed; in particular no
additional Delta/12 term is to be added to the formulas below.

Exact rational group arithmetic gives

    X(12P)=87327675552949665650305747958870915094857521
            /188867320255610078124863022039459027838224,
    d(12P)=434588679391916602332,
    t(12P)=−X(12P)/Y(12P)
          =4061199830818098935573944975439863124829948
            /31468337648588340182480022250211888181370801.

At ell=2,3,5,11 the respective valuations are

    v_ell(X(12P)) = −4, −6, 0, 0,
    v_ell(Y(12P)) = −6, −9, 0, 0,
    v_ell(d(12P)) =  2,  3, 0, 0,
    v_ell(psi_12(P)) = 74,123,90,0.

The negative X valuations at 2 and 3 put 12P in the formal group. At 5
and 11 the checker reduces the actual coordinates and verifies that at
least one partial derivative of the Weierstrass equation is nonzero.
Thus 12P lies in E_0 at every bad prime, proved directly without assuming
any Tamagawa numbers. Its parameter at 13 has valuation one.

Dividing the difference of the last two valuation lists by 144 proves

    lambda_away(P) = −1/2 log13(2) −5/6 log13(3) −5/8 log13(5)
                   =62652731752+O(13^10).

There are no other nonzero finite contributions for the integral point P.
The coefficient at 11 is zero. The displayed sum is 65 mod169, so it is
provably nonzero. The away contribution to h_W(P,P) is twice this sum.

#### 5. Away-from-13 local terms for every rational point on E

Let Q=nP+epsilon T != O_E, with epsilon in {0,1}, and let d=d(Q).
Define these periodic correction tables, indexed starting at zero:

    kappa2 = [0,1/2],
    kappa3 = [0,5/6,4/3,3/2,4/3,5/6],
    kappa5 = [0,5/8,1/2,5/8].

Then the exact identity is

    lambda_away(Q) = log13(d)
      −kappa2[(n+epsilon) mod2] log13(2)
      −kappa3[(n+3epsilon) mod6] log13(3)
      −kappa5[(n+2epsilon) mod4] log13(5).

Here d can contain 13: its 13-power contributes zero to log13(d). No
factorization of d is required to evaluate the formula.

For the proof, the correction lambda_ell(Q)/log13(ell)−v_ell(d(Q)) is
constant on each component of the Neron model. The checker proves directly
that the class of P has order 2,6,4 at 2,3,5 respectively: these multiples
are nonsingular, and all earlier positive multiples are singular. It also
checks that P+T is nonsingular at 2, 3P+T at 3, and 2P+T at 5. Hence the
component classes of T are respectively P,3P,2P. The local multiplication
relation applied to psi_n(P), for n up to 6, gives exactly the tables above.
At 11 both P and T are nonsingular, so all rational points are in E_0 there.

This argument concerns component classes attained by E(Q). It does not
claim to determine the entire geometric component group at every prime.

##
### 5.30. Complete away-prime targets and certified infinity projection

The follow-up proofs [COMPLETE_AWAY_TARGETS.md](COMPLETE_AWAY_TARGETS.md) and
[INFINITY_PROJECTION.md](INFINITY_PROJECTION.md) close the two local gaps at 2 and 3
and compute the infinity de Rham projection modulo 13³. Both received independent
mathematical and computational review.

For the published correspondence and rho=h13−global, the exact away target sets are
Ω1={0,−1/6,−2/3} log13(3) and
Ω5={0,1/2,2/3,−4/3} log13(3).
Their residues modulo2197 are {0,1482,1534} and {0,2145,663,871}.
The dyadic contribution vanishes by an analytic two-disc argument; no full wild
regular model is claimed. The 3-adic proof includes exhaustive local accessibility
and an independent intersection calibration fixing the scalar and sign.

The certified infinity projections and holomorphic corrections use the fixed
compatible complement, with a rigorous omitted-Frobenius-tail bound. Executable
certificates are scripts/check_complete_away_targets.py and
scripts/check_infinity_projection.py.

These results concern only the two sixth-mode genus-three curves. The moving
normalized QC functions, certified root isolation and final sieve remain open.
No classification of rational points or solution of the original PTE problem is claimed.

## 6. The local 13-adic double integral and its error bound

Put t=−X/Y and write X=u(t)/t². The equation gives

    u²+(−560t²−1)u+72900t^4=0, u(0)=1,
    u=(1+560t²+sqrt((1+560t²)²−291600t^4))/2,
    omega=w(t)dt, w=1−t u'/(2u).

The series u and w are even, have 13-integral coefficients, and converge
for t in 13Z13. Let H and g be the uniquely normalized odd/even primitives

    H'=(−X+280)w, H=1/t+O(t),
    g'=Hw−1/t, g(0)=0,
    log sigma(t)=log13(t)+g(t).

This gives (d/dlog_omega)^2 log sigma=−X+280. It is enough to work with
log sigma; no assertion of integral coefficients for a canonical unit-root
sigma function is being made. The terms needed here are

    g(t)=−140t² +(169525/3)t^4 −(310492000/9)t^6
         +(522868628750/21)t^8 −(177457835600000/9)t^10 + ... .

Because Xw is 13-integral apart from its t^−2 leading term, the two
termwise integrations give the conservative coefficient bound

    v13([t^n]g) >= −2 floor(log_13(n+1)).

For v13(t)>=1, every omitted term with n>=12 consequently has valuation
at least n−2floor(log_13(n+1))>=10. The same inequality tends to infinity,
so it also proves convergence and a bound on the entire tail. Truncation
at degree ten therefore has absolute error O(13^10). The checker computes
through degree twenty and also verifies agreement with degree ten.

The rational logarithms are computed by removing the power of 13 from q
and using log13(q)=log13(q^12)/12 for the resulting 13-adic unit. The
ordinary log(1+z) series has v13(z)>=1; its k-th term has valuation at least
k−v13(k). The checker retains k=1,...,12. All omitted terms have valuation
at least 12, and its two guard digits accommodate finite-series divisions.
Thus each displayed logarithm is rigorously known modulo 13^10.

Balakrishnan–Besser, Theorem 4.1 identifies the local self-pairing for this
splitting with tau=−2log sigma. Corollary 4.3 then gives

    h_W(P,P) = (2/144)(log13(d(12P))−log sigma(t(12P))).

The values used, all modulo 13^10=137858491849, are

    log13(d(12P)) = 82740515494,
    log13(t(12P)) = 53322772061,
    g(t(12P))     = 47336002272.

Division by 144 loses no 13-adic digits. This proves the stated height.
The checker independently obtains exactly the same value through 24P.
It also uses 6P, where the finite contribution has the additional term
−(1/2)log13(5); this third calculation agrees only with that term retained.

#### 7. How much of the genus-three height is now available

Let f:X_delta->E0->E be the elliptic quotient and keep the direct-sum
splitting specified in the de Rham certificate. Functoriality and degree(f)=2
give, for U=f*P in J(X_delta)(Q),

    h_X(U,U)=2h_W(P,P),
    h_X(U,ZU)=8h_W(P,P)=59279523173+O(13^10),

where Z=3f*f_*−2 acts by 4 on this factor. These statements use the
polarization-compatible pullback and the ordinary bilinear height, not a
bare transport of matrix entries through an unspecified isogeny basis.

This supplies the rank-one elliptic block of the global height and every
elliptic finite local term. It does not supply the rank-two genus-two block
or its away-from-13 component corrections. A full curve-level QC local
function additionally requires an explicit basepoint and the normalization
of the height attached to Z; the elliptic tables must not be labeled as
the complete away contribution for X_delta. In particular nothing here
justifies setting the genus-three contributions at 11,19,53 or other bad
primes to zero.

#### References and standard inputs

- J. Balakrishnan and A. Besser, *Coleman–Gross height pairings and the
  p-adic sigma function*, Theorem 4.1 and Corollary 4.3:
  https://arxiv.org/abs/1201.6016 . These apply to the specified complement;
  their additional unit-root specialization is not used.
- J. Silverman, *Computing heights on elliptic curves*, Mathematics of
  Computation 51 (1988), 339–358, for finite local height normalization,
  component corrections, and multiplication formulas. The normalization
  without the discriminant constant is explicitly fixed in Section 4 above.
- The existing exact rank-one generator and genus-three de Rham certificates
  supply only the previously stated group decomposition and splitting data.

The machine checks certify algebra, valuations, finite component tests,
series coefficients, and modular arithmetic. The cited local-height and
Coleman–Gross theorems identify these computations with the asserted
heights; they are standard mathematical inputs, not claims proved by the
Python program.


### 5.18 Mixed genus-two divisors and exhaustive good-prime height terms

This section removes self-intersection ambiguity from the remaining
genus-two height calculation and computes all its finite contributions
away from the model's bad primes and 13. It does not compute the complete
height matrix. The exact certificate is
`scripts/check_genus_two_mixed_heights.py` (Python and SymPy).

Use the monic models and classes P,Q of 5.13 and 5.16:

    f_k=x(x-k)(x+5k)(x^2-68kx-275k^2), k=-15,-3.

Write A_P,A_Q for the effective supports of their Mumford representatives.
Introduce R=P+sQ, where s=1 for k=-15 and s=-1 for k=-3, with effective
support A_R given by

    k=-15: u_R=x^2-276x+11025, v_R=406350-7794x;
    k=-3:  u_R=x^2+65x+528,   v_R=-1076x-12804.

These are exact Jacobian identities, not merely compatible points on the
curve. To verify them, put a=x(P), b=y(P), and interpolate

    V=s*v_Q+(b-s*v_Q(a))*u_Q/u_Q(a).

Then (f_k-V^2)/((x-a)u_Q) has monic normalization u_R, and
-V modulo u_R equals v_R. The checker verifies every division and curve
identity. Each u is squarefree, gcd(u,v)=1, and the three u polynomials are
pairwise coprime over Q. Thus all three supports and their hyperelliptic
conjugates are disjoint where required for the mixed pairings below.

#### Recovering the global height block without self-pairings

Let iota(x,y)=(x,-y), and D_A=A-iota(A). The divisor of u_A is
A+iota(A)-2deg(A)infinity, so D_A represents twice the Jacobian class
[A-deg(A)infinity]. Write

    U=h(D_P,D_Q), Vh=h(D_P,D_R), W=h(D_Q,D_R).

Bilinearity of the global Coleman-Gross pairing gives exactly

    h(P,Q)=U/4,
    h(P,P)=(Vh-sU)/4,
    h(Q,Q)=s(W-U)/4.

Division by four loses no 13-adic precision. This reduces the global
matrix to three pairings of disjoint divisors and avoids choosing tangent
vectors for self-intersections. The formulas apply to the global height;
individual local summands depend on divisor representatives. In
particular the finite sums below must not be mislabeled as intrinsic
local heights of the original P,Q generators.

For each Mumford pair (u,v), an explicit differential with residue divisor
D_A is

    omega_A = n_A(x) dx/(u(x)y), n_A=v*u' mod u.

Its residues at A and iota(A) are respectively +1 and -1. The roots of u
are simple and non-Weierstrass, and deg(n_A)<deg(u) guarantees regularity
at the unique infinity. The numerators for P,Q,R are

    k=-15:
      320625/32,
      10732500x+611752500,
      59705100-1338444x;

    k=-3:
      10692,
      4266108-18468x,
      44332x+303996.

These are residue-correct differentials, not yet Coleman-Gross normalized
ones. Their holomorphic corrections must be determined for the splitting
compatible with the chosen Richelot transport and genus-three quotient.
No third-kind integration or such transport normalization is asserted here.

#### Every good-prime contribution away from 13

Let S={2,3,5,11,19,53}. Exact discriminants show that both monic models have
good reduction outside S. All displayed coefficients are integral outside
2; the effective supports are finite on these models, so they acquire no
intersection at infinity at a prime outside S.

The exact resultants in the order (P,Q),(P,R),(Q,R) are

    k=-15: 2268225/16, -21375/16, 48556065600;
    k=-3:  6318, 3762, 8301366.

Their numerator factorizations are

    k=-15:
      3^2*5^2*17*593,
      -3^2*5^3*19,
      2^6*3^6*5^2*7*19*313;

    k=-3:
      2*3^5*13,
      2*3^2*11*19,
      2*3^5*19*29*31.

For a good prime ell!=13 an intersection is possible only if ell divides
the respective resultant. Every resulting good-prime valuation is one.
The gcd of the two reduced u polynomials is linear; its ordinate is
nonzero, so x is a local parameter and the intersection multiplicity is
one. If the two ordinates agree, the signed intersection of D_A,D_B is
+2; if they are opposite it is -2. Using the height character
chi_ell(ell)=-log13(ell) from 5.17 gives the exhaustive sums over
ell outside S union {13}:

    k=-15:
      (D_P,D_Q): -2log13(17)-2log13(593);
      (D_P,D_R): 0;
      (D_Q,D_R): +2log13(7)+2log13(313).

    k=-3:
      (D_P,D_Q): 0;
      (D_P,D_R): 0;
      (D_Q,D_R): -2log13(29)-2log13(31).

For independently reproducible sign checks, the common x-coordinate and
the two ordinates modulo ell are

    k=-15, P/Q:
      ell=17:  (1,6,6);
      ell=593: (501,513,513);
    k=-15, Q/R:
      ell=7:   (3,2,5);
      ell=313: (311,299,14);
    k=-3, Q/R:
      ell=29:  (28,17,17);
      ell=31:  (11,5,5).

These are complete contributions on the indicated good-prime set, not
finite-search evidence. They show concretely why computing only bad-prime
contributions is insufficient for these chosen height representatives.

#### Semistable vertical corrections at 11,19,53

At these three primes each reduced curve is geometrically irreducible
with precisely one ordinary node. The node x-coordinate is respectively
0, k, 34k modulo the prime. The discriminant valuations are 2,2,1, giving
the node thicknesses 2,2,1. This follows locally by separating the unique
double-root factor: all other roots stay distinct, and the two-root
discriminant valuation equals the smoothing-parameter valuation.

The tangent coefficient after division by the square of the node factor
is respectively (9,2,31) for k=-15 and (5,14,10) for k=-3. Each is nonzero.
The assertion is geometric; split tangent directions are not assumed.

Resolving a thickness-two node introduces one exceptional component;
the dual graph has two vertices joined by two edges. The hyperelliptic
involution preserves both components, although it can interchange
branches. Consequently the horizontal divisor D_A=A-iota(A) has degree
zero on each component. Its vertical correction can be taken to be zero.
At thickness one the total model is already regular and its special fiber
has only one component, so again there is no vertical correction. This
argument applies to all three D_A, including ramified closed-point
supports, by functoriality of intersection with the involution-invariant
components.

This does NOT prove that their local height contributions vanish.
Horizontal intersections must still be calculated on the regular model.
In particular A_Q for k=-15 is ramified at 53: its u discriminant has
53-adic valuation one, so treating that support as two unramified sections
would be invalid.

#### Exact remaining inputs

The three normalized 13-adic mixed integrals and the horizontal terms at
2,3,5,11,19,53 remain uncomputed; regular-model corrections at 2,3,5 also
remain. At k=-3, P and Q share a reduction disc at 13, as witnessed by
their resultant's 13-factor. Disjoint support over Q does not permit a
disjoint-residue-disc integration shortcut.

The global reconstruction and the good-prime arithmetic have been
independently checked. The intersection interpretation uses the standard
finite-prime definition of the Coleman-Gross pairing on a regular model.
Neither a complete genus-two height block, the genus-three QC function,
its roots, nor a final rational-point sieve has been produced.



### 5.19 Complete mixed finite terms outside 2,3,5,13

For the fixed divisors of 5.18, this section finishes the horizontal
intersection calculation at 11,19,53. Together with 5.18 it computes every
finite local contribution outside {2,3,5,13}. The remaining three finite
primes and normalized 13-adic integrals are not computed.

The exact checker is `scripts/check_genus_two_semistable_heights.py`.
It uses the adjacent mixed-height certificate and existing SymPy.
The calculation was independently checked using explicit blowup charts
and quotient-module intersection lengths.

#### Regular-model intersections, not just vertical corrections

The vertical corrections for D_A=A-iota(A) vanish at these primes by
5.18. Nevertheless their horizontal intersections must be evaluated on
the regular model. All common reduced supports are as follows, with
tuples (x,y_A,y_B):

    k=-15, ell=19:
      P/R: (4,0,0), a node;
      Q/R: (6,8,11), a smooth point with opposite ordinates.
    k=-3, ell=11:
      P/R: (0,0,0), a node.
    k=-3, ell=19:
      P/R: (14,14,5), a smooth point with opposite ordinates;
      Q/R: (16,0,0), a node.

There are no other common reductions at 11,19,53. Completeness is checked
by polynomial gcds of each pair of reduced u polynomials; every nonconstant
gcd has degree one. In particular no unlisted common point over a finite
extension is being omitted.

The two smooth intersections have nonzero ordinates. Each corresponding
resultant has 19-adic valuation one, so x is a local parameter and the
intersection multiplicity is exactly one. The opposite ordinates give
D_A dot D_B=-2. With chi_19(19)=-log13(19), each contributes
+2log13(19).

At a shared node of thickness two, use the blowup chart

    X=(x-a)/ell, Y=y/ell,

where a=0 at 11 and a=k at 19. Every support passing through a shared
node has a simple u-root modulo ell. Hensel lifting modulo ell^2 therefore
determines its unramified section and its reduction on this chart.
The exact coordinates (X,Y) modulo ell are

    k=-15, ell=19, a=-15:
      P: (18,9), R: (6,4);
    k=-3, ell=11, a=0:
      P: (3,4), R: (4,10);
    k=-3, ell=19, a=-3:
      Q: (16,13), R: (10,3).

These sections meet the affine part of the exceptional component. In
each potentially intersecting pair the X-coordinates differ. Consequently
their strict transforms are disjoint, also after either hyperelliptic
conjugation. All node pairs thus contribute zero on the regular model.

At 53, Q for k=-15 has ramified support at the thickness-one node. Neither
P nor R shares any reduced support with Q there, so there is no horizontal
intersection. No unramified-section assumption is made for that support.

#### Independent module-length check

Before resolving a node, the same-sheet and opposite-sheet intersection
lengths can be checked as the lengths of

    Z_ell[x]/(u_A,u_B,v_A-v_B),
    Z_ell[x]/(u_A,u_B,v_A+v_B).

In the free module Z_ell[x]/u_A, multiplication by u_B and by v_A minus
or plus v_B gives a rectangular presentation matrix. The minimum valuation
of its maximal minors is its cokernel length. All denominators in this
calculation are powers of two, hence units at the three primes used.

The resulting same/opposite triples, in prime order (11,19,53), are

    k=-15:
      P/Q: (0,0,0), (0,0,0);
      P/R: (0,1,0), (0,1,0);
      Q/R: (0,0,0), (0,1,0).
    k=-3:
      P/Q: (0,0,0), (0,0,0);
      P/R: (1,0,0), (1,1,0);
      Q/R: (0,1,0), (0,1,0).

These are singular-model lengths where a node occurs, not automatically
regular-model intersection numbers. The explicit blowup calculation above
removes those shared node intersections and proves the claimed heights.
In each node pair the same/opposite corrections cancel in the signed
pairing, agreeing with the independent module calculation.

#### The complete certified finite sum outside 2,3,5,13

At 11 and 53 all three mixed local pairings are zero for both curves.
At 19 the only nonzero pairings are Q/R for k=-15 and P/R for k=-3,
each equal to +2log13(19). Adding the good-prime terms from 5.18 gives

    k=-15:
      (D_P,D_Q): -2log13(17)-2log13(593);
      (D_P,D_R): 0;
      (D_Q,D_R): +2log13(7)+2log13(19)+2log13(313).
    k=-3:
      (D_P,D_Q): 0;
      (D_P,D_R): +2log13(19);
      (D_Q,D_R): -2log13(29)-2log13(31).

These are exact identities for the sum over all finite primes other than
2,3,5,13, with the fixed representatives of 5.18. They are neither complete
global heights nor the full away-from-13 sums. In particular, no vanishing
assumption is made at 2,3,5.

Still needed are regular-model horizontal and vertical terms at those
three primes, normalized third-kind integration at 13 with certified
precision, and the compatible splitting under Richelot transport.
Only then can these mixed pairings recover the global genus-two block.
The genus-three QC function, root isolation and final completeness sieve
remain uncomputed; the original size-eleven problem is still open here.



### 5.20 Odd bad-prime terms: all finite contributions outside 2 and 13

This section completes the fixed mixed divisors' local heights at 3 and 5.
Only k=-3, P/Q contributes at 3, with value -log13(3); only k=-15, P/R
contributes at 5, with value -2log13(5). All other mixed terms there vanish.
Together with 5.18–5.19 this computes every finite contribution outside
{2,13}, but not the full away-from-13 or global height.

The certificate `scripts/check_genus_two_odd_local_heights.py` uses only
existing SymPy. The charts, support positions and intersection lengths
were independently derived and checked. No software was installed.

#### Character normalization and regular neighborhoods

Put z=x/k, w=y/k^2 and

    F(z)=z(z-1)(z+5)(z^2-68z-275), w^2=kF(z).

When needed at ell=3 or 5, use K=Q_ell(pi), pi^2=ell, and extend the
character by chi_K=chi_ell composed with the norm. Then
chi_K(pi)=-log13(ell), since the norm of pi is -ell and log13(-1)=0.
The local height after base change is twice the original local height;
divide the final intersection coefficient by two.

The displayed charts arise on the normalization of blowups of the
z-line separating the indicated branch discs. At every support point
used, their displayed equations are smooth over the valuation ring.
A proper regular model dominating these charts may be chosen; further
resolution away from these smooth neighborhoods does not change their
intersections or components. We need not assert an entire minimal
dual graph.

For every D_A=A-iota(A), the calculations below verify degree zero on
each component of this regular model. Consequently its vertical
correction may be taken to be zero. This does not remove the need to
calculate horizontal intersections.

#### Prime 3: charts

Write k=3u, u=-5 or -1, and Y=w/pi. Thus Y^2=uF(z).
The cluster chart z=1+3t, Y=9Yc is

    Yc^2=u(1+3t)t(t+2)(t^2-22t-38),

with reduction Yc^2=u t(t-1)(t-2)^2. The outer normalization and this
cluster normalization are geometrically irreducible and preserved by
iota. The intermediate chart z=1+pi*r, Y=pi^2*j is

    j^2=u(1+pi*r)r(r+2pi)(r^2-22pi*r-114),

with reduction j^2=u r^4. Finally the twin chart t=11+pi*X, Yc=pi*V is

    V^2=u(1+3(11+pi*X))(11+pi*X)(13+pi*X)(X^2-53).

Its reduction is V^2=2u(X^2+1), a smooth, geometrically irreducible conic
preserved by iota. Every identity follows by substitution using pi^2=3.

For k=-15, the supports are:
- P meets the smooth outer component at z=0.
- With a^2=265, the two Q points have
  r=11pi plus or minus a, j=2650pi plus or minus 302a.
  Since a^2=1 modulo 3, these meet the two intermediate components
  j=plus or minus r^2, one point on each. The involution exchanges the
  components, but A_Q and iota(A_Q) each have degree one on each.
- With b^2=11, R has twin-chart coordinates
  X=pi(-24 plus or minus b)/5,
  V=(-2754 plus or minus 866b)/25.
  Thus X=0 and V^2=2 modulo pi. Both are smooth points of the same
  invariant component; their unramified quadratic support is retained.

The three effective supports occupy different components and avoid their
intersections. All anti-divisor component degrees are zero, including
the exchanged intermediate components. Hence all three local mixed
heights at 3 vanish for k=-15.

#### Prime 3: the nonzero intersection for k=-3

Let s^2=-29, s=1 modulo 3. The complete support positions are:
- P: t=-4, X=-5pi, V=44 in the twin chart.
- One Q point: t=3-s, X=(-8-s)/pi, V=-38(1+s), also in that chart.
- The other Q point: t=3+s, a smooth cluster point at t=1,Yc=0.
- One R point is at smooth outer infinity. The other is a smooth cluster
  point with t=1 modulo 3 and 9t^2-59t+38=0.

Every component met here is preserved by the involution, so again every
D_A has zero degree on every component.

P and the twin Q point reduce to (X,V)=(0,2). Their difference is
X_Q-X_P=(7-s)/pi. Since (7-s)(7+s)=78 and 7+s is a unit,

    v_K(X_Q-X_P)=1.

The ordinate V is nonzero, so X is a local parameter. The same-sheet
intersection has length one; the opposite-sheet intersection is absent.
The signed intersection of D_P,D_Q is two over K. Norm-character
normalization and division by two therefore give

    h_3(D_P,D_Q)=-log13(3).

For the cluster Q/R pair, the ordinates before scaling are

    w_Q=1026(s-1), w_R=4(269t_R-266).

Their 3-adic valuations are 4 and 3. Because Yc=w/(9pi), their K-valuations
are 3 and 1. At this smooth branch point Yc is a local parameter, so both
same-sheet and opposite-sheet differences have valuation one. The signed
contributions cancel. P/R occupy different components and do not intersect.
This proves that the other two mixed local heights vanish.

These valuations are exact. The certificate Hensel-lifts the simple roots
modulo 3^12 and tests valuations strictly below that precision; the omitted
root tails cannot change them.

#### Prime 5

For k=-3, the original monic model reduces to

    y^2=x^3(x+3)(x+4) modulo 5.

Its only singular affine point is x=y=0. None of the three u polynomials
vanishes at x=0 modulo 5, and every pairwise resultant is a 5-adic unit.
Thus the supports are disjoint smooth points on one invariant component.
All three local mixed heights vanish, without resolving remote singularities.

For k=-15, use pi^2=5 and eta=w/pi, so eta^2=-3F(z).
The two Q points have

    z=34 plus or minus pi*sqrt(159),
    eta=4770pi plus or minus 906sqrt(159),

and reduce to z=4, eta=plus or minus 2. They are smooth outer points.
The outer component is geometrically irreducible and invariant; after
removing its repeated factor its equation is
(eta/z)^2=-3z(z-1)(z+2) modulo 5.

R has one point at the unique smooth outer infinity. In coordinates
v=1/z, h=eta/z^3 the equation is

    h^2=-3v(1-v)(1+5v)(1-68v-275v^2),

whose derivative with respect to v is a unit at v=h=0.

The inner chart t=z/5, xi=eta/pi^3 has equation

    xi^2=-3G(t),
    G=t(5t-1)(t+1)(5t^2-68t-55).

The identity F(5t)=125G(t) proves the transformation. Its reduction is
xi^2=t^2(t+1), a geometrically irreducible invariant component. P and the
inner R point meet its smooth point t=3,xi=1, not its node at t=0.
These are all support points. Thus every anti-divisor has zero degree
on every component, including any components added by further resolution.

For the actual intersection put W=9xi=y/625. Then

    (t_P,W_P)=(-3/4,513/32),
    25t_R^2+92t_R+49=0, t_R=3 modulo 5,
    W_R=(23382t_R+16254)/25.

Hensel lifting gives t_R=103 modulo 125 and W_R=4 modulo 5, equal to W_P.
Moreover evaluating the R polynomial at t_P gives -95/16, of valuation
one; its derivative is a unit. Hence v_K(t_R-t_P)=2.
The nonzero ordinate makes t a local parameter. The same-sheet intersection
has length two over K; the opposite-sheet intersection is absent.
The signed intersection is four, and descent gives

    h_5(D_P,D_R)=-2log13(5).

All other pairs occupy disjoint neighborhoods and have zero contribution.
The checker uses precision 5^8, beyond all divisions and valuations needed
for these assertions.

#### Complete certified finite sums outside {2,13}

For the fixed divisor representatives of 5.18 the sums are now

    k=-15:
      (D_P,D_Q): -2log13(17)-2log13(593);
      (D_P,D_R): -2log13(5);
      (D_Q,D_R): 2log13(7)+2log13(19)+2log13(313).
    k=-3:
      (D_P,D_Q): -log13(3);
      (D_P,D_R): 2log13(19);
      (D_Q,D_R): -2log13(29)-2log13(31).

They are representative-dependent partial sums, not complete global
heights. The dyadic terms, normalized 13-adic integration and compatible
splitting transport remain uncomputed. Consequently the full genus-two
height block, genus-three QC function, root isolation and final
completeness sieve still remain unfinished. No original PTE solution or
general exclusion is claimed.


### 5.21 Dyadic terms and complete away-from-13 height inputs

All six dyadic mixed local heights for the fixed divisors of 5.18 are zero.
This completes every non-13 finite contribution needed for those three
mixed pairings on each genus-two curve. It does not compute their normalized
13-adic integrals, and hence does not complete the global height block.

The certificate is `scripts/check_genus_two_dyadic_heights.py`.
The regular coordinates, component degrees, exact quadratic-field
valuations and signed intersection cancellation were independently checked.
Only existing SymPy is used.

#### k=-3: smooth neighborhoods suffice

For a=1 and a=3, set x=a+4t,y=4Y. The polynomial f_-3(a+4t)/16 has integral
coefficients and reduces to t+1 modulo 2. Therefore the resulting affine
model Y^2=f_-3(a+4t)/16 is smooth everywhere on its special fiber: the
derivative with respect to t is a unit. Its geometrically irreducible
component is preserved by iota, which sends Y to -Y.

P=(33,10692) lies in the a=1 chart. The odd root of u_R lies in the a=3
chart: u_R(1)=594 has valuation one and u_R'(1)=67 is odd. Its ordinate
is divisible by four. The other R root lies at the smooth point x=0 of
the original integral model. The polynomial u_R reduces to x(x+1), so
these are unramified local sections.

For Q, writing x=1+d gives a quadratic whose constant and linear
coefficients have valuations one and one. Its Newton polygon gives
v_2(d)=1/2 for both roots. Thus neither Q point lies in any of the above
smooth neighborhoods. P and R occupy distinct neighborhoods as well.

These charts may be included in one proper regular model by dominating
the corresponding blowups and resolving elsewhere. D_P and D_R have
degree zero on every component, since each point and its hyperelliptic
conjugate meet the same invariant component. Their vertical corrections
can therefore be taken to be zero. All their horizontal intersections
with the other divisors are empty. Using symmetry for the Q/R pairing
proves that all three local heights vanish. We do not need to determine
Q's own vertical correction.

#### k=-15: the two pairings involving P

Here v_2(x(P))=-2 and v_2(y(P))=-5. In the infinity chart v=1/x,h=y/x^3,
the equation is

    h^2=v+960v^2-124200v^3+2565000v^4+69609375v^5.

Its derivative with respect to v at infinity is a unit, so this is a
smooth neighborhood. Both P and iota(P) meet that same invariant
component. The monic u_Q,u_R have integral odd roots and reduce to
(x+1)^2; their supports remain finite and cannot meet this infinity
neighborhood on a dominating regular model.

Thus D_P has zero component degree, its vertical correction is zero,
and the P/Q and P/R horizontal intersections are empty. Both local
heights vanish.

#### k=-15: the remaining Q/R pairing on a wild regular chart

Let L=Q_2(r), r^2=11, and pi=1+r. This is a ramified quadratic extension
with uniformizer pi. Normalize v_L(pi)=1, so v_L(2)=2.
The field also contains the Q support: 795/11 is 1 modulo 8 and is a
square in Q_2. Choose c with c^2=795/11 and c=1 modulo 8. This last
choice is possible because the ratio is 1 modulo 16.

Write

    x_Q=-510 plus or minus 15cr,
    x_R=138 plus or minus 27r,
    a=138+27r.

Set X=(x-a)/pi^2 and Y=y/pi^3. Exact arithmetic in Q[r]/(r^2-11) gives
the valuations of the Taylor coefficients f_-15^(j)(a)/j!, for j=0,...,5:

    [6,4,2,2,0,0].

After multiplication by pi^(2j-6), these are

    [0,0,0,2,2,4].

The first three coefficients reduce to one. Consequently the integral
model

    Y^2=f_-15(a+pi^2 X)/pi^6

reduces to Y^2=1+X+X^2. The derivative with respect to X is one.
Thus the entire affine special fiber is smooth and geometrically
irreducible. Its component is preserved by iota.

All Q and R support points and their hyperelliptic conjugates lie in
this chart. For R this follows directly from its displayed roots and
v_L(y_R)=3. For Q, its x-distances from a have valuations four and two;
also v_L(y_Q)=3, checked by the exact quadratic-field norm. Thus their
X,Y coordinates are integral. Every D_Q,D_R component degree is zero,
so vertical corrections vanish.

For completeness the intersection lengths can be obtained without an
approximate quadratic root. Matched signs of Q/R give

    x_Q-x_R=-648 plus or minus (15c-27)r,

of 2-adic valuation two; mismatched signs have valuation one because
15c+27 has valuation one. Mismatched points therefore have distinct
reductions in X and cannot intersect.

For a matched pair put q=x_Q and d=x_R-x_Q, so v_2(q)=0, v_2(d)=2.
The ordinate formulas give

    y_R-y_Q=5796q-7794d+1971000,
    y_R+y_Q=-21384q-7794d-1158300.

In the first line the three term valuations are 2,3,3; in the second
they are 3,3,2. Thus both differences have valuation exactly two, or
v_L=4. Both points reduce to the same smooth point of the displayed
chart. Since the derivative with respect to X is a unit, Y is a local
parameter. Its same-sheet and opposite-sheet differences both have
valuation 4-3=1. The signed intersections cancel, for each matched
pair. Hence the local height over L is zero; local-height base change
gives zero over Q_2 as well.

This uses a regular model and signed intersections, not a claim that
naive singular-model resultant valuations give local heights.

#### Complete away-from-13 input

Let U_away,V_away,W_away denote the sums of the finite local mixed
pairings for (D_P,D_Q),(D_P,D_R),(D_Q,D_R), respectively, over ALL
primes other than 13. Combining 5.18–5.21 gives the exact formulas

    k=-15:
      U_away=-2log13(17)-2log13(593),
      V_away=-2log13(5),
      W_away=2log13(7)+2log13(19)+2log13(313).

    k=-3:
      U_away=-log13(3),
      V_away=2log13(19),
      W_away=-2log13(29)-2log13(31).

With the same branch log13(13)=0 as before, the certified residues
modulo 13^10 are

    k=-15: (76380091502, 69520774258, 79405079481);
    k=-3:  (24898284736, 34648651414, 71972655326).

For each unit q in these formulas the checker uses
log13(q)=log(q^12)/12. The log(1+z) terms have valuations at least
j-v_13(j); terms after j=12 have valuation at least 12. Calculating z
modulo 13^12 and using exact rational sums therefore certifies the
residues modulo 13^10. A separate computation at precision 13^20
reduces to the same values. No floating-point approximation is used.

If U_13,V_13,W_13 are the still-uncomputed normalized 13-adic mixed
integrals, the global reconstruction of 5.18 applies to
U=U_13+U_away, V=V_13+V_away, W=W_13+W_away. The away sums depend on
these divisor representatives and must not be substituted for intrinsic
local generator heights under a different representative choice.

All finite-prime inputs away from 13 for these pairings are now complete.
The missing arithmetic is normalized 13-adic third-kind integration with
the compatible splitting and Richelot transport. The genus-two global
height block, full genus-three QC function, certified roots and final
completeness sieve remain unfinished. The original distinct-entry
size-eleven problem is not solved.


### 5.22 Explicit Richelot de Rham transport and the normalized splitting

The previously missing splitting transport is now explicit. This is an
exact algebraic calculation, not a normalized Coleman integral.
The certificate is `scripts/check_richelot_de_rham_transport.py` (SymPy).
An independent derivation checked the correspondence, trace, exact
primitives and cup products.

Put A=x^2+5, B=4x^2+30x+80, C=x^2+275, delta=1 or 5, and k=-15/delta.
The curves are delta*y^2=ABC and Y^2=kF(z), with F from 5.20.
A Richelot correspondence realizing the isogeny phi:J(C_delta)->J(D_k)
of 5.12 is, on its ordinary affine open,

    E=A(x)(z^2-68z-275)+18B(x)z=0,
    yY=kA(x)(z^2-68z-275)(x-z).

The bracket identity of 5.12 is
sum G_i(x)h_i(z)=-(x-z)^2, with h_i=L_i/Delta.
It verifies the squared ordinate equation on E=0. We use the displayed
ordinate sign. Changing the sign changes the isogeny by -1, not the
transported complementary subspace. The correspondence extends to the
Jacobian isogeny supplied by the Richelot theorem cited in 5.12.

#### Full pullback matrix, including exact corrections

On D put beta_j=z^j dz/Y, j=0,...,3. On C use the basis c0,c1,c2,c3 of
5.15. The rows of the following matrix express phi*(beta_j) modulo exact
differentials in that ordered C basis:

    T=[[2,0,0,0],
       [0,2,0,0],
       [0,-32/3,-8/15,0],
       [-2200/9,-2464/9,-1114/45,-16/45]].

This is not only a holomorphic tangent-space calculation. To derive it,
differentiate E, so dz/dx=-E_x/E_z, substitute the ordinate relation, and
take the trace in the quadratic algebra Q(x)[z]/E. For the first two
forms the trace coefficients multiplying dx/y are 2 and 2x. For the
last two they are

    -2(4x^3+265x^2+1100x-1375)/(x^2+5),

    2(291x^5+3220x^4+154650x^3+880000x^2
      +474375x-1512500)/(x^2+5)^2.

Subtracting the corresponding T rows leaves the exact differentials

    d(2delta*y/(15A)),
    d(delta*y*H/A^2),

where H=-176/9-214x/9+256x^2/45+2x^3/45.
For example, d(delta*y*R) equals
(R'ABC+R(ABC)'/2)dx/y. Thus both reductions are polynomial/rational
identities, checked without period or point-search assumptions.

A separate Laurent calculation at D's unique infinity gives its cup matrix

    CD=(1/k)*[[0,0,0,4/3],
              [0,0,4,512/3],
              [0,-4,0,-736],
              [-4/3,-512/3,736,0]].

With CC from 5.15 one has T*CC*T^transpose=2*CD and det(T)=512/675.
This verifies the polarization multiplier of the (2,2)-isogeny and
the chosen cup/sign convention.

#### Transported complement and explicit normalization projection

The chosen C complement is span(c2,c3+525c1); adding 45c2/8 to its
second generator gives exactly the equivalent complement of 5.15.
Its inverse image under phi* is

    W_D=span(beta2+16beta1/3,
             beta3+392beta1/9+1100beta0/9).

Direct multiplication verifies both that its image is the chosen C
complement and that it is isotropic for CD.

For the monic model X=kz, v=k^2Y, put eta_j=X^j dX/(2v).
Since eta_j=k^(j-1)beta_j/2, this becomes

    k=-15: span(eta2-80eta1, eta3-412500eta0+9800eta1);
    k=-3:  span(eta2-16eta1, eta3-3300eta0+392eta1).

The full monic transport matrix is diag(k^(j-1)/2)*T.
Its determinant is 32k^2/675, a 13-adic unit. Its entries and those of
its inverse are 13-integral; this basis transport loses no 13-adic digits.

In particular, once the Coleman-Gross de Rham class of a third-kind
form has actually been computed as sum a_j eta_j, the holomorphic
differential to subtract is explicitly

    (a0-(1100k^3/9)a3) eta0
      +(a1-(16k/3)a2-(392k^2/9)a3) eta1.

This formula fixes the normalization target. It does not compute the
a_j or the integral of the normalized form; residue-correct forms alone
do not supply these missing data.

#### Polarization-compatible height weighting on the genus-three cover

Let psi be the polarized dual Richelot isogeny, so psi*phi=[2], and
g:X_delta->C_delta the degree-two quotient. Functoriality of heights
with these compatible splittings gives

    h_C(psi D,psi E)=2h_D(D,E),
    h_X(g*psi D,g*psi E)=4h_D(D,E).

The correspondence Z acts by -2 on this factor. Hence

    h_X(g*psi D,Zg*psi E)=-8h_D(D,E).

This supplies the genus-two block's weighting, complementary to the
elliptic block's +8 weighting in 5.17. It does not assert explicit
transported divisor coordinates or a completed global height matrix.

Splitting transport is therefore no longer an unspecified input.
Normalized third-kind integrals, their global-index data, explicit
divisor/endpoint computations and the full QC function remain missing.


### 5.23 Certified Frobenius matrices at 13

This section computes actual de Rham Frobenius matrices, not only their
characteristic polynomials. The exact-rational checker is
`scripts/check_genus_two_frobenius13.py` (existing SymPy).
Independent modular and f-adic-first implementations gave the same
matrices, and an independent exact-rational computation checked every
entry.

In the monic basis eta_i=x^i dx/(2y), i=0,...,3, columns are Frobenius
images. The matrices, each with entrywise error O(13^3), are

    k=-15:
      [[949,39,842,458],
       [481,1456,1678,583],
       [234,1248,1524,661],
       [1092,1391,833,465]];

    k=-3:
      [[1248,871,933,2000],
       [1989,741,1422,1207],
       [741,351,673,1186],
       [1911,377,229,1732]].

#### Exact finite reduction algorithm

Let f=f_k, p=13, and R=f(x^13)-f(x)^13. The polynomial R is divisible
by 13 and has degree at most 64. Use the usual Frobenius lift
x->x^13, y->y^13(1+R/y^26)^(1/2). The pullback of eta_i is the sum over
n>=0 of

    13 binom(-1/2,n) x^(13(i+1)-1) R^n dx/(2y^(13+26n)).

To reduce N(x)dx/(2y^(2m+1)), compute B=N*(f')^(-1) modulo f, with
degree B<5, and A=(N-Bf')/f. Then replace the numerator and exponent by

    N <- A+2B'/(2m-1), m <- m-1.

This identity follows by subtracting an exact differential. Once m=0,
a polynomial of degree e>=4 is reduced using d(x^(e-4)y), whose leading
coefficient in the dx/y convention is (2e-3)/2. Each operation in the
checker is over Q, so no modular working-precision or rounding issue
is hidden in these divisions.

#### Certified infinite-tail bound

Put R=13E. The binomial coefficient binom(-1/2,n) is 13-integral, so
the n-th term has coefficient valuation at least n+1 before reduction.
For bounding denominators, expand its numerator f-adically; since f
is monic, each block coefficient remains integral and has degree at
most four.

For a block with a negative odd power of y, apply Kedlaya's denominator
bound: a degree-at-most-2g numerator divided by y^(2m+1) needs at most
floor(log_p(2m+1)) powers of p to clear reduction denominators.
Here that loss is at most floor(log_13(13+26n)).
This is the explicit standard theorem input:
[Kedlaya, Counting Points on Hyperelliptic Curves using
Monsky-Washnitzer Cohomology, Section 3, Lemma 2](https://arxiv.org/html/math/0105031).

The polynomial dx/(2y) part has degree at most

    51+64n-5(6+13n)=21-n.

We bound its horizontal denominators directly. For e=4,...,21 the
only factors 2e-3 divisible by 13 are 13 and 39. Each degree is
eliminated at most once, so at most two 13-adic digits can be lost.
No positive-power denominator lemma is needed.

The reduced n-th term therefore has valuation at least

    n+1-max(floor(log_13(13+26n)),2).

For every n>=4, both quantities inside the maximum are at most n-2:
13+26n<13^(n-1), starting at n=4 and preserved as n increases.
Thus every omitted term has valuation at least three, and the bound
tends to infinity. The entire omitted tail is O(13^3).
Exact reduction of n=0,1,2,3 is therefore sufficient to certify the
displayed matrices.

#### Independent consistency and precision checks

Both matrices have characteristic polynomial
T^4+6T^2+169 modulo 13^3, agreeing with the independent finite-field
counts. This is a consistency check, not the proof of truncation.

The monic cup matrix has nonzero upper-triangular entries

    C03=1/3, C12=1, C13=128k/3, C23=-184k^2.

It follows either by direct Laurent residues at infinity or from 5.22.
The checker verifies M^transpose*C*M=13C modulo 13^3.

There is also an independent twist check. The coordinate identification
over Q_13(sqrt(5)) scales the differential basis by
D=diag(1/5,1,5,25), apart from the common sqrt(5) factor. Since 5 is
a nonsquare modulo 13, Frobenius changes that factor's sign. Accordingly

    D*M_-15+M_-3*D=0 modulo 13^3,

which the matrices satisfy.

The operator M-I has determinant 176, a 13-adic unit, so its inverse
can be formed without 13-adic precision loss. In contrast,
det(M-13I)=13^2*176. Its reduction has rank two, and its Smith
valuations are 0,0,1,1. Thus inverse entries may have valuation -1.
If the matrix itself is known only to O(13^3), a general inverse is
guaranteed only to O(13): the perturbation contains an inverse on
both sides. We do not silently budget just one digit for that
approximate-matrix inversion.

#### Remaining integral computation

These are cohomology matrices. To compute Coleman integrals one must
also retain and evaluate the exact primitives removed during reduction,
perform the necessary tiny integrals on the actual endpoints, and handle
the quadratic supports and shared reduction discs correctly. The
third-kind de Rham classes and their normalization then require the
compatible data of 5.22.

No normalized mixed third-kind integral, global genus-two height block,
complete genus-three QC function, root isolation or final point sieve
has yet been computed. The original distinct-entry size-eleven problem
remains unsolved.


### 5.24 Actual Coleman basis integrals and third-kind normalization

This section computes the basis integrals needed to normalize all six
third-kind forms of 5.18. It is no longer only a cohomology-matrix
calculation. The certificate `scripts/check_coleman_basis13.py` uses
exact rational reductions and finite unramified etale arithmetic, with
existing SymPy only. Two independent implementations agree.

Write I_A,j=integral over D_A of eta_j, eta_j=x^j dx/(2y).
The following vectors have entrywise error O(13^3):

    k=-15:
      P: (1768,533,2062,1190);
      Q: (195,117,182,1343);
      R: (1963,650,14,22).
    k=-3:
      P: (1079,533,852,402);
      Q: (598,572,1772,1289);
      R: (481,2158,1297,214).

#### Retaining the reduction primitives

The reductions in 5.23 now retain their odd exact primitives. In a vertical
step N=A*f+B*f', the primitive added is

    -B/((2m-1)y^(2m-1)).

Removing a polynomial leading term of degree e contributes

    lead(N)*x^(e-4)*y/(2e-3).

Differentiation verifies both signs and the factors of two. If g denotes
the column vector of accumulated primitives and M has columns equal to
Frobenius images, then

    phi*eta=M^transpose*eta+dg.

Oddness under the hyperelliptic involution fixes the otherwise irrelevant
primitive constants.

For a Mumford support (u,v), use the finite etale algebra
S=Z_13[alpha]/u(alpha), with b=v(alpha). Every u is squarefree modulo 13,
every coefficient is 13-integral, and b is a unit. This includes split
quadratic algebras as well as the unramified quadratic field.

Hensel-lift alpha^13 to the unique root sigma(alpha) of u in its residue
disc. Put a=sigma(alpha), b_sigma=v(a), h=alpha^13-a.
The tiny integral from sigma(A) to phi(A), through precision 13^3, is

    a^j*h/(2b_sigma)
      +(j*a^(j-1)/(2b_sigma)
        -a^j*f'(a)/(4b_sigma^3))*h^2/2.

The first derivative term with j=0 is zero. The ordinate branch is fixed
by b_sigma and its reduction, which agrees with b^13. Since h is divisible
by 13, every omitted integrated power has valuation at least three,
including possible denominators divisible by 13.

Taking traces sums the conjugate support points. Frobenius functoriality
and the anti-invariant divisors give the exact linear equation

    (Id-M^transpose) I_A
       =2 Tr_S(g(alpha,b))-2 Tr_S(tiny).

There is no division by the support degree. The two factors already
account for A-iota(A). Since det(Id-M)=176 is a 13-adic unit, solving
preserves the certified precision.

The primitive tail is also O(13^3), not just its cohomology class.
The finite-branch proof of the denominator estimate in 5.23 bounds the
primitive coefficients by the same power of 13; polynomial reductions
lose at most two digits. All endpoints have integral x and unit y.
Thus the omitted n>=4 primitives evaluate to O(13^3).
The checker additionally uses four guard digits for finite endpoint
evaluation: primitive coefficients are scaled by 13^4, checked integral,
evaluated modulo 13^7 and divided back. Every division is checked.

#### Checks beyond holomorphic logarithms

The first two coordinates agree with 5.16, including the factor two for
D_A. They obey I_R=I_P+sI_Q in their holomorphic coordinates.

The second-kind coordinates need not be additive on Jacobian classes.
An independent principal-divisor check uses

    k=-15: V=33x^2/2+3240x-448875/2, s=1;
    k=-3:  V=-10x^2+426x+7524, s=-1.

For q=(y-V)/(y+V), exact interpolation gives
div(q)=D_P+sD_Q-D_R. Laurent expansion at infinity gives

    Res(log(q)*eta)=(0,0,33,-88821/4), k=-15;
    Res(log(q)*eta)=(0,0,-20,12076/3), k=-3.

These reduce modulo 13^3 to I_P+sI_Q-I_R. Thus the computed second-kind
coordinates also pass an independent reciprocity check.

#### Sign-correct normalization

For omega_A=(v*u' mod u)dx/(u*y), let Psi(omega_A)=sum a_j eta_j.
Use the signed cup matrix C from 5.23. The global-symbol convention is
minus Res(omega*integral eta). Its relation to Psi is given in
[Gajovic–Muller, Definition 2.8, Proposition 2.9 and equation (6.1)](https://arxiv.org/html/2307.15787);
see also the [Balakrishnan–Besser sign errata](https://math.bu.edu/people/jbala/cg_heights_errata.pdf).

There is no infinity correction here. In t=x^(-1/2), omega_A=O(t^4)dt,
whereas a primitive of eta_j has order at least 3-2j apart from its
constant. Their product has no residue for j<=3. The symbol vector is
therefore -I_A. With our antisymmetric cup convention this yields

    a=C^(-1) I_A,

not its negative. Explicitly,

    a0=552k^2 I1+128k I2-3I3,
    a1=-552k^2 I0-I2,
    a2=-128k I0+I1,
    a3=3I0.

Applying the transported complement of 5.22, the holomorphic correction
to subtract from omega_A is lambda0*eta0+lambda1*eta1, where

    lambda0=-1100k^3 I0/3+552k^2 I1+128k I2-3I3,
    lambda1=-16k I1/3-I2.

All coefficients are 13-integral and det(C)=1/9, so this introduces no
additional precision loss. The certified correction pairs modulo 13^3 are

    k=-15: P (1545,1032), Q (1301,390), R (1238,1455);
    k=-3:  P (2021,1085), Q (1409,789), R (614,276).

The forms are now normalized to the specified precision. Integrating
these normalized third-kind forms against another divisor is a separate
step, undertaken for three pairings below.


### 5.25 Three normalized mixed heights and partial global height matrices

This section computes three actual mixed heights, not just normalization
data. The certificate `scripts/check_normalized_mixed_heights13.py` uses
the adjacent basis-integral checker and existing SymPy. Independent
implementations, using two coordinate descriptions of the endpoint
algebras, agree on every intermediate matrix and integral.

All congruences below have absolute error O(13^3), for the splitting and
logarithm branch already fixed. The normalized local results are

    k=-15:
      h_13(D_P,D_Q)=1755;
      h_13(D_P,D_R)=884.
    k=-3:
      h_13(D_P,D_R)=1781.

#### Raw third-kind integrals via explicit even models

For the rational point P=(a,b), set

    u=1/(x-a), v=y*u^3/b,
    G_P(u)=u^6*f(a+1/u)/b^2.

Since f(a)=b^2, this is a monic degree-six model. Direct differentiation
gives the exact pullback identity

    b dx/((x-a)y)=-u^2 du/v=-2 eta_even,2,

where eta_even,j=u^j du/(2v). Thus a raw omega_P integral is minus twice
the corresponding even-model basis integral. No further factor arises
from D_P already representing twice a Jacobian class.

Both even models have good reduction at 13. Their reductions are

    k=-15: G=u^6-u^5+4u^3+3u^2+u;
    k=-3:  G=u^6-u^4+4u^3-2u^2+4u.

Each is squarefree. For the three selected target divisors, x-a is a
13-adic unit, so all transformed endpoints have integral unit u,v.
In particular the shared-disc k=-3 P/Q pair is NOT included here.

For an effective degree-two support (u_B,v_B), the transformed monic
support polynomial is the monic normalization of
u^2*u_B(a+1/u); the ordinate is v_B(a+1/u)*u^3/b, reduced modulo it.
The checker verifies the transformed curve equations and etale support
conditions rather than assuming them.

#### Certified even-degree integration

The even affine Frobenius calculation uses five basis forms,
eta_even,0 through eta_even,4. It retains the third-kind direction
associated with the two infinities. Its characteristic polynomial is
(T-13)(T^4+6T^2+169), and Id-M is invertible modulo 13.
The checker also verifies that the infinity-residue row r satisfies
r*M=13r, so that the residue direction and column convention are tested.

Vertical reduction and its primitive are unchanged from 5.24.
For polynomial degree e>=5, the removed primitive is

    lead(N)*u^(e-5)*v/(2e-4).

Here is the required even-degree truncation justification. The numerator
difference G(u^13)-G(u)^13 has degree at most 77 and is divisible by 13.
An f-adic negative block has numerator degree at most five. At the six
distinct finite branch points, use y as uniformizer and integrate its
polar part. Clearing at most floor(log_13(2m+1)) digits makes those
primitive polar coefficients integral. The values of the leading
degree-at-most-five primitive numerator at all six roots then recover
that polynomial by a unit Vandermonde matrix. Subtract and descend
through the remaining powers. This proves the same bound for the
entire negative-block primitive and its reduced form. It does not
blindly apply the odd-degree lemma outside its hypotheses.

The polynomial part has degree at most

    64+77n-6(6+13n)=28-n.

Horizontal reduction divides by 2e-4 for e=5,...,28. Only 26 and 52
supply factors of 13, so at most two digits are lost. As in 5.23,
every n>=4 term and its primitive evaluate to O(13^3) at the selected
finite unit endpoints. Oddness again fixes primitive constants.
Exact reduction of n=0,...,3, followed by the etale trace/tiny-integral
method of 5.24, is therefore sufficient.

The computed even-basis integral vectors over the transformed D_B are

    k=-15, B=Q: (416,2002,1222,323,714);
    k=-15, B=R: (468,234,910,1142,762);
    k=-3,  B=R: (2015,325,1079,1213,1872).

Their first two coordinates independently satisfy the exact pullback
relations eta_even,0=b(a*eta0-eta1), eta_even,1=-b*eta0.
The raw integrals J_PB=integral over D_B of omega_P are consequently

    k=-15: J_PQ=1950, J_PR=377;
    k=-3:  J_PR=39.

#### Applying normalization and completing the available global entries

The normalization in 5.24 gives

    h_13(D_P,D_B)=J_PB-lambda_P,0*I_B,0-lambda_P,1*I_B,1.

It yields the three local heights stated above, with no additional
factor two or four. Add the complete away-from-13 terms of 5.21:

    k=-15:
      h(D_P,D_Q)=1950,
      h(D_P,D_R)=1248;
    k=-3:
      h(D_P,D_R)=1274.

Now apply the GLOBAL Jacobian reconstruction of 5.18. For k=-15,
R=P+Q, so the following genuine global height entries are certified:

    h_D(P,P)=923+O(13^3),
    h_D(P,Q)=1586+O(13^3).

For k=-3, R=P-Q, the certified global combination is

    h_D(P,P)-h_D(P,Q)=1417+O(13^3).

Division by four is a 13-adic unit and loses no precision. These are
partial global height matrices; the remaining mixed pairings have not
been computed. In particular h_D(Q,Q) is not supplied, and the
shared-disc k=-3 P/Q integral still requires separate treatment.

The full genus-two blocks, genus-three QC function, root isolation and
final rational-point sieve remain incomplete. None of these arithmetic
congruences constructs or excludes an original size-eleven PTE pair.


### 5.26 Completing both global genus-two height matrices

The two remaining Q/R pairings and a moved-divisor P/T pairing now
complete both global Coleman-Gross height matrices modulo 13^3=2197.
The certificates are `scripts/check_remaining_qr_heights13.py` and
`scripts/check_complete_height_blocks13.py`, using existing SymPy.
The algorithms, denominator bounds, finite terms and reconstruction
were independently reviewed and reproduced. No software was installed.

These matrices refer to the P,Q classes, monic eta basis and compatible
splittings fixed in 5.18 and 5.22. They are

    k=-15: H = [[923,1586],[1586,689]] modulo 2197;
    k=-3:  H = [[2093,676],[676,1027]] modulo 2197.

This completes these global height blocks at the stated precision.
It does not compute the genus-three local QC functions or their roots.

#### Q/R integration with split or unramified coefficient fields

For a root a of u_Q, put b=v_Q(a), u=1/(x-a), v=y*u^3/b.
Then v^2=G(u)=u^6*f(a+1/u)/b^2 is monic of degree six, and
the differential with residues (a,b)-(a,-b) is -2 eta_even,2.
Summing over both roots of u_Q gives exactly omega_Q of 5.18.
Both u_Q and u_R have unit discriminants at 13; their resultant and
ordinates are units. Thus all transformed endpoints are finite units.

For k=-3 the roots a lie in Q13 and reduce to 7 and 11. Their
residues modulo 2197 are 2035 and 102, with ordinates 969 and 336.
The five-component integrals over D_R on the two even models are

    a=2035: [520,1872,377,308,1352];
    a=102:  [663,962,1339,953,1491].

The raw contributions are 1443 and 1716, whose sum is 962 modulo 2197.
Independent endpoint calculations used both the original support algebra
and the transformed support polynomial with separately lifted Frobenius.
The two methods give the same vectors. The first two entries also agree
with the exact holomorphic pullbacks
eta_even,0=b(a eta_0-eta_1), eta_even,1=-b eta_0.

For k=-15 use the unramified quadratic field

    K=Q13[a]/(a^2+1020a+81225), sigma(a)=-1020-a.

The two roots of u_R lie in K; one is (1905+2110a) modulo 2197.
Retaining the coefficient-field Frobenius is essential: the Frobenius
numerator is G^sigma(u^13)-G(u)^13, not G(u^13)-G(u)^13.
If M has columns the reduced Frobenius images, the endpoint equation is

    sigma(I)=M^t I+rhs,
    rhs=2 sum g(endpoints)-2 sum tiny integrals.

Conjugating and substituting gives

    (Id-sigma(M)^t M^t) I = sigma(M)^t rhs+sigma(rhs).

All divisions in this linear solve are by units. Its Frobenius-square
characteristic polynomial is

    (T-169)(T^2+6T+169)^2,

so its determinant at T=1 is a 13-adic unit. An independent check
also verifies rM=13 sigma(r), with the infinity-residue row r from
5.25. The resulting integral vector, as constant/a coefficient pairs, is

    [(793,195),(1950,1196),(1508,1378),(21,87),(802,1322)].

The raw contribution from one source point is (1378,1638).
Taking its field trace gives 1716 modulo 2197.
Both holomorphic pullback identities and the original semilinear
equation hold. The coefficient arithmetic and structural Frobenius
checks were independently rerun.

The modular implementation is certified, not merely checked at two
precisions. For the retained terms n=0,...,3, vertical reduction starts
at m=6+13n<=45. Among its odd denominators 1,...,89, only 13,39,65
lose a digit. Horizontal reduction has degree at most 28 and its
denominators 2e-4 lose at most two more digits, at 26 and 52.
Thus a reduction path loses at most five digits. Starting modulo 13^20
with numerators scaled by 13^6 leaves at least nine unscaled digits;
the quadratic-field run uses 13^22. Coefficient inversions, Hensel
lifts, polynomial remainders and endpoint evaluations use only units.
The output is required only modulo 13^3. Larger guard runs agree,
but this explicit budget, not agreement alone, proves reliability.

The omitted n>=4 Frobenius terms have the six-branch-point primitive
bound proved in 5.25, since the transformed models have good reduction
and endpoints are units. Tiny terms of integrated degree at least
three have valuation at least three. These bounds also apply over
the unramified quadratic coefficient field.

Subtracting the already certified holomorphic normalizations gives

    k=-15: h13(D_Q,D_R)=1716-1301*1963-390*650=2119 mod2197;
    k=-3:  h13(D_Q,D_R)=962-1409*481-789*2158=2119 mod2197.

Adding the respective complete away terms 390 and 871 gives global
W=312 for k=-15 and W=793 for k=-3.

#### Avoiding the shared 13-adic disc by an exact divisor move

For k=-3, direct P/Q integration remains a shared-disc problem.
It is unnecessary for the matrix: use T=4P-Q instead.
Repeated exact Cantor addition of P to R=P-Q gives

    u_T=x^2-14122312959396*x/446872195225
            +119757412617129/446872195225;
    v_T=(94388884819291490706*x-1395182795689585717794)
            /298727359424984125.

Every addition is certified by interpolation V through P and the
current reduced divisor, exact division of f-V^2 by (x-33)u,
and reduction of -V modulo the new monic u. Hence this is an exact
Jacobian relation, not a congruence.

The coefficients of u_T,v_T are 13-integral; discriminant(u_T),
resultant(u_T,v_T), and u_T(33) are all 13-units. Both finite-unit
and non-Weierstrass conditions are checked. In particular, the earlier
candidate 2P-Q would NOT work with this unit-endpoint algorithm:
its ordinate resultant has valuation one. Disjoint x reductions alone
would not justify that algorithm.

Using the same rational-P even model as in 5.25 gives

    I_even(D_T)=[273,1859,143,1464,1266] modulo2197,
    raw integral = -2*143 = 1911 modulo2197.

An independent exact-rational primitive implementation gives the same
vector, with primitive traces [1898,1001,1404,1973,839] and tiny traces
[1105,676,1391,923,559]. The holomorphic integrals on the odd model are
4 I_P-I_Q=[1521,1560] modulo2197. This use of the class relation is
valid for holomorphic integrals; no additivity of second-kind coordinates
on Jacobian classes is assumed. The fixed P correction therefore gives

    h13(D_P,D_T)=1911-2021*1521-1085*1560=663 modulo2197.

#### All away contributions for the new representative

The old away terms cannot be reused for T. The complete new sum is

    away(D_P,D_T)=log13(3)-2log13(7)-2log13(61)
                    +2log13(569)+2log13(1188727)
                 =1560 modulo2197.

Here is the exhaustive regular-model calculation. The x resultant is

    u_T(33)=2*3^5*7*61*569*1188727/(5^2*133697^2).

At the good primes 7,61 the ordinates agree, while at 569,1188727
they are opposite. Each intersection has multiplicity one and nonzero
ordinate. Thus their contributions have respectively signs -2,-2,+2,+2.

The only denominator primes are 5 and 133697. At either prime one T
point goes to infinity and the other is integral. Their finite reduced
x coordinates are respectively 4 and 41333, different from P's 3 and33.
There is no intersection, including with P's conjugate.
At 11,19,53 the reductions of T and P are disjoint. The component
degree of D_P is zero by the regular models of 5.18--5.21.

At 2, in the primitive polynomial u_T(1+d), its constant and linear
coefficients have valuation one and the leading coefficient is a unit.
Both roots have v2(d)=1/2, outside the smooth P neighborhood x=1+4t
of 5.21. D_P has zero component degree, and no horizontal intersection.

At 3 use x=-3-9t and the twin chart X=(t-11)/sqrt(3), V=y/243
from 5.20. After removing the constant factor81 the T polynomial is

    G(t)=446872195225*t^2+1867060681194*t+2051187669782.

It reduces to t^2-1. P has t=-4; G(-4) has valuation one and
G'(-4) is a unit. The nearby T root lies in the same twin chart, with

    V_T=-38*(91996963761492681*t+181757575281293068)
                 /298727359424984125.

At t=2 modulo3 its ordinate is 1, opposite to V_P=44=2 modulo3.
The X difference has valuation one over Q3(sqrt(3)).
Thus the signed intersection is -2 over that quadratic extension;
norm-character descent gives +log13(3), not -log13(3).
The other T point lies in the distinct smooth cluster neighborhood.
D_P's component is invariant and its vertical correction vanishes.
This proves the full away sum, including all denominator and bad primes.

#### Global reconstruction

For k=-15, U=1950,V=1248,W=312 in the formulas of 5.18 give
the displayed first matrix. For k=-3 put Z=h(D_P,D_T)=663+1560=26.
Since D_T represents twice (4P-Q), global bilinearity gives

    h(P,P)=(Z-V)/12,
    h(P,Q)=(Z-4V)/12,
    h(Q,Q)=h(P,Q)-W/4,

with V=1274,W=793. These give the displayed second matrix.
The denominators4 and12 are 13-units. This bypasses direct shared-disc
integration without changing the global height or omitting finite terms.


### 5.27 Explicit global height quadratic forms and remaining QC gap

Let log_eta denote the two holomorphic Coleman logarithms for eta_0,eta_1.
For the genuine Jacobian generators P,Q, put
L0=(log_eta(P),log_eta(Q))/13. Since D_P,D_Q represent twice those
classes, their integral vectors must be divided by two first.
The certified matrices modulo169 are

    k=-15: L0=[[68,92],[105,89]];
    k=-3:  L0=[[126,23],[105,22]].

Both determinants are6 modulo13. With H from5.26 define
A=L0^(-t)(H/13)L0^(-1). Explicitly,

    k=-15: A=[[37,131],[131,135]] modulo169;
    k=-3:  A=[[90,114],[114,36]] modulo169.

The identities 13 L0^t A L0=H modulo2197 were independently checked.
The determinants of H/13 are respectively7 and5 modulo13; in particular
both global height blocks are nondegenerate at the computed precision.

For a rational Jacobian class D, write z=log_eta(D)/13. Then

    h(D,D)=13*z^t*A*z modulo13^3.

This identity holds on all J(Q), not just the displayed subgroup:
the rank-two upper bounds and independent P,Q from5.13 show they span
J(Q) tensor Q13. Torsion has both logarithm and height zero.
Good reduction and #J(F13)=176, prime to13, imply
log_eta(J(Q13)) is contained in13 Z13^2, so z is integral for every
rational class. Thus the asserted error is uniform on those classes.

Precision is explicit: H modulo13^3 gives H/13 modulo13^2;
the logarithms modulo13^3 give L0 modulo13^2. Its unit determinant
makes inversion lossless. Therefore A modulo13^2 determines the
displayed height to absolute precision13^3 on integral z.
Equivalently, in unscaled logarithms its coefficient matrix is A/13,
known only modulo13, not modulo13^3. No precision is gained by
rewriting the formula.


#### The explicit global term on the genus-three targets

Use the quotient notation from5.14 and5.22: f:X->E, g:X->C,
Richelot phi:J_C->J_D with polarized dual psi, and Z acting as
4 on the elliptic factor and -2 on the genus-two factor.
For a rational basepoint B0 and a rational point B on X, put

    D_B=[B-B0], e_B=f_*D_B, d_B=phi*g_*D_B.

Adjunction, the orthogonal projector decomposition
D_B=(f^*f_*D_B+g^*g_*D_B)/2, and the polarization multiplier
phi^dagger phi=[2] give the exact GLOBAL identity

    h_X(D_B,ZD_B)=2h_E(e_B,e_B)-(1/2)h_D(d_B,d_B).

It recovers the weights +8 and -8 on the pulled-back generators
of5.22; it uses genuine Jacobian classes, not doubled anti-divisors.

In affine cover coordinates (z,a,b) with y=ab, define normalized
Coleman integrals from B0 to B by

    e=(1/13) integral dz/(2b),
    w1=(1/13) integral dz/(k*a*b),
    w2=(1/13) integral z*dz/(a*b).

The elliptic map X_E=2b+2z^2+280, Y_E=2zX_E pulls
dX_E/(2Y_E) back to dz/(2b). The independently checked Richelot
transport of5.22 gives phi^*eta_0=c_0/k and phi^*eta_1=c_1,
where c_0=dz/y,c_1=z dz/y. Thus these are precisely the logarithms
of e_B and d_B, divided by13; no extra factor two applies.

The elliptic coefficient from5.17 satisfies
13h_E(P,P)/log_E(P)^2=6699432669=72 modulo169.
Consequently the explicit global quadratic term is

    h_X(D_B,ZD_B)
      =13*(144e^2-(1/2)*(w1,w2)*A_k*(w1,w2)^t) modulo2197.

The orders of the good reduction groups are12 for E and176 for J_D,
both prime to13, so e,w1,w2 are integral. All divisions by2 are units.
This proves the stated precision for rational points. The expression
extends analytically to local points, but this extension is not the
local-height side of the quadratic-Chabauty equation.


These are the genus-two GLOBAL quadratic height expressions needed
for the genus-three product construction. The compatible genus-two
weight is -8 as proved in5.22. They are not by themselves the functions
on residue discs that quadratic Chabauty requires. Still missing are
the basepoint-normalized genus-three local height functions and
representative corrections, the complete allowed non-13 local-value
sets for arbitrary curve points, certified power-series tails and root
isolation on all twelve discs, and the final rational-point sieve.
No finite QC root set or original eleven-distinct-integer construction
or general impossibility theorem is claimed.



### 5.28 A correctly centered QC divisor and explicit local differential

This section constructs the actual curve-level height expression.
The global quadratic term of5.27 alone is not a quadratic-Chabauty
equation: a basepoint-dependent linear correction has to be included.
The correction is explicit here, using only the already certified
height matrices. No numerical local13 series or roots are claimed.

Write X_delta as W^2=F_delta(t) using5.14, let iota(t,W)=(t,-W),
and retain the elliptic involution sigma of5.15. Let H=infinity+
+infinity-, R_iota=Fix(iota), R_sigma=Fix(sigma). The exact identities are

    div(dt/W)=2H, div(W)=R_iota-4H,
    div(M(t))=R_sigma-2H.

In particular K_X=2H with this differential. The two points in H are
rational; H/2 below denotes a rational-coefficient divisor, not a
claimed rational point.

Use actual rational basepoints

    B_1=(15/2,9025/16), B_5=(5/6,9025/144).

Both are outside the two ramification loci, and sigma(B_delta)=infinity+.
Their orbit under iota,sigma consists of B,iota(B),infinity+,infinity-.

#### The actual correspondence divisor

Use the proper graph correspondence

    Gamma=Gamma_iota-3 Gamma_sigma.

It induces -Z, not Z, on the Jacobian. Its fiber degrees are -2 and
its diagonal intersection has degree8-3*4=-4; thus the required
trace-zero intersection is -4-(-2)-(-2)=0.
It meets the diagonal properly, avoiding a self-intersection convention.

For a varying point P define the degree-zero divisor

    D_Gamma(B,P)=R_iota-3R_sigma-iota(P)+3sigma(P)
                             -iota(B)+3sigma(B).

This is the divisor obtained from the diagonal restriction minus the
two fiber restrictions, not simply Z[P-B]. Its class equals

    Z[P-B]+C_B,  C_B=2Z[B-H/2].

Indeed, R_iota is linearly equivalent to4H, R_sigma to2H,
iota(B) toH-B, and sigma(H) toH. These identities give
C_B=2B+6sigma(B)-4H in the Jacobian.

Put j_H(P)=[P-H/2] in J tensor Q and q(A)=h_X(A,ZA).
Symmetry and self-adjointness of Z now give the exact global formula

    h_X([P-B],[D_Gamma(B,P)])=q(j_H(P))-q(j_H(B)).

The right side contains the previously missing linear correction when
expressed in logarithms from B to P. No new global height coefficients
are needed.

The correspondence/divisor construction and finite-local-value theorem
used here are Balakrishnan--Dogra, *Quadratic Chabauty and Rational
Points I*, Definition6.2 and Theorems6.3 and1.2:
[primary text](https://arxiv.org/html/1601.00388).
Their theorem applies to this corrected divisor, with the established
rank and logarithm hypotheses. It does not justify using the uncorrected
self-pairing instead.

#### Centered logarithms and the global term

Since 2j_H(P)=[P-iota(P)], its logarithm is half the integral from
iota(P) to P. In the affine quotient coordinates (z,a,b) of5.14 set

    e_H(P)=(1/26) integral_iota(P)^P dz/(2b),
    w_H(P)=(1/26) (integral_iota(P)^P dz/(k*a*b),
                         integral_iota(P)^P z*dz/(a*b)).

Then, with A_k from5.27,

    q(j_H(P))=13*(144 e_H(P)^2-(1/2)w_H(P)^t A_k w_H(P))
                  modulo13^3.

The corrected global term is this expression minus its value at B.
Alternatively, translate the B-to-P logarithms of5.27 by the fixed
centered logarithms of B. It is incorrect to omit that translation.
The factors1/2 and the precision follow from5.27; normalized logarithms
are integral, and A_k is known modulo13^2.

#### An explicit residue-correct local kernel, including infinity

For a finite point Q=(t_Q,W_Q), define

    nu_Q=(1/2)*(1+W_Q/W)*dt/(t-t_Q).

Its residue divisor is Q-H/2. The apparent pole at iota(Q) cancels,
and each infinity has residue-1/2.
For Q=infinity^epsilon, epsilon=+1 or-1, put

    nu_Q=-epsilon*c*t^3*dt/(2W),
    c=1 for delta1, c=5625 for delta5.

The residues are +1/2 on that infinity and -1/2 on the other.
This infinity formula is necessary: sigma(B)=infinity+ for our
actual basepoints, so a finite-point formula alone would not suffice.

Put g=W/M^3. A differential with exactly residue divisor D_Gamma is

    omega_(B,P)=dlog(g)+3nu_sigma(B)-nu_iota(B)
                         +3nu_sigma(P)-nu_iota(P).

The principal part contributes R_iota-3R_sigma+2H; the four nu terms
contribute the rest minus2H, giving exactly the stated residue divisor.
For our actual basepoint, sigma(B)=infinity+, so the required residue
there is3; it must not be incorrectly canceled. There are no hidden
poles at B or P on the disjoint-support open.
This gives a fully explicit rational differential before normalization.

Let Psi be the Coleman--Gross de Rham projection and pr_hol projection
onto holomorphic forms along the compatible complement fixed in5.22.
The normalized differential is

    eta_(B,P)=omega_(B,P)-pr_hol(Psi(omega_(B,P))),

and the actual local height is integral_B^P eta_(B,P).
This specifies the normalization; the moving de Rham projection and
the resulting certified power series have NOT yet been computed.
The logarithmic-exact term can also be integrated as
log13(g(P)/g(B)), but that term must not be silently omitted.

Define the basepoint-normalized function

    rho_B(P)=h13(P-B,D_Gamma(B,P))
                       -q(j_H(P))+q(j_H(B)).

For rational P it belongs to the negative sum of the away local values.
The theorem makes that target set finite. Formulae alone do not enumerate
the target set or isolate roots of rho_B.

#### Rational exceptional points and a normalization stress test

For direct disjoint-support evaluation exclude the ramification loci,
H, B and iota(B); the limiting value at B is zero.
The two ramification loci have no rational points: F_delta(t)>0 for
real t, since N,D have no common root, and the discriminants of M are
-380 and-1520. Thus the omitted rational points are precisely the known
four-point orbit, not an unexamined set of possible rational exceptions.

For delta1 this orbit lies over z=infinity and maps to O or(0,0) on E.
For delta5 it lies over z=5/2, with elliptic images(180,900)=P_E or
(405,2025)=-P_E+(0,0). These have integer indices0 or plus/minus1,
already excluded as valid center lifts by the exact collision test5.9.
This does not classify any other rational points on X.

One can check the diagonal trivialization separately. The fiber-degree-zero
cycle Delta+3Gamma_sigma-2(H cross X+X cross H) restricts to div(M^3).
Proper rigidification at H/2 instead has diagonal trivializer
(sigma(t)-t)^3, since

    div(sigma(t)-t)=R_sigma-H-sigma(H).

The two choices differ by explicit vertical principal divisors and must
not be mixed. The graph kernel above fixes the actual local-height sign
through D_Gamma, without guessing it from a diagonal line bundle.

There is a concrete reason normalization matters. At the good prime7,
Hensel lifting at infinity supplies Q7 points with t=7^(-n) for every
n>=1. The leading coefficients of W and M are units. Thus
v7(W/M^3)=2n, and omitting its character term would leave the unbounded
sequence -2n log13(7), up to a fixed constant. This is a statement about
local points, not infinitely many rational points over Q.
The corrected D_Gamma family has the finite-value property instead.


### 5.29 Away-height targets: only primes 2 and3 remain

For both curves, the corrected curve-level local height of5.28 is
identically zero at every finite prime other than2,3,13.
In particular its local target is {0} at5,11,19 and53, even though
those primes occur in the displayed model discriminants.
The targets at2 and3 are not yet computed.

The exact certificate is
`scripts/check_qc_normalization_and_local_targets.py`.
Both ramified and unramified model calculations, accessible components,
and the graph-symmetry argument received independent review.
This result concerns the corrected CURVE-level height. It does not
erase the nonzero genus-two generator-pair contributions in5.18--5.26.

#### Good primes and the graph-height input

The discriminants of the degree-eight polynomials are

    delta1: 2^52*3^33*5^30*11^2*19^12*53;
    delta5: 2^108*3^33*5^58*11^2*19^12*53.

Their leading coefficients add no exceptional primes. Therefore the
curves have good reduction outside S={2,3,5,11,19,53}; the correctly
normalized local height vanishes there by the theorem cited in5.28.

For bad-prime computations we use Theorem3.2 of Balakrishnan, Dogra,
Muller, Tuitman and Vonk, *Quadratic Chabauty for modular curves:
algorithms and examples*:
[primary article](https://www.cambridge.org/core/journals/compositio-mathematica/article/quadratic-chabauty-for-modular-curves-algorithms-and-examples/C317E63BB0128EB3F2163621CA63E66E).
It expresses the normalized local height through a graph potential,
specified by its Laplacian and value zero at the basepoint vertex.

The needed symmetry consequence is elementary. If an involution
preserves the metric graph and commutes with the correspondence, it
preserves the graph measure: signed edge permutations preserve the
projection onto graph homology, and component traces are invariant
under conjugation. Hence j composed with sigma minus j is harmonic,
thus constant. Applying sigma twice forces that constant to be zero.
Therefore j is constant on each sigma-orbit of vertices, even when
sigma does not fix the basepoint. If all accessible vertices are in
one such orbit, their normalized height is zero.

#### Primes11 and53: a single component

The following are all double-root nodes, listed as
(t, F(t)/ell, F''(t)/2) moduloell:

    delta1, ell11: (5,4,5), (6,3,4);
    delta5, ell11: (3,3,3), (4,9,4);
    delta1, ell53: (38,12,12);
    delta5, ell53: (20,50,49).

Every last two entries is a unit. The residual polynomial after removing
the double factors is squarefree and nonconstant, hence not a square
over the algebraic closure. Infinity is smooth. Thus the displayed proper
models are regular semistable with one geometric component and ordinary
nodes of thickness one. A section cannot specialize to such a node:
after an unramified splitting its local equation is xy=ell, whereas
both coordinates would have positive integral valuation. All rational
local points therefore hit the single vertex, so the target is {0}.

#### Prime19: two exchanged elliptic components

Let a=-2 for delta1 and a=4 for delta5, and write
F(a+x)=sum(c_j x^j). In both cases the coefficient valuations are

    [4,3,2,1,0,0,0,0,0].

The inner chart t=a+19u, V=W/19^2 has smooth quartic reduction

    delta1: 8u^4-2u^3-3u^2+8u-9, discriminant11 mod19;
    delta5: 8u^4+u^3+2u^2+2u+1, discriminant7 mod19.

The outer normalizations are respectively

    x^4+3x^3-3x^2+4x+8, discriminant7;
    x^4-8x^3-8x^2+x+8, discriminant11.

An annular chart is x*r=19, w=W/x^2, with

    w^2=sum_(j<4)(c_j/19^(4-j))*r^(4-j)
                      +c_4+sum_(j>4)c_j*x^(j-4).

It is integral by the valuation list and agrees with the original
equation on the generic fiber. At x=r=0 it has w^2=8 modulo19.
After an unramified quadratic extension, 2w is a unit and the two
nodes have exactly the local equations x*r=19. These charts form
the normalization of the corresponding blowup of the t-line.
They give a proper regular semistable model with precisely two
smooth genus-one components joined by two thickness-one nodes.

The involutions on the inner charts are

    sigma1(-2+19u)=(15u-10)/(2u-1),
    sigma5(4+19u)=(25u+4)/(30u+5).

Their nonconstant reductions show that sigma exchanges the inner
and outer components. It commutes with Z, so the graph argument
makes both vertex values equal; normalization makes both zero.
All Q19 sections hit vertices, not thickness-one nodes. Thus the
target is {0}, without asserting potential good reduction at19.

#### Prime5: two different accessible orbits, one for each twist

Use K=Q5(pi), pi^2=5. On X1 set t=pi*u, W=25V.
The resulting monic polynomial G=F1(pi*u)/pi^8 is integral and reduces to

    (u-1)^2(u+1)^2(u^2+2u-1)(u^2-2u-1).

There are exactly two ordinary nodes, at u=plus/minus1; all other
branch roots are simple and infinity is smooth. Its discriminant has
K-valuation 2*30-8*7=4. The involution acts as

    u -> (3pi*u-32)/(2u-3pi),

and exchanges the nodes. Their thicknesses are therefore equal;
the total discriminant valuation shows that each has thickness two.
An independent direct check gives v_K(G),v_K(G'),v_K(G'')=(2,2,0)
at each node. The critical point moves only by O(pi^2), so its
smoothing value retains valuation two, confirming the thickness.
Resolving each node introduces a single exceptional component.
The main component is preserved and the exceptional pair is exchanged.

For X1(Q5), finite t has integer 5-adic valuation, so
v_K(u)=2v5(t)-1 cannot be zero. Such points cannot reduce to u=plus/minus1;
the points at infinity are smooth as well. Thus every X1(Q5) point
meets the main component. The target is zero by normalization.

For X5(Q5), use the K-isomorphism of covers a1=pi*a5 with the same z,b.
Its X1 parameter satisfies u=a5-2z/pi. The equation
5a5^2=4z^2+30z+80 forces v5(z)>=1: for v5(z)<=0 the right side has
even valuation, a contradiction. Then a5 is a unit and a5^2=1 mod5,
so u reduces to plus/minus1. There are no Q5 points above z=infinity,
since their slope would require sqrt5. This also covers the rational
hyperelliptic infinities of X5, whose z-coordinate is finite.

Over K a section through a thickness-two node, locally xy=pi^2,
has both valuations one. After resolution it meets the unique
exceptional component, not a component intersection. All X5(Q5)
points therefore meet the two exchanged exceptional vertices.
They form a single sigma-orbit containing the basepoint, so their
normalized heights are also zero. Descent by the local extension
degree preserves this vanishing. The argument does NOT confuse
the main-component orbit of X1 with the exceptional orbit of X5.

#### What remains before root isolation

The complete negative away sum for the corrected family is now drawn
from the two uncomputed finite value sets at2 and3. At3 there are
nested branch clusters and multiple potentially accessible components;
the involution preserves rather than exchanges the two initial clusters.
The prime2 also requires a genuinely wild regular-model analysis.
Neither target set is inferred from generator heights or bounded search.

At13,5.28 gives the exact rational kernel and normalization prescription,
but the moving de Rham projections, basepoint constants, certified
series and tail bounds are still missing. Therefore no numerical QC
equation, complete isolated-root list, or final Mordell--Weil sieve
has been produced. The original distinct-entry size-eleven task
remains unresolved.


### 5.30 Complete away targets and infinity projection

[COMPLETE_AWAY_TARGETS.md](COMPLETE_AWAY_TARGETS.md) completes the local sets
left open in5.29. The exact QC targets are
Omega1={0,-1/6,-2/3}log13(3) and
Omega5={0,1/2,2/3,-4/3}log13(3).
[INFINITY_PROJECTION.md](INFINITY_PROJECTION.md) computes the de Rham
projection at infinity and its compatible holomorphic normalization modulo13^3.
Both have independently reviewed proofs and executable certificates.

### 5.31 Exclusion of the sixth shift-chain mode

**Theorem.** Under the established reductions and arithmetic inputs of
5.7--5.30, the sixth shift-chain type of5.5 has no valid distinct-entry
size-eleven PTE lift.

**Proof.** The complete proof, analytic estimates, dependency ledger and
independently reviewed exact calculations are in
[SIXTH_MODE_QC_EXCLUSION.md](SIXTH_MODE_QC_EXCLUSION.md).
Here is the finite conclusion and its scope.

The mod240 elliptic restriction and squareclass parity retain four residue
discs on each X_delta. Moving-source differentiation computes their
third-kind normalization. On the B_delta disc write t=tB+13u. The corrected,
centered QC function satisfies
rho1/13=6u+6u^2 and rho5/13=12u+6u^2 modulo13.
All roots of the permitted target equations are simple: respectively four
and eight on these two discs. Integral analytic tails prove completeness;
this is not a bounded point search. The explicit graph-divisor cocycle and
invariant isotropic splitting transport these roots to the other six discs,
giving48 local QC roots.

The exact elliptic group and necessary center-square tests at61,101,107
restrict n modulo3120 to {0,1,1560} for epsilon0 and {0,1560,3119} for
epsilon1. On B1 parity forces n=0 modulo13; on B5 it forces n=12.
The elliptic logarithm gives n=11u and n=12+11u respectively modulo13.
Only u=0 survives on each disc. Simple-root uniqueness identifies it with
B_delta. Symmetry leaves only the eight known orbit points, all excluded by
the exact collision tests in5.9. This proves the theorem.

Certificates are scripts/check_sixth_mode_qc_roots13.py (existing SymPy) and
scripts/check_sixth_mode_candidate_sieve.py (standard-library Python).
The latter independently verifies the elliptic formal logarithm and uses
the full certified elliptic group; no unstated genus-three saturation is
needed.

This supersedes the unfinished sixth-mode QC/root/sieve status in earlier
sections. It does not classify all rational points on X1 or X5, since the
remaining center-square conditions are used in the final sieve. It does
not exclude the other17 shift-chain types or solve the unrestricted
size-eleven PTE existence problem.

### 5.32 Finiteness for all eighteen symmetric shift-chain types

**Theorem 14 (complete normalized shift-lift finiteness).** There are only
finitely many pairs of disjoint rational sets `A=-A`, `B=-B`, each of size
ten, with equal moments through degree nine, for which the shift `t=1`
and cancellation construction of Section 5.1 produces two sets of eleven
distinct entries. Consequently, allowing an arbitrary common center and
nonzero rational shift yields only finitely many affine-equivalence
classes of successful inputs and outputs for this construction.

This theorem covers all eighteen types, including the twelve types whose
finiteness was previously unproved. It does not enumerate the finite set,
provide an effective height bound, or decide whether any valid size-eleven
output exists. In particular, it excludes an infinite family of
inequivalent solutions obtained by this construction; the unrestricted
size-eleven problem remains open.

**Proof.** Section 5.1 gives an exhaustive finite list of maximal chain
types. We number them in the order produced by its partition enumeration
(increasing `m`, number of A chain pairs, then the two length lists).
The older “fifth mode” is type 2 below; the older “sixth mode” excluded
in Section 5.31 is type 10, not type 6 of this table.

| Type | m | A pair lengths | B pair lengths | Necessary curve genus | Argument |
|---|---:|---|---|---:|---|
| 1 | 1 | (4) | (1,1,1,2) | 2 | Sections 5.2–5.3 |
| 2 | 1 | (1,3) | (1,1,3) | 3 | Section 5.4 |
| 3 | 1 | (1,3) | (1,2,2) | 3 | Below |
| 4 | 1 | (2,2) | (1,1,3) | 2 | Below |
| 5 | 1 | (2,2) | (1,2,2) | 2 | Below |
| 6 | 1 | (1,1,2) | (1,4) | 4 | Below |
| 7 | 1 | (1,1,2) | (2,3) | 4 | Below |
| 8 | 1 | (1,1,1,1) | (5) | 4 | Below |
| 9 | 2 | (3) | (1,1,1,2) | 2 | Sections 5.2–5.3 |
| 10 | 2 | (1,2) | (1,1,3) | 3 | Section 5.5; exclusion in 5.31 |
| 11 | 2 | (1,2) | (1,2,2) | 3 | Below |
| 12 | 2 | (1,1,1) | (1,4) | 2 | Below |
| 13 | 2 | (1,1,1) | (2,3) | 2 | Below |
| 14 | 3 | (2) | (1,1,1,2) | 2 | Sections 5.2–5.3 |
| 15 | 3 | (1,1) | (1,1,3) | 2 | Below |
| 16 | 3 | (1,1) | (1,2,2) | 2 | Below |
| 17 | 4 | (1) | (1,1,1,2) | 2 | Sections 5.2–5.3 |
| 18 | 5 | () | (1,1,1,1,1) | 3 | Below |

For each type every valid configuration maps, with finite fibers, to a
fixed smooth projective curve over Q of genus at least two. Faltings's
theorem [6, Theorem 7] makes its rational-point set finite. The finite
union over the eighteen types is finite. Normalization by the common
symmetry center and the nonzero shift then proves the affine statement.

Here are the remaining twelve reductions, including all exceptional
branches. Every lower-case center variable denotes a square of a rational
center; `S,T` and `U,V` denote the sum and product of the specified pair
of squared centers. Write

\[
 M_L(k,z)=\sum_{j=0}^{k}\binom{2k}{2j}z^{k-j}
                    \sum_{d\in O_L}d^{2j},\qquad
 C_m(k)=\tfrac12\sum_{d\in O_{2m}}d^{2k}.
\]

For two squared centers use `p_0=2`, `p_1=S`,
`p_k=S p_(k-1)-T p_(k-2)`; substitute these into `M_L` when taking a
pair of centers. For three or four singleton squares, Newton's identities
determine the unique monic cubic or quartic from its initial power sums.
All displayed eliminations start with the exact four even moments of
orders 2,4,6,8. Resultants stated up to proportionality differ only by a
nonzero rational constant.

#### Type 3: (1,(1,3),(1,2,2))

Let a,b be the A length-one and length-three squared centers; let c be
the B singleton square and S,T the sum/product of the two B length-two
squares. The first two equations determine

\[
 a=2S-3b+c-5/4,
\quad T=-\frac{16S^2-96Sb+32Sc-64S+96b^2-48bc+156b-20c+27}{32}.
\]

The resultant in S of the two remaining equations factors as

\[
 (4b-25)(4b-9)(4b-1)(4c-1)^3
 ((b-c-1)^2-4c)^2 H_3(b,c),
\]

where

\[
 H_3=1024b^2-1520bc+636b+640c^2-1860c+2025.
\]

The c=1/4 factor puts a B vertex in the centered A chain. The b=1/4
and b=9/4 factors repeat A vertices. At b=25/4 the prescribed A chains
have additional unit edges, so they cannot be the specified maximal
chains; equivalently they contradict the degree-nine edge bound for a
nontrivial disjoint seed. Finally the squared factor says that the B
singleton is one of the vertices of an A length-three chain. Thus a
valid input lies on H3(b,c)=0.

Writing b=u²,c=v² and homogenizing gives a smooth plane quartic. A useful
smoothness certificate is that every nonempty principal minor of

\[
\begin{pmatrix}2048&-1520&636\\-1520&1280&-1860\\636&-1860&4050\end{pmatrix}
\]

is nonzero: they are 2048,1280,4050,311040,7889904,1724400,−2747105280.
Indeed, at a singular point of a quartic quadratic in U²,V²,W², the
squares of its nonzero coordinates are killed by the corresponding
principal submatrix. This is impossible here. The quartic has genus
three. For each b,c the third moment equation is a cubic in S with
nonzero constant leading coefficient. It gives at most three S, then
fixes a,T and all remaining squared centers. Thus the fibers are finite.

#### Type 4: (1,(2,2),(1,1,3))

Let S,T be the sum/product of the A length-two squares; let c be the B
length-three square, and U,V the sum/product of its two singleton squares.
The first two equations give

\[
 U=2S-3c-3/4,\qquad
 V=(8S^2-48Sc-24S+16T+48c^2+66c+9)/8.
\]

The resultant in T of the remaining equations is proportional to

\[
 (4c-9)(4c-1)(-6S+20c+25)(-2S+4c+1)
 [4S^2-16Sc-20S+16c^2+24c+25].
\]

Values c=1/4,9/4 give collisions with the central chain. For the linear
branch S=2c+1/2 the third equation forces T=(4c−1)²/16; the A centers
are ±(v−1/2),±(v+1/2), where c=v², and their length-two chains overlap
the B length-three chain. On the quadratic branch,
S=2v²+5/2±2v and T=(2v∓1)²(2v±3)²/16. These centers also overlap the
B chain. The coefficients of T in these equations vanish only at c=−11/4
or c=1/4,9/4. The former is not a rational square; the latter are already
excluded. Hence no positive-dimensional exceptional branch was discarded.

Every valid input must therefore satisfy

\[
 S=(20c+25)/6,\quad T=(240c^2+40c+167)/48,\quad
 U=(44c+91)/12,\quad V=(256c^2+97c+466)/36.
\]

The coefficient determining T is a nonzero multiple of 4c+11. Since
c,T,V are all rational squares, so is their product. This gives a
rational point on

\[
 Y^2=3c(240c^2+40c+167)(256c^2+97c+466).
\]

The right side has degree five and is squarefree, so the curve has genus
two. All center squares are determined by c and two quadratic equations;
the rational-point fibers are finite.

#### Type 5: (1,(2,2),(1,2,2))

Let `c=w^2` be the B singleton-center
square. Let `S,T` be the sum and product of the two A length-two
center-squares, and `U,V` those of the two B length-two center-squares.
The first two moments give

```
U = (8S-4c+1)/8,
V = (-64Sc+16S+128T+48c^2-56c+11)/128.
```

The third residual is
`-3(4c-1)(-96Sc+88S+128T+80c^2-200c+173)/512`.
The factor `c=1/4` puts the B singleton in the centered A chain. Otherwise
solve for `T`; the fourth residual becomes
`(4c-1)(4c-9)(4c-25)(16S-20c+85)/4096`.
At `c=9/4` the A squared-center polynomial has root 1, so that its length-two
chain overlaps the centered A chain. At `c=25/4` it has root 4, so that
an A vertex coincides with the B singleton. The remaining branch is

```
S=(20c-85)/16,  T=(80c^2-840c+589)/256,
U=(12c-83)/16,  V=(4c-49)(4c-9)/256.
```

Because `c,T,V` are rational squares, this gives a rational point on

```
Y^2=c(80c^2-840c+589)(4c-49)(4c-9).                 (R5)
```

The right side is squarefree of degree five, so its smooth projective
model has genus two. Each `c` fixes both center quadratics and hence
only finitely many configurations.

#### Type 6: (1,(1,1,2),(1,4))

Let a be the A length-two square, S,T its singleton sum/product, b the B
singleton square and c the B length-four square. The first two equations
give

\[
 S=-2a+b+4c+17/4,
\quad T=(12a^2-8ab-32ac-28a+16bc+17b+24c^2+8c+16)/4.
\]

The resultant in b is proportional to

\[
 (a-c)((a-c-1)^2-4c) E_6(a,c)H_6(a,c),
\]

where

\[
\begin{split}
 E_6&=3a^2-8ac-7a+6c^2+2c+4,\\
 H_6&=a^3-5a^2c-9a^2+10ac^2+15ac+24a
                         -10c^3+10c^2-20c-16.
\end{split}
\]

The first two factors make the A length-two chain overlap a B length-four
chain. The third equation is linear in b with coefficient proportional
to E6. If E6=0, its constant coefficient also has to vanish. Their
resultant in a is a nonzero multiple of c(c−4)²(c−1)³. For c=0 the
reflected B chains coincide; for c=1 or4 they meet the centered A chain.
Thus E6 is nonzero on every valid input, b is uniquely determined, and
H6(a,c)=0.

The projective cubic H6=0 is smooth: its affine singularity ideal
(H6,∂aH6,∂cH6) is the unit ideal over Q, and its intersection with the
line at infinity is given by
A³−5A²C+10AC²−10C³, with discriminant −200 at C=1.
Its three affine intersections with a=0 are simple, since
H6(0,c)=−2(5c³−5c²+10c+8) has discriminant −1483200.
At each of its three distinct points at infinity, the function a=A/Z
has a simple pole; A is nonzero there. The double cover a=u² therefore
has six simple branch places. It is geometrically connected (a has odd
valuation at these places), and its smooth completion has genus four
by Riemann–Hurwitz. Every valid rational configuration supplies a
rational point on this cover. With a,c fixed, the formula for b and then
S,T give finite fibers.

#### Type 7: (1,(1,1,2),(2,3))

Let a and S,T describe the A length-two pair and singleton squares;
let b,c be the B length-two and length-three squares. The first two
moments fix

\[
 S=-2a+2b+3c+7/4,
\quad T=\frac{48a^2-64ab-96ac-32a+16b^2+96bc+32b+48c^2-12c+9}{16}.
\]

The resultant in b of the remaining equations is proportional to

\[
 (4c-9)(4c-1)^2(8a-12c-5)
 [16a^2-32ac-8a+16c^2-8c+1]^2 H_7(a,c),
\]

where

\[
 H_7=256a^3-1280a^2c-64a^2+2560ac^2-2560ac-2976a
                         -960c^3+3920c^2+780c+3159.
\]

The c factors put B vertices in the central A chain. The squared
factor is 16[a−(v+1/2)²][a−(v−1/2)²], c=v², which forces a cross-color
overlap. On the remaining linear branch a=(12c+5)/8, the third moment,
away from c=1/4,9/4, forces b=(20c+19)/24. The fourth moment then is a
nonzero multiple of (4c−9)(4c−1)²(4c+11). Its only remaining possibility
c=−11/4 is not a rational square. Hence every valid input satisfies
H7(a,c)=0.

This cubic is smooth over the algebraic closure: the affine singularity
ideal is the unit ideal; the cubic at infinity has discriminant
−1630745395200, and H7(0,c) has discriminant −1165197783859200.
Neither affine intersection set has repeated roots, and A≠0 at infinity.
As for type6, the cover a=u² has six simple branch places and genus four.
For fixed a,c, the third moment is a polynomial of degree at most two
in b. All its coefficients can vanish only when
8a−12c−5=0 and (4c−1)(4c−9)=0, which are already excluded. Thus the
map to this genus-four curve has finite fibers even at a vanishing
quadratic leading coefficient.

#### Type 8: (1,(1,1,1,1),(5))

Put `b=v^2` for the B length-five
center. The four A singleton-center squares have prescribed power sums
`p_k=M_5(k,b)-C_1(k)` for `1<=k<=4`, so their monic quartic `K_b` is
uniquely determined. Its discriminant is

```
disc(K_b) = 5(4b-1)^2 h8(b)/2^24,
h8(b) = 26214400b^10 - 904396800b^9 + 13285785600b^8
 -130282291200b^7 + 1026147491840b^6 -6126693150720b^5
 +21780357808640b^4 -27633876692480b^3 +12983200818768b^2
 +1848981232152b +6457704345.
```

The case `b=1/4` gives a B vertex in the centered A chain and is invalid.
Complete rational splitting makes the quartic discriminant a rational
square, giving `Y^2=5h8(b)`. This squarefree degree-ten model has genus
four. The uniquely determined quartic and finitely many square-root
choices give finite fibers.

#### Type 11: (2,(1,2),(1,2,2))

Let a,b be the A singleton and length-two squares, c the B singleton
square, and S,T the sum/product of the B length-two squares. Then

\[
 a=2S-2b+c-2,
\quad T=-(2S^2-8Sb+4Sc-11S+6b^2-4bc+11b-4c+9)/4.
\]

Both remaining moments contain a factor −S+b+1, which makes a=c and
therefore repeats a vertex between colors. Remove that factor. Their
resultant in S is proportional to

\[
 (4c-9)(4c-1)
 [16b^2-32bc-8b+16c^2-8c+1]H_{11}(b,c),
\quad H_{11}=16b^2+80bc-140b-80c^2+280c-281.
\]

The first two factors put a B singleton in the centered length-four A
chain, and the quadratic factor puts it in the A length-two chain.
Thus H11(b,c)=0. Its even plane quartic, with b=u²,c=v², has derivative
matrix

\[
\begin{pmatrix}32&80&-140\\80&-160&280\\-140&280&-562\end{pmatrix},
\]

whose nonempty principal minors are 32,−160,−562,−11520,−37584,11520,829440.
They are all nonzero, so the quartic is smooth of genus three. The third
moment after cancellation is a quadratic in S with nonzero constant
leading coefficient. Hence all fibers are finite.

#### Type 12: (2,(1,1,1),(1,4))

Put `b=v^2` for the B length-four
center and `c=w^2` for its singleton center. The three A singleton squares
have `p_k=M_4(k,b)+c^k-C_2(k)`. Their fourth Newton residual is
`-4b(b-4)(b-1)(b+4c-10)`. The possibilities `b=0,1,4` force a B vertex
into the centered A chain. Consequently `c=(10-b)/4`. The resulting
cubic has discriminant `-h12(b)/64`, where

```
h12(b)=50b^6-1650b^5+26675b^4-187600b^3+49540b^2-34390b-81.
```

Thus every valid configuration maps with finite fibers to the smooth
genus-two curve `Y^2=-h12(b)`.

#### Type 13: (2,(1,1,1),(2,3))

Put `b=v^2` and `c=w^2` for the B
length-two and length-three centers. Now
`p_k=M_2(k,b)+M_3(k,c)-C_2(k)`. The fourth Newton residual is
`-(b-1)(4c-9)(4c-1)(6b+4c-25)/8`. The first three factors force B/A
collisions. Hence `b=(25-4c)/6`. The cubic discriminant is `-h13(c)/62208`,
where

```
h13(c)=921600c^6-24832000c^5+192608000c^4-644627200c^3
       +936509680c^2-438571880c-82872167.
```

Since `62208=3*144^2`, a square discriminant yields the squarefree
sextic curve `Y^2=-3h13(c)`, again of genus two with finite fibers.

#### Type 15: (3,(1,1),(1,1,3))

Let `c` be the B length-three center
square, `S,T` the symmetric functions of the A singleton squares, and
`U,V` those of the B singleton squares. Eliminating the first two
moments gives

```
U=S-3c+27/4,
V=(-48Sc+108S+16T+96c^2-228c+27)/16.
```

The third residual is
`-3(4c-9)(-96Sc+12S+48T+160c^2+160c-47)/64`.
Here `c=9/4` makes the B chain overlap the centered A chain. Eliminating
`T`, the fourth residual becomes
`(4c-25)(4c-9)(4c-1)(-8S+20c+15)/64`.
The three constant choices of `c` all cause B/A collisions. The remaining
branch is

```
S=(20c+15)/8,   T=(160c^2-20c+49)/96,
U=(69-4c)/8,    V=(2c-23)(4c-31)/48.
```

The square condition on `cTV` gives the squarefree genus-two curve

```
Y^2=2c(2c-23)(4c-31)(160c^2-20c+49).              (R15)
```

The same two quadratic polynomials prove finite fibers.

#### Type 16: (3,(1,1),(1,2,2))

Now `c` is the B singleton-center
square, `S,T` refer to the A singleton squares, and `U,V` to the two B
length-two center-squares. The first two moments give

```
U=(4S-4c+31)/8,
V=(-16S^2-32Sc+296S+64T+48c^2-296c-73)/128.
```

Write `D=4S+4c-51`. For `D!=0`, the third moment gives

```
T=(64S^3+64S^2c-944S^2+192Sc^2-1376Sc+1996S
   -320c^3+2320c^2-1996c-441)/(64D).
```

The fourth residual then factors as

```
(4c-25)(4c-1)(-4S+4c+1)(-4S+4c+9)(-4S+4c+49)
 (12S+20c-195)/(1024D^2).
```

At `D=0`, the third residual is
`-3(4c-25)(4c-21)(4c-1)/64`. The values `c=1/4,25/4` are direct
B-singleton/A-central collisions; `c=21/4` is not a rational square.
Thus a valid configuration has `D!=0`. The same collisions exclude the
first two factors of the fourth residual. For `S=c+j/4`, with
`j=1,9,49`, substitution gives the A squared-center polynomial
`(X-c)(X-j/4)`. Hence an A singleton coincides with the B singleton;
these branches are invalid. What remains is

```
S=(195-20c)/12,  T=(160c^2-1335c+2304)/36,
U=4(9-c)/3,      V=2(c-6)(4c-9)/3.
```

Again `cTV` is square, giving

```
Y^2=6c(c-6)(4c-9)(160c^2-1335c+2304).             (R16)
```

Its right side is squarefree of degree five, so this is genus two.
The fixed center quadratics give finite fibers, including possible
multiple-root parameter values (which cannot yield a valid seed).

#### Type 18: (5,(),(1,1,1,1,1))

Double all centers. The centered A chain then has positive squared
vertices 1,9,25,49,81. Equality of the first four even moments fixes
the first four elementary symmetric coefficients of the five B center
squares. Their product is z² for a rational z. Consequently they are
the roots of

\[
 Q_z(X)=X^5-165X^4+8778X^3-172810X^2+1057221X-z^2.
\]

Complete rational splitting forces the discriminant of this polynomial
to be a rational square. That discriminant is

\[
\begin{split}
 d(z)={}&3125z^8+108276547500z^6-1117199475459133506z^4\\
       &+533350028141359954915500z^2
         +2470517618419162907684321203125.
\end{split}
\]

It is squarefree of degree eight. Thus Y²=d(z) has genus three.
For each z the polynomial fixes the five center squares, giving finite
fibers. The auxiliary square z is essential: discriminant alone as a
function of an unrestricted constant coefficient would only give a
quartic hyperelliptic model of genus one.

#### Verification and scope

The exact certificate `scripts/check_all_shift_chain_types.py` rebuilds
the moment equations, verifies their elimination identities, and checks
the exceptional branches, squarefree hyperelliptic polynomials, smooth
plane quartics, and the six simple branch places of each genus-four
cover. The finite-fiber arguments above are part of the mathematical
proof. No search bound, approximate genus computation, GRH assumption,
or conjecture about a Tate–Shafarevich group enters Theorem 14.

The conclusion changes the scope of the remaining work: the entire
symmetric-size-ten shift route is a finite collection of fixed
rational-point problems, with explicitly necessary lifting conditions.
The genus bounds and Faltings do not replace those point calculations.
Section 5.31 separately excludes type 10; the theorem here does not
claim exclusions for the other seventeen types.


### 5.33 Exact bounded center search in four remaining types

For types 4, 5, 15 and 16, write the last squared center as
`c=(n/d)^2`, with `n>=0`, `d>0`, and `gcd(n,d)=1`. An exact search over
`n<=1000000`, `d<=1000000` leaves only the following configurations
in which both center quadratics split completely into rational squares:

| Type | n/d | First center pair | Second center pair | Invalidity |
|---|---|---|---|---|
| 4 | 1/2 | 2, 1 | 5/2, 3/2 | repetitions and cross-color overlap |
| 5 | 7/2 | 3, 1 | 2, 0 | repetitions and cross-color overlap |
| 5 | 9/2 | 4, 2 | 3, 1 | cross-color overlap |
| 15 | 1/2 | 3/2, 1/2 | 5/2, 3/2 | repetitions and cross-color overlap |
| 15 | 7/2 | 9/2, 7/2 | 3/2, 1/2 | cross-color overlap |
| 16 | 3/2 | 7/2, 1/2 | 3, 0 | repetitions and cross-color overlap |

Signs and pair order do not change the reflected chains. Substitution
in the original chain lists verifies all four moment equations, but
none is a disjoint, distinct-entry seed. The exceptional branches were
already excluded in Section 5.32. Thus there is no valid input of these
four types in this **center-height box**. This is not a height bound
on all rational points or on the final integer PTE sets.

Let `S,T` and `U,V` be the two center-quadratic coefficient pairs in
Section 5.32. Put `D=12,16,24,12` for types 4,5,15,16 respectively.
The scaled center squares `X=D*d^2*u^2` are roots of

```
X^2-Sbar X+Tbar,
Sbar=D*d^2*S(n^2/d^2),  Tbar=D^2*d^4*T(n^2/d^2),
```

and likewise for `U,V`. All barred coefficients are integers. Rational
roots of a monic integral polynomial are integers. Therefore both
original centers are rational if and only if

```
r^2=Sbar^2-4Tbar,
Xplus=(Sbar+r)/2, Xminus=(Sbar-r)/2,
D*Xplus and D*Xminus are nonnegative integer squares.
```

Indeed a rational number whose square is an integer is integral, and
the centers are `sqrt(D*Xplus)/(D*d)` and
`sqrt(D*Xminus)/(D*d)`. Integer square-root decisions are exact;
floating point supplies only an initial root guess.

For each of the 32 primes from 7 through 149, the sieve checks whether
each reduced center polynomial has two roots in `{D*r^2:r in F_p}`.
No prime divides `D`. The homogeneous coordinates `(n^2,d^2)` include
the denominator-zero chart. The pair `(n,d)=(0,0) mod p` is discarded
only because it violates primitivity. A Chinese-remainder wheel for
7,11,13,17 enumerates every numerator surviving those four tests.
Subsequent prime tests and the exact rational-square test only discard
impossible configurations. The input height cap of ten million bounds
all signed 128-bit intermediates strictly below overflow.

The completed million-height run visited 60,468,577,849 wheel candidates;
six primitive pairs survived all modular tests, and all six passed the
exact center-square test. The independent Python checker rebuilds the
four even-moment identities, checks the actual C++ coefficient tables,
compares a direct rational-arithmetic enumeration through height 100,
and reconstructs the original chains to certify every collision above.
This is a finite computation, not an effective replacement for Faltings
or an unbounded rational-point classification.

Reproduction:

```sh
mkdir -p build
g++ -O3 -std=c++17 -fopenmp scripts/search_four_shift_types.cpp -o build/search_four_shift_types
build/search_four_shift_types 1000000 12 > build/four_types_1000000.txt
python scripts/check_four_shift_search.py
python scripts/check_four_shift_search.py build/four_types_1000000.txt
```

The optional last two C++ arguments select an inclusive denominator
interval without changing the numerator bound. Only a run containing
a `COMPLETE` line counts as finished; the whole-box checker also
requires denominators 1 through the stated height. Binaries and search
output stay in ignored `build/`.

### 5.34 A rank obstruction on the type-16 necessary curve

The genus-two curve R16 from Section 5.32 is

    Y²=f(c)=6c(c−6)(4c−9)(160c²−1335c+2304).

Its Jacobian has rational Mordell–Weil rank at least two. Consequently
the classical Chabauty criterion `rank J(Q)<genus` cannot be applied
directly to this necessary curve.

**Proof.** The curve has rational points

    P=(3/2,567/2),  Q=(48/5,127008/25),

as direct substitution verifies. Put h=3840, X=hc, Z=h²Y. The new
equation is Z²=F(X), with F(X)=h⁴f(X/h) monic of degree five. It has
three rational roots h·0, h·6, h·9/4 and one irreducible quadratic
factor: the latter's discriminant before scaling is307665, not a square.
The rational two-torsion therefore has dimension three, generated by
the three branch classes T_r=[(hr,0)−∞]. This also follows directly
from the Galois action on even subsets of the six branch points.

Use the injective odd-degree two-Kummer map from Section 5.13 into
the square classes of Q[X]/F. For a nonbranch point with coordinate hc,
its three rational components are h(c−r), r=0,6,9/4. For T_r, the
component at r is F′(hr)=h³f′(r), and the other components are
h(r−r′). Thus their square classes, projected to the three rational
factors, are exactly

| Class | Component at 0 | Component at 6 | Component at 9/4 |
|---|---:|---:|---:|
| T₀ | 15 | −10 | −15 |
| T₆ | 10 | 6 | 1 |
| T₉/₄ | 15 | −1 | −6 |
| [P−∞] | 10 | −30 | −5 |
| [Q−∞] | 1 | 6 | 1 |

Apply, in order, the five square-class characters
`v₂(component 3), v₃(component 1), v₅(component 2),
v₂(component 1), v₃(component 2)`, all modulo two. The resulting
five-by-five matrix, in the displayed row order, is

    [0 1 1 0 0]
    [0 0 0 1 1]
    [1 1 0 0 0]
    [0 0 1 1 1]
    [0 0 0 0 1]

and has odd determinant. Hence these five classes are independent in
J(Q)/2J(Q), even after discarding the quadratic component. Since
dim J(Q)[2]=3, the Mordell–Weil theorem gives rank J(Q)≥5−3=2.
This proves the claim.

The exact checker is `python scripts/check_shift_curve_rank_obstructions.py`.
This is a rank lower bound, not a rank upper bound or a saturated basis.
The two displayed abscissas are not rational squares, so these points
do not themselves lift to rational PTE centers. Covering methods,
higher-depth methods or other rational-point arguments remain possible;
the obstruction only rules out a direct application of the classical
rank-smaller-than-genus criterion to R16.

### 5.35 Higher rank obstructions for types 5 and 15

The other necessary genus-two curves in the bounded search satisfy

    rank J(R4)(Q)≥1, rank J(R5)(Q)≥3, rank J(R15)(Q)≥4.

Thus the direct classical Chabauty rank criterion fails for R5 and
R15 as well as R16. The R4 lower bound does not decide its applicability.

**Proof for R5 and R15.** Use the monic change X=hc, Z=h²Y and
the rational components of the odd-degree Kummer map from Section5.34.
Both curves have three rational linear factors and one irreducible
quadratic; consequently their rational two-torsion dimension is three.

For R5, h=1280 and the component order is r=0,49/4,9/4. Take
P₁=(1/4,192), P₂=(25/4,1920), P₃=(81/4,27648).
Their projected square classes, together with the three branch classes,
are

| Class | Component 1 | Component 2 | Component 3 |
|---|---:|---:|---:|
| T₀ | 2945 | −5 | −5 |
| T₄₉/₄ | 5 | 2 | 2 |
| T₉/₄ | 5 | −2 | 7 |
| [P₁−∞] | 5 | −15 | −10 |
| [P₂−∞] | 5 | −30 | 5 |
| [P₃−∞] | 5 | 10 | 10 |

For R15, h=2560 and the component order is r=0,23/2,31/4. Take
P₁=(1/4,135), P₂=(9/4,1221), P₃=(25/4,2205),
P₄=(2/5,4662/25). The projected square classes are

| Class | Component 1 | Component 2 | Component 3 |
|---|---:|---:|---:|
| T₀ | 3565 | −115 | −310 |
| T₂₃/₂ | 115 | 17871 | 6 |
| T₃₁/₄ | 310 | −6 | −341 |
| [P₁−∞] | 10 | −2 | −3 |
| [P₂−∞] | 10 | −370 | −55 |
| [P₃−∞] | 10 | −210 | −15 |
| [P₄−∞] | 1 | −111 | −6 |

All displayed points and square classes follow by exact substitution.
Write (j,p) for v_p on component j modulo two, and (j,−1) for its
sign bit. Apply the following ordered columns to the respective tables:

    R5:  (1,5),(1,19),(2,−1),(2,2),(2,3),(2,5);
    R15: (1,2),(1,5),(1,31),(2,−1),(2,3),(2,5),(3,2).

Both square matrices have determinant1 over the integers, hence are
invertible modulo two. The tables therefore give six and seven
independent classes in J(Q)/2J(Q), respectively. Subtracting the
rational two-torsion dimension3 proves the rank lower bounds3 and4.

**The R4 lower bound.** Here h=184320 and f has the rational root0
and two irreducible quadratics, so dim J(Q)[2]=2. Use T₀, the
rational branch-divisor class T_q for q(c)=240c²+40c+167, and
P=(1/4,270). The class T_q is the sum of the two branch points
minus2∞; its double is the divisor of the corresponding monic quadratic.
In the component at the rational root0 their Kummer values are

    h³f′(0),  h²·167/240,  h/4,

with square classes1167330,2505,5. The characters v₂₃₃,v₁₆₇,v₅
give an odd-determinant three-by-three matrix, proving independence
in J(Q)/2J(Q). Hence rank≥3−2=1.

The exact checker is `python scripts/check_shift_curve_rank_obstructions.py`.
It rebuilds the monic models, factor counts, points, all Kummer components
used above and the invertible minors. These are rank lower bounds,
not upper bounds, saturated bases or rational-point classifications.
Several displayed points lie on excluded exceptional or colliding
configurations. None supplies a valid new size-eleven output.

### 5.36 Rational Néron–Severi rank and the quadratic-Chabauty criterion

All four Jacobians J(R4), J(R5), J(R15), J(R16) have rational
Néron–Severi rank exactly one. Consequently the standard sufficient
quadratic-Chabauty inequality does not apply directly to R5, R15 or
R16 either. This does not assert that their depth-two sets are infinite
or that every variant of nonabelian Chabauty fails.

**Proof.** The following exact good-reduction counts determine the
Frobenius polynomials P_p(T)=T⁴−a₁T³+a₂T²−pa₁T+p². Here
a₁=p+1−#C(F_p), a₂=(#C(F_{p²})−p²−1+a₁²)/2.

| Curve | p | #C(F_p) | #C(F_{p²}) | a₁ | a₂ | Real fixed field |
|---|---:|---:|---:|---:|---:|---|
| R4 | 17 | 18 | 238 | 0 | −26 | Q(√15) |
| R4 | 23 | 34 | 538 | −10 | 54 | Q(√17) |
| R5 | 11 | 12 | 118 | 0 | −2 | Q(√6) |
| R5 | 13 | 16 | 154 | −2 | −6 | Q(√33) |
| R15 | 17 | 20 | 290 | −2 | 2 | Q(√33) |
| R15 | 41 | 38 | 1614 | 4 | −26 | Q(√7) |
| R16 | 11 | 12 | 118 | 0 | −2 | Q(√6) |
| R16 | 13 | 20 | 170 | −6 | 18 | Q(√17) |

Each P_p is irreducible over Q and p does not divide a₂, so each
reduction is a simple ordinary abelian surface over F_p. The same
specialization argument used in Section5.14 applies: rational
endomorphisms inject into the quartic field Q(π_p)=End⁰_{F_p}(J_p),
compatibly with the Rosati involution π_p↦p/π_p. Its fixed field is
generated by u=π_p+p/π_p, satisfying

    u²−a₁u+(a₂−2p)=0.

The discriminants a₁²−4(a₂−2p), in table order, are
240,68,96,132,132,448,96,68, giving the displayed quadratic fields.
For each curve the two real quadratic fields are nonisomorphic.

Since End⁰_Q(J) embeds in a field, it is commutative and has no zero
divisors, hence is itself a field. Its Rosati-fixed subfield embeds in
both displayed quadratic fields and must be Q. The identification
of NS(J)_Q with Rosati-symmetric rational endomorphisms therefore
gives ρ_Q=1. This is rational, not geometric, Néron–Severi rank.

For all four curves the usual inequality is consequently

    rank J(Q) < g+ρ_Q−1 = 2.

The lower bounds in Sections5.34–5.35 violate it for R5, R15 and
R16. For R4 its applicability still requires a rank-one certificate;
the interval in Section5.37 does not supply one.
The checker `scripts/check_shift_curve_rank_obstructions.py` independently
counts over F_p and F_{p²}, checks good reduction, ordinarity and exact
irreducibility, and verifies the distinct fixed fields. No approximate
Frobenius roots or numerical endomorphism estimates enter this proof.

### 5.37 Software-assisted Selmer upper bounds and exact ranks

An exact two-descent computation gives the following bounds and ranks:

| Curve | dim Sel²(J/Q) | dim J(Q)[2] | Independently proved lower bound | Rank conclusion |
|---|---:|---:|---:|---|
| R4 | 5 | 2 | 1 | 1≤rank≤3 |
| R5 | 6 | 3 | 3 | rank=3 |
| R15 | 7 | 3 | 4 | rank=4 |
| R16 | 5 | 3 | 2 | rank=2 |

The upper bounds are software-assisted: `scripts/check_shift_curve_selmer.m`
uses Magma's documented odd-degree two-Selmer implementation, explicitly
selecting `Al:="TwoSelmerGroupOld"`. It was executed successfully in
Magma2.30-1 through the official calculator. The returned group orders
are32,64,128,32 and the rational two-torsion orders are4,8,8,8. The
script asserts every value. The Python checker supplies the independent
lower bounds, not an independent implementation of the full descent.

No GRH assumption is used in these upper bounds. The script explicitly
sets `SetClassGroupBounds(1000)`. Every nontrivial field in the descent
algebra Q[c]/f(c) is quadratic, with fundamental discriminants

    R4: −231,−155; R5: 505; R15: −215; R16: 34185.

Their unconditional Minkowski bounds are below1000 (indeed below93).
For a real quadratic field the bound is √D/2; for an imaginary one it
is2√|D|/π. The script checks the actual maximal-order discriminants
and Minkowski bounds before running descent, and the Python checker
independently reconstructs the same fundamental discriminants from
the quadratic factors. Thus the imposed class-generator bound is
unconditional. Magma's printed suggestion to use GRH for performance
does not change this setting or enable that assumption.

The exact Kummer sequence

    0 → J(Q)/2J(Q) → Sel²(J/Q) → Sha(J/Q)[2] → 0

gives rank≤dim Sel²−dim J(Q)[2]. Combining with Sections5.34–5.35
proves the displayed exact ranks for R5,R15,R16; in these three cases
it also gives Sha(J/Q)[2]=0. No finiteness of the full Tate–Shafarevich
group, saturated Mordell–Weil basis or rational-point completeness
is inferred. R4 remains a rank interval, not a rank-one certificate.

The software interfaces and their correctness/rigour conventions are
documented in the primary Magma handbook:
[two-Selmer groups](https://docs.magma-maths.org/ArithmeticGeometry/HyperellipticCurves/jacobians-number-fields.html)
and [class-group bounds](https://docs.magma-maths.org/GlobalFields/NumberFieldsAndOrders/ideal-class-group.html).
The reproducible Magma program is a separate optional check; no Magma
installation is required to run the Python and C++ certificates above.

### 5.38 The closed two-vertex Richelot graph for type 4

Let C be the type-4 necessary curve

    C: y²=3x(240x²+40x+167)(256x²+97x+466),

and put

    D: y²=−6477(128x²−233)(240x²−167)
                  (13040x²+138176x+2441).

**Proposition.** Over Q, each of J(C) and J(D) has exactly one
Richelot kernel. Its quotient, with the induced principal polarization,
is the other Jacobian. Neither vertex has an indecomposable double
Richelot isogeny. The two vertices are not isomorphic as principally
polarized abelian varieties. Consequently their connected component
in the rational Richelot/double-Richelot graph consists of these two
vertices and their mutual edge. In particular no further Jacobian is
obtained by a polarization-compatible isogeny of 2-power degree.

Here compatibility means that pulling back the target principal
polarization gives a positive integer multiple of the source principal
polarization. We use the factorization theorem for such 2-power
isogenies into rational Richelot and double Richelot isogenies, as
stated in the [Magma handbook's Richelot section](https://magma.maths.usyd.edu.au/magma/handbook/text/1615).
Its graph-traversal caveat about encountering non-Jacobian surfaces
does not apply: the calculation below exhausts the outgoing kernels
at both vertices and none has a non-Jacobian quotient.

**Proof.** Write A=(x,240x²+40x+167,256x²+97x+466), and let Δ be
the determinant of the coefficient rows in the order (x²,x,1).
For [u,v]=u'v−uv', the Richelot dual of y²=cA₁A₂A₃ has polynomial

    −(c/Δ)[A₂,A₃][A₃,A₁][A₁,A₂].

Direct calculation gives Δ=−69088 and the three brackets

    13040x²+138176x+2441,  2(128x²−233),  167−240x².

The displayed polynomial for D is 8636² times this dual polynomial.
Conversely, take B=(128x²−233,240x²−167,
13040x²+138176x+2441). Its determinant is −4773151744 and its
brackets are

    138176(240x²+40x+167),
    −69088(256x²+97x+466),
    69088x.

With c=−6477 the resulting dual polynomial is 17272² times the
polynomial for C. Both determinants are nonzero and both curve
polynomials are square-free. These are rational Jacobian quotients,
with the indicated square factors accounting for their twists.

Richelot kernels correspond to Galois-invariant partitions of the six
Weierstrass points into three unordered pairs. On C, the points 0 and
∞ are rational, while the two quadratic factors have discriminant
square classes −155 and −231. These classes are independent in
Q*/Q*². The Galois action is therefore generated by the independent
transpositions within those two pairs. In an invariant matching, a
transposition fixing the other four points forces its two exchanged
points to be paired together: otherwise the partner of a fixed point
would change. Thus both quadratic pairs are forced, and 0 must be
paired with ∞. There is exactly one invariant matching.

For D the three quadratic discriminant square classes are

    466, 2505, 32925841=13·19·133303.

They are independent: the prime 233 occurs only in the first class,
167 only in the second, and 13 only in the third. The Galois group
acts by three independent transpositions, so the same argument forces
each of the three displayed quadratic pairs. Again the matching is
unique. The faithful action of S₆ on genus-two two-torsion identifies
the two-torsion Galois image orders as 4 for J(C) and 8 for J(D).
They cannot be isomorphic over Q, even before imposing polarization.

It remains to exclude a rational double Richelot isogeny whose two
individual steps are not rational. Its kernel K has abstract type
Z/4 × Z/2 × Z/2. The annihilator of K∩J[2] in J[2] is a rational
line generated by P. There are three maximal isotropic planes
containing P, each eligible as a first Richelot kernel. If no first
step is rational, Galois must permute these three planes without a
fixed point, hence transitively. At both vertices the action on J[2]
factors through a 2-group. Its orbit lengths are powers of two, so an
action on three objects has a fixed point. This contradiction excludes
indecomposable double Richelot edges. The factorization theorem now
shows that every compatible 2-power route stays at these two
vertices. ∎

The exact checker `scripts/check_type4_richelot_graph.py` verifies both
determinants, square factors, square-freeness and discriminant classes.
It also enumerates all fifteen matchings, independently checking the
unique fixed matching at each vertex. The optional Magma checker
`scripts/check_type4_richelot_graph.m` verifies both directed neighbors,
the absence of double edges, and the complete `TwoPowerIsogenies`
output. It passed in Magma2.30-1; no class-group assumption, rank
algorithm or point-search bound enters this result.

This closes the search for additional vertices in this rational
2-power graph. It does not determine the rank of J(C), which remains
between one and three, or classify C(Q). Higher descent on either of
the two vertices and odd-degree isogenies are not excluded. In
particular, finding this dual does not by itself improve Section5.37's
Selmer bound or settle the type-4 lifting conditions.

## 6. Open boundary

The original distinct-entry size-11 problem remains unsolved here: no pair
of disjoint eleven-element integer sets and no general impossibility proof
has been established. Section 5 gives a distinct-entry shift-lift reduction and
normalized finiteness for all eighteen chain types (Theorem 14), and
Section 5.31 excludes valid distinct-entry lifts in the older sixth mode
(type 10 in the complete table). The remaining seventeen types' required
rational-point sets and lifting conditions are not classified. Section 3
restricts symmetric distinct-entry constructions. Section 4 applies to
the unrestricted multiset formulation, and proves at least fourteen distinct values across both supports,
or at least eight per side under symmetry. The remaining splitting and
higher-genus torsion and rational splitting problems in Sections 4.5–4.10
are unresolved here. Theorem 12 excludes two exceptional parameter loci
and gives an exact quintic formulation. Theorem 13 gives an exact
simultaneous-splitting Pell criterion for total support fourteen. Neither
result increases the support bound or establishes existence.

## Appendix: a splitting pencil and its 2-Selmer group

[PENCIL_SELMER.md](PENCIL_SELMER.md) studies the distinct-entry construction
P(x)(Q(x)+t), with P=x⁵−110x³+2629x−2520 and
Q=x⁶−110x⁴+2629x²+2520x. Its sextic discriminant is 64f(t),
giving a genus-two necessary rational-point condition. A complete staged
2-descent, with an exact unconditional class-number and global-basis certificate, gives
Selmer dimension three and 1≤rank J(Q)≤3. The torsion-free conclusion
and rank lower bound one are unconditional. Exact certificates, the
local-image completeness argument and reproducibility instructions are
in the linked document and scripts/pencil_selmer/. This does not prove
rank one, classify rational points, or settle splitting of a valid fiber.

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

[3] B. Mazur, *Modular curves and the Eisenstein ideal*, Publications
Mathématiques de l'IHÉS 47 (1977), 33–186, Introduction, Theorem (8).
[DOI and journal record](https://doi.org/10.1007/BF02684339).

[4] Z. Scherr and K. Thompson, *Quartic Integral Polynomial Pell Equations*,
arXiv:2307.04566, Section 2.
[Record](https://arxiv.org/abs/2307.04566),
[full text](https://arxiv.org/html/2307.04566v1).

[5] J. S. Balakrishnan, F. Najman, A. Shnidman, and A. V. Sutherland,
*Rational torsion on simple genus two Jacobians*, arXiv:2608.28543,
Introduction, Section 1.2, and Section 4.4, equation (10).
[Record](https://arxiv.org/abs/2608.28543),
[full text](https://arxiv.org/html/2608.28543v1).


[6] G. Faltings, *Endlichkeitssätze für abelsche Varietäten über Zahlkörpern*,
Inventiones Mathematicae 73 (1983), 349–366, Theorem 7;
erratum, Inventiones Mathematicae 75 (1984), 381.
[DOI](https://doi.org/10.1007/BF01388432),
[original paper](https://pazuki.perso.math.cnrs.fr/index_fichiers/Faltings83.pdf),
[erratum](https://doi.org/10.1007/BF01388572).
