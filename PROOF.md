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

These results concern lifting a symmetric ideal size-ten pair, not all
size-eleven pairs. They provide exact reductions and a finiteness result
for four subclasses, but no size-eleven witness or general exclusion.
No novelty claim is made for the lifting or finiteness mechanisms.

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
square, complete-splitting, and noncollision conditions. The other fourteen
chain types and constructions outside this shift-lift route are not settled
by these reductions.
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
The remaining twelve types and other construction routes are not settled
by the reductions in this section; no eleven-distinct-integer solution
or general impossibility theorem is established.

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

#####12 Genus-two Richelot targets: certified torsion and positive rank

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

#### The precise remaining certification

The simpler j=2 duals provide four rational branch points and rational
2-torsion of rank three. This is a concrete input for a genuine two-descent
or Richelot-isogeny Selmer computation. It is not the computation itself.
To use ordinary Chabauty on the selected original pair, the missing steps
are:

1. A certified rank upper bound at most one for each of J(D_1) and J(D_5)
   (equivalently the selected original Jacobians). An exact two-Selmer
   computation, Richelot-isogeny descent, and if needed a Cassels-Tate pairing
   are possible methods. The bounds might instead turn out to exceed one;
   no outcome has been assumed.
2. Explicit sufficient Mordell-Weil data on the original Jacobians, with
   the index and relevant saturation requirements proved for the chosen
   point-completeness method. A non-torsion dual point alone does not supply
   all of this information.
3. A complete rational-point certificate for each original G(2,delta),
   using justified Chabauty-Coleman plus a Mordell-Weil sieve or another
   effective method, rather than a finite search.
4. Testing every remaining z against all center-square conditions, the
   singleton-square conditions and the collision exclusions from the PTE
   reduction.

No Selmer group, rank upper bound, saturated Mordell-Weil group, or complete
rational-point set has been computed here for any of the eight curves.
The original distinct-eleven-integer PTE problem remains unresolved.

## 6. Open boundary

The original distinct-entry size-11 problem remains unsolved here: no pair
of disjoint eleven-element integer sets and no general impossibility proof
has been established. Section 5 gives a distinct-entry shift-lift reduction and
normalized finiteness for six of eighteen chain types; it does not classify
the required genus-two and genus-three rational points. Section 3
restricts symmetric distinct-entry constructions. Section 4 applies to
the unrestricted multiset formulation, and proves at least fourteen distinct values across both supports,
or at least eight per side under symmetry. The remaining splitting and
higher-genus torsion and rational splitting problems in Sections 4.5–4.10
are unresolved here. Theorem 12 excludes two exceptional parameter loci
and gives an exact quintic formulation. Theorem 13 gives an exact
simultaneous-splitting Pell criterion for total support fourteen. Neither
result increases the support bound or establishes existence.

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
