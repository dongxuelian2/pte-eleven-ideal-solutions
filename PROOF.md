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

## 4. General multiset restrictions

Unlike Section 3, this section assumes neither symmetry nor distinct entries.
All statements apply to rational multisets and therefore to integer multisets.
They are structural restrictions, not a solution of the existence problem.
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


## 5. Open boundary

The general size-11 existence problem remains open here. Section 3
restricts symmetric distinct-entry constructions. Section 4 applies to
the unrestricted multiset formulation, and proves at least fourteen distinct values across both supports,
or at least eight per side under symmetry. The remaining splitting and
higher-genus torsion problems in Sections 4.5–4.7 are unresolved here.

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
