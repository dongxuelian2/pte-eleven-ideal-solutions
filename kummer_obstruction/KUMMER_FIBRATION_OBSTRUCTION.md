# A structural obstruction to elliptic-fibration shortcuts

## Statement and scope

Let J be the genus-two Jacobian in the certified PTE pencil computation. For each of its seven nonzero rational 2-Selmer classes ε, let Xε be the minimal resolution of the associated twisted Kummer quartic (the precise twist, not an arbitrarily chosen geometrically isomorphic quartic).

**Theorem.** Xε admits no genus-one fibration defined over Q. More strongly, no such fibration exists after any finite odd-degree extension of Q, or after any extension of degree less than ten. A section is not assumed in this statement.

The degree bound is sharp: for each ε there is an extension L/Q of degree exactly ten over which Xε admits a genus-one fibration over P¹_L. Thus ten is the minimum degree even when a split base is required. A section is not asserted.

This is an unconditional structural obstruction, not a bounded point search. It does not exclude rational points on Xε, does not compute a Cassels–Tate pairing, does not determine the Jacobian rank, and does not solve the original size-eleven PTE problem. It rules out a specific proposed way of replacing the current Kummer-point problem by a rational elliptic fibration. The degree-ten construction is an existence proof from a divisor class, not an explicit equation or an effective point-recovery algorithm.

## 1. The curve's geometric endomorphism ring

Put h(A)=A^5−660A^3+70983A−102060. Modulo 11 it is A^5−2. In that quotient Frobenius sends A to 4A; the order of 4 in F11* is five. All five roots therefore have Frobenius orbit length five, and h is irreducible modulo 11 and over Q.

Modulo 59,

h=(A+5)(A−9)(A−4)(A²+8A+23).

The quadratic discriminant is 31 modulo 59, a nonsquare. Both reductions are squarefree. Consequently the Galois group of h contains a 5-cycle and a transposition and is S5. One elementary justification is that a transitive group of prime degree is primitive, and the graph consisting of the conjugates of a transposition has a nontrivial invariant component system unless it is connected. Connected transpositions generate the full symmetric group.

For K=Q(a), h(a)=0, the curve's monic defining quintic F has root

θ=(−165a^4+70983a²−255150a−750870)/2.

The identity F(θ)=0 is verified exactly in the accompanying checker and in the earlier certificate. Since θ is not rational and [K:Q]=5, Q(θ)=K. Hence h and F have the same normal closure and the same S5 action on their five conjugates. This is a field identification; it does not assert that the curves y²=h and y²=F are isomorphic.

Zarhin's theorem, in characteristic zero and degree at least five with Galois group Sn or An, now gives End(J_Qbar)=Z. See Theorem 2.1 of [Hyperelliptic Jacobians without Complex Multiplication](https://arxiv.org/html/math/9909052v2), Mathematical Research Letters 7 (2000), 123–132. Thus NS(J_Qbar) has rank one. The standard Kummer resolution formula gives geometric Picard rank 1+16=17 for Xε.

## 2. Every nonzero twist has a transitive node action

The rational Weierstrass point identifies J[2] with

V={ (v1,...,v5) in F2^5 : v1+...+v5=0 },

with S5 acting by coordinate permutations. V is irreducible under both S5 and A5: a nonzero vector has weight two or four, and A5 is transitive on each of these two weight sets. Weight-two vectors span V, while sums of pairs of distinct weight-four vectors give weight-two vectors.

Moreover H¹(S5,V)=0. Indeed F2^5=F2·(1,1,1,1,1) ⊕ V, since five is odd. Shapiro identifies H¹(S5,F2^5) with H¹(S4,F2), which is one-dimensional. Under this identification the inclusion of the constant module induces restriction H¹(S5,F2)→H¹(S4,F2). The sign homomorphism restricts nontrivially, so this map is an isomorphism. The direct-summand cohomology then forces H¹(S5,V)=0. No semisimplicity assertion in characteristic two is being used.

Let N be the kernel of the mod-two representation G_Q→S5. Inflation–restriction injects H¹(G_Q,V) into H¹(N,V)^{S5}, because the preceding H¹(S5,V) vanishes. N acts trivially on V, so the restriction of any nonzero class ε is a nonzero homomorphism N→V. Its image is an S5-stable nonzero subspace; irreducibility implies that its image is all of V.

For the certified twisted Kummer construction, the sixteen nodes form the affine V-torsor specified by ε: they are the points above zero on the corresponding 2-covering before quotient and resolution. The rational Weierstrass point provides an untwisted rational theta-characteristic reference. Explicitly, descent acts by x ↦ ε(σ)+σx; the fixed set of inversion is J[2], so the nodes carry precisely this affine action. Thus the Galois image in the affine group contains all translations V and projects onto S5. It is V⋊S5, and in particular acts transitively on the sixteen nodes.

This applies to each of the seven nonzero certified classes. It is not a statement about the zero twist, whose node action has a rational fixed node.

## 3. The invariant intersection form is anisotropic

Let H be the pullback of a hyperplane section of the Q-defined Kummer quartic, and let E1,...,E16 be the exceptional curves over its nodes. In NS(Xε,Qbar)⊗Q,

H²=4, H·Ei=0, Ei·Ej=0 for i≠j, Ei²=−2.

These seventeen classes are independent by their nondegenerate intersection matrix. Since the geometric Picard rank is seventeen, they span NS(Xε,Qbar)⊗Q.

Galois fixes H and transitively permutes the Ei. Therefore

(NS(Xε,Qbar)⊗Q)^{G_Q}=QH ⊕ QE, where E=E1+...+E16,

and its intersection matrix is diag(4,−32). We deliberately use Galois invariants rather than assume every invariant class descends as a line bundle; no Brauer descent issue is needed for the obstruction.

If aH+bE has square zero, then a²=8b². For rational a,b this forces a=b=0. A genus-one fibration defined over Q would supply a nonzero Galois-invariant geometric fiber class with square zero, a contradiction. This excludes genus-one fibrations even without a rational section or with a nonsplit genus-zero base.

## 4. Odd-degree extensions and extensions of degree less than ten

Let G=V⋊S5 be the image on the nodes, and let L/Q be finite. Its restricted image H_L⊂G has index dividing [L:Q].

If [L:Q] is odd, normality of V gives

[V:V∩H_L]=[VH_L:H_L],

which is both a power of two and a divisor of the odd number [G:H_L]. Hence V⊂H_L.

More generally, suppose [G:H_L]<10. Put U=H_L∩V and let P be the image of H_L in S5. Then

    [G:H_L]=[S5:P] [V:U].

If U were proper, [V:U]≥2, so [S5:P]<5. The only possibilities are P=S5 and P=A5. Indeed an index-three or index-four subgroup would give a transitive homomorphism S5→S3 or S5→S4. Simplicity of A5 and 60>|S4| force A5 into its kernel. The image then has order at most two, which cannot be transitive on three or four letters; the unique index-two subgroup is A5.

The subspace U is P-stable. Irreducibility under S5 and A5 forces U=0, giving [V:U]=16, a contradiction. Thus V⊂H_L. In particular this holds whenever [L:Q]<10, because the restricted-image index divides the extension degree. This adds degrees four, six and eight to the previous quadratic case.

In both cases the nodes remain transitive, so the invariant intersection form remains diag(4,−32). The same isotropic-class argument rules out a genus-one fibration over L.

The argument uses only a rational Weierstrass point, quintic Galois group S5, geometric Picard rank seventeen, and a nonzero true translation-twist class. It therefore applies more generally to such genus-two Jacobians. The Galois-module and Picard-space framework is standard; see Skorobogatov–Zarhin, [Kummer varieties and their Brauer groups](https://arxiv.org/abs/1612.05993), Sections 2 and 4. The low-degree conclusion is an elementary subgroup-index refinement of that framework. It uses invariance of the geometric hyperplane class, not descent of every invariant class to a line bundle, so no Brauer splitting assumption is inserted. A nonzero Mordell–Weil class modulo two is included: its underlying Jacobian torsor may be trivial, but its node J[2]-torsor is nontrivial.

## 5. A degree-ten extension attains the bound

Write the affine node group as G=V⋊S5 as above. Choose the fifth coordinate and put

    U={v∈V : v5=0},    P=Stab_S5(5)≅S4,    B=U⋊P.

The subspace U has dimension three, so |B|=8·24=192 and [G:B]=1920/192=10. Let M/Q be the finite Galois extension cut out by the node action and put L=M^B. Thus [L:Q]=10. On the sixteen nodes B has exactly two orbits

    O0={v∈V : v5=0},    O1={v∈V : v5=1},

each of size eight: translations by U are already transitive on either set. Consequently the exceptional divisor Z=Σ_(v∈O0) E_v is defined over L. Since the given quartic is defined over Q, its hyperplane pullback H is an actual divisor class over Q. The divisor class

    D=2H−Z ∈ Pic(Xε,L)

satisfies

    D²=4H²+Z²=16−16=0,    D·H=8.

It is therefore nonzero and isotropic. No assertion that H−Z/2 is integral is needed.

We now use [Hassett–Tschinkel, *Rational points on K3 surfaces and derived equivalence*, Lemma 15, Section 2.1](https://www.math.brown.edu/bhassett/papers/derivedsurvey/DerivedSurvey17.pdf). Over a characteristic-zero field, a K3 surface with a nonzero square-zero divisor admits a genus-one fibration, allowing a nonsplit conic as base. The proof moves the divisor to the nef cone by Picard–Lefschetz reflections; the unique nef representative is Galois invariant, and a multiple descends and induces the fibration. Thus an isotropic vector alone is not being silently identified with an already nef fiber. Applying this result to D proves existence over L. Together with Section 4, it proves that the minimum possible extension degree is exactly ten.

This argument is not restricted to the seven computed classes. Let C/Q have a rational Weierstrass point and a separable quintic model with Galois group S5. For every nonzero true class ε∈H¹(Q,Jac(C)[2]), the standard translation-twisted Kummer has the same minimum degree ten. Zarhin's theorem supplies Picard rank seventeen, and Section 2 supplies the full affine node group. No Selmer local-solubility condition is required. If only the geometric hyperplane class is specified, rather than a Q-defined quartic model, choose a positive integer n for which nH is defined over Q: take the norm of a line bundle representing H from a finite Galois field of definition, using invariance of H and Pic⁰=0. Then replace D by 2nH−nZ; its square remains zero and its intersection with H is 8n. Thus no Brauer splitting assumption is required in this general form either.

For an arbitrary nonzero H¹ class, this establishes the exact minimum degree with a genus-zero base allowed. For a true **2-Selmer class**, including each of the seven classes here, its degree-ten base is in fact split. Let Yε be the associated 2-covering of J. At every place v of Q, the Selmer condition gives Yε(Q_v)≠∅. The covering is smooth of dimension two, so a local point has a positive-dimensional real or p-adic analytic neighborhood. One can therefore choose a point outside the finite fixed locus of its involution. Its image in the smooth part of the quotient gives a point of Xε(Q_v).

For every place w of L above v, base change along Q_v⊂L_w then gives Xε(L_w)≠∅. If Xε,L→B is the genus-one fibration just constructed, these points give B(L_w)≠∅ at every finite and infinite place w. The smooth projective genus-zero curve B is a conic over the number field L. By the Hasse–Minkowski theorem it has an L-point and hence B≅P¹_L; see the [number-field conic local-global principle](https://magma.maths.usyd.edu.au/magma/handbook/text/1537). All completions, including real places, have been included; no stronger hypothesis on L is needed. This proves the split-base assertion for the seven Selmer classes. More generally the same argument applies whenever the covering is everywhere locally soluble.

Neither the local points nor the splitting of the base give a section of the fibration, an L-point of Xε, or a point on its generic fiber. For arbitrary H¹ classes without local solubility, no split-base assertion is made. The ingredients are standard Kummer geometry, the cited K3 divisor criterion and the conic Hasse principle; the contribution here is their sharp subgroup-and-divisor consequence, not a new general fibration criterion.

## Reproducibility and review

`check_fibration_obstruction.py` verifies the exact polynomial identities, squarefree Frobenius factorizations, the identification F(θ)=0, the elementary module calculations, and the degree-ten subgroup, orbit and divisor-intersection calculations. The cohomological and geometric steps above use the stated standard theorems and are mathematical arguments, not proof-assistant formalizations. The previously certified model labels ensure these statements apply to the actual seven computed quartics.
