# A structural obstruction to elliptic-fibration shortcuts

## Statement and scope

Let J be the genus-two Jacobian in the certified PTE pencil computation. For each of its seven nonzero rational 2-Selmer classes ε, let Xε be the minimal resolution of the associated twisted Kummer quartic (the precise twist, not an arbitrarily chosen geometrically isomorphic quartic).

**Theorem.** Xε admits no genus-one fibration defined over Q. More strongly, no such fibration exists after any finite odd-degree extension of Q, or after any quadratic extension of Q. A section is not assumed in this statement.

This is an unconditional structural obstruction, not a bounded point search. It does not exclude rational points on Xε, does not compute a Cassels–Tate pairing, does not determine the Jacobian rank, and does not solve the original size-eleven PTE problem. It rules out a specific proposed way of replacing the current Kummer-point problem by a rational elliptic fibration. No claim is made for other even-degree extensions.

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

with S5 acting by coordinate permutations. V is irreducible: a nonzero vector has weight two or four; its S5 orbit spans V in either case.

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

## 4. Odd-degree and quadratic extensions

Let G=V⋊S5 be the image on the nodes, and let L/Q be finite. Its restricted image H_L⊂G has index dividing [L:Q].

If [L:Q] is odd, normality of V gives

[V:V∩H_L]=[VH_L:H_L],

which is both a power of two and a divisor of the odd number [G:H_L]. Hence V⊂H_L.

If [L:Q]=2, the image has index one or two. Since the S5-module V is nontrivial and irreducible, [S5,V]=V. Thus G_ab=S5_ab=C2, and every index-two subgroup contains V. Again V⊂H_L.

In both cases the nodes remain transitive, so the invariant intersection form remains diag(4,−32). The same isotropic-class argument rules out a genus-one fibration over L.

## Reproducibility and review

`check_fibration_obstruction.py` verifies the exact polynomial identities, squarefree Frobenius factorizations, the identification F(θ)=0, and the elementary module/intersection calculations. The cohomological and geometric steps above use the stated standard theorems and are mathematical arguments, not proof-assistant formalizations. The main theorem and extension argument were independently reviewed. The previously certified model labels ensure these statements apply to the actual seven computed quartics.
