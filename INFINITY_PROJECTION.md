# Certified infinity projection for the sixth-mode curves

This is a local arithmetic input to PR #1, not a rational-point classification or a solution of the original distinct-integer PTE problem. The moving QC functions, root isolation, and final sieve remain unfinished.

## Models and normalization

Use the models and compatible complement of PROOF.md §§5.28–5.29. Explicitly
F=(N²+5D²)(N²+275D²), W²=F, Y=W/c, f=F/c², where
- δ=1: N=80−t², D=4t−30, M=2t²−30t+160, c=1;
- δ=5: N=25t²−120t+80, D=10t²−8, M=−60t²+100t−48, c=5625.

Set ωj=t^j dt/(2Y), ρ=ω3. If rj is the coefficient of u^j in (u⁸f(1/u))^(−1/2), use the residue-zero basis
η=(ω0,ω1,ω2,ω4−r1ρ,ω5−r2ρ,ω6−r3ρ).
The triples (r1,r2,r3) are
(0,−2080,33600) and (208/75,6752/1875,−386816/84375).

Let H=∞+ +∞−. Our third-kind convention is νQ with residue divisor Q−H/2.
Thus ν∞ε=−ερ. It is not the differential for ∞+−∞−, which would be −2ρ.

## Certified result

All vectors below are modulo 13³=2197:
- δ=1: Ψρ=(1433,392,1539,1586,728,1248);
- δ=5: Ψρ=(1194,1844,368,1690,533,1898).

Projection onto the holomorphic span (ω0,ω1,ω2) along the prescribed compatible complement gives
- δ=1: λρ=(289,925,590);
- δ=5: λρ=(427,1168,628).

Consequently λν∞+ is respectively (1908,1272,1607) and (1770,1029,1569).
These constants have been reproduced by independent implementations, not inferred from a finite sample of curve points.

## Frobenius computation and precision proof

Write R=f(t¹³)−f(t)¹³. The pullback of ωi is expanded as
13 t^(13(i+1)−1) Σ(n≥0) binom(−1/2,n) R^n dt/(2Y^(13+26n)).
The script retains n=0,…,4 and performs polynomial vertical and horizontal reduction.
It uses modulus 13²⁵ and an auxiliary numerator scale 13⁸.

For retained terms, at most four vertical denominator divisions (13,39,65,91) and three horizontal divisions (26,52,78) lose precision. Hence at least ten digits survive before truncating the analytic tail.

Here is the separate infinite-tail justification. For a differential
A(t)dt/(2Y^(2m+1)), finite branch principal parts can be removed in one odd primitive
B(t)/Y^(2m−1), deg B<8m. Since f has good reduction at 13, the finite étale interpolation and inverse power-series substitutions are integral. Integrating local principal parts loses at most floor(log₁₃(2m−1)) digits. The remaining polynomial numerator has degree at most max(deg A−8m,6).

At infinity, integrate its Laurent principal part and take the polynomial part of the primitive divided by Y. This gives the sharper horizontal loss floor(log₁₃(max(1,d−3))) for numerator degree d; one must not add the valuations of every elimination pivot to this sharp bound.

For term n, m=6+13n and deg A≤90+103n, since deg R≤103. Thus d≤max(42−n,6). For n≥5 the horizontal loss is at most one. The term valuation is therefore at least
n+1−floor(log₁₃(11+26n))−1 ≥4.
It tends to infinity with n. The whole omitted tail is O(13⁴).

In the affine basis (η,ρ), the computed Frobenius has residue row (0,0,0,0,0,0,13). Write its proper block as FdR and its final proper column as v. Frobenius compatibility of Ψ gives
(FdR−13I) Ψρ=v.
This is the even-degree projection formula in Gajović–Müller, Proposition 5.3:
https://arxiv.org/html/2307.15787

For both curves the determinant has valuation 3, the inverse entries have minimum valuation −1, and the solved vector is integral. The certified matrix modulo 13⁴ thus determines Ψρ modulo 13³. Stability of inversion follows since the approximate inverse times the matrix error is O(13³).

Additional checks: the affine characteristic polynomial is
(T−13)(T²−2T+13)(T⁴+6T²+169) modulo 13⁴.

## Exact complement verification

With κ=(N'D−ND')/M, reduce the following coefficient functions multiplying dt/W:
κMN²/D²,
κ(N³+15N²D/4)/D²,
κ(N⁴+4125N²D²/32+525ND³)/D³.
These are the previously fixed elliptic and genus-two complement classes.

Exact rational Hermite reduction gives residue zero. The determinants of the lower 3×3 blocks in the η basis are respectively 3/760 and 4271484375/9728, both 13-adic units. The cup matrix, computed from Laurent primitives at both infinities, has determinant 1/2304; the three complement classes pair to zero exactly. Thus this is an integral Lagrangian complement, and projection loses no additional precision.

The executable certificate is scripts/check_infinity_projection.py. It recomputes Frobenius, checks the affine residue row and characteristic polynomial, solves the projection, and performs the exact complement reduction.

## Uniform kernel near a source at infinity

For a moving source Q=(a,b) in the original W-model, put h=1/a and V=h⁴b. Instead of the finite-source kernel use
ν∞(h,V)=−h dt/(2(1−ht))−V t³dt/(2W(1−ht)).
It has the same residue divisor as νQ and differs from it by a holomorphic polynomial differential. Their normalized third-kind differentials therefore agree exactly. At h=0, V=εc it becomes −ερ.

This removes the artificial holomorphic blowup when σP tends to infinity. It does not by itself calculate the full P-dependent de Rham projection or certify any zero of the QC functions.
