# Moving normalization and QC root isolation at 13

This continuation uses arithmetic inputs at PR head `39bf9e2b7b19656486898bd17c1d3d526a4fe941`. The new identities, series, root isolation, symmetry, and candidate sieve received independent review; the upstream results listed below remain dependencies rather than being silently re-proved here. The candidate sieve below identifies every simultaneous-center-lift candidate as a colliding point; this addresses only the sixth shift-chain mode, not general PTE.

## Moving projection

Use the monic model Y²=f(t)=F(t)/c² and ωj=t^j dt/(2Y) of INFINITY_PROJECTION.md. The projection λ is holomorphic projection of Coleman–Gross Ψ along the previously certified invariant isotropic complement. For Q=(a,y), νQ=(1+y/Y)dt/(2(t−a)). Its source derivative, modulo the exact differential of −1/(2(t−a))−Y/(2y(t−a)), is

R_a(t)dt/(4yY), R_a=((f′(a)+f′(t))(t−a)−2(f(t)−f(a)))/(t−a)².

R_a is polynomial of degree ≤6 and residue zero. Ψ commutes with this source derivative: pair with each fixed second-kind de Rham test differential using the global-symbol characterization. On a fixed residue disc the Coleman primitives are analytic. Differentiate the symbol in the translated source parameter, commuting coefficient extraction of the residue and analytic differentiation. The resulting symbol is exactly that of ∂aν; nondegeneracy of the de Rham cup pairing proves the asserted equality. This uses more than mere Q13-linearity of Ψ.

For the infinity chart h=1/a,v=h⁴y, v²=revf(h), use the uniform kernel
νI=−h dt/(2(1−ht))−v t³dt/(2Y(1−ht)).
Its derivative, modulo d_t[−t/(2(1−ht))−h³Y/(2v(1−ht))], is RI dt/(4vY), where
RI=[−revf′ t³(1−ht)−2revf t⁴+h³f′(t)(1−ht)+2h⁴f(t)]/(1−ht)².
Again RI is polynomial and residue zero.

The published λρ constants determine λ at infinity. At B, σ*ν∞+ has residue (B−ιB)/2, so νB−σ*ν∞+=½dlog(t−tB)+a holomorphic polynomial differential. The logarithmic term has Ψ zero; the polynomial difference is computed exactly. σ preserves the prescribed complement, hence projection commutes with σ. Negative-sheet values are negatives because the two kernels sum to a logarithmic differential.

All source-derivative coefficients are 13-integral (the integrated projection coefficients can have integration denominators): the curve coefficients, projection matrix and the unit square-root series are integral. For source parameter in 13Z13, integrated terms of order n≥3 have valuation ≥n−v13(n)≥3. Thus constant+linear+quadratic projections modulo13³ suffice. Values and executable exact checks are in `scripts/check_sixth_mode_qc_roots13.py`.

## Local height on B discs

Set q=t(P)−tB and r=t(integration variable)−tB. The opposite-sheet finite kernel has the nonsingular expression
νιP = [(f(t)−f(tP))/(t−tP)]dt/[2Y(t)(Y(t)+Y(P))].
This permits an integral two-variable expansion throughout r,q∈13Z13. The σP source is treated by the infinity kernel, with h(q)=L(tP)/A(tP), v(q)=9025Y(P)/A(tP)^4. Every denominator is a unit on this bidisc.

Integrate the published corrected kernel
 dlog(W/M³)+3νσB−νιB+3νσP−νιP
minus its holomorphic normalizer from r=0 to r=q. The result is
- δ1: H_B(P)=991q+2083q² mod13³;
- δ5: H_B(P)=1627q+81q² mod13³.

The exact formulas for the constant, r-linear and q-linear coefficients of this integrand appear in `scripts/check_sixth_mode_qc_roots13.py`. For the raw kernel, omitted integrated terms of total order n≥3 have valuation ≥n−v13(n)≥3. A projection coefficient may itself contain 1/j, followed by integration in r with denominator i+1. Here n=i+j+1 and the bound is n−v13(j)−v13(i+1)≥3 for n≥3,j≥1. Thus the doubly integrated normalizer has the same required tail bound. These inequalities also tend to infinity, proving convergence after scaling q,r by13. These are local heights, not raw unnormalized third-kind integrals.

## Centered global term

Compute the cup matrix in the η basis by Laurent primitives at both infinities. Its sign is checked directly by
 cup(Ψ(∂aνQ),ωj)=−a^j/(2y(a)).
Consequently the centered Abel–Jacobi logarithms of B are −cup(ΨνB,hol). Holomorphic correction pairs trivially with holomorphic forms, so these follow from the certified Ψρ and σ pullback. In the three quotient coordinates of PROOF §5.27 their residues modulo2197 are
- δ1: (0,1313,832);
- δ5: (1339,1339,1079).
All are divisible by13, as independently required by the prime-to-13 reduction group orders.

For raw quotient logs L=(Le,Lw), put J=diag(144,−A_k/2), using A_-15=[[37,131],[131,135]] and A_-3=[[90,114],[114,36]]. The global quadratic expression is LᵀJL/13. If the log increment is d1q+d2q²+..., its change through the required precision is
[2LBᵀJd1 q+(d1ᵀJd1+2LBᵀJd2)q²+2d1ᵀJd2 q³]/13.
The omitted term involving LB times q³ has valuation at least3. For a product of log increments of degrees a,b with n=a+b≥4, its value after dividing by13 has valuation at least n−v13(a)−v13(b)−1≥3. The cubic cross term is explicitly retained. The estimates tend to infinity and apply to the infinite analytic tails, not just finite polynomials. Errors in LB modulo13³ and J modulo13² likewise give value error≥13³.

Thus, with q=13u and rho=H_B−(global change),
- δ1: rho/13 =149u+136u² modulo169;
- δ5: rho/13 =129u+110u²+117u³ modulo169.
In particular reductions modulo13 are 6u+6u² and12u+6u².

The entire analytic functions rho(tB+13u)/13 have integral coefficients on Z13. Terms of degree≥3 vanish modulo13; this follows from the same integral-kernel bounds and the single division by13 in the global quadratic form, not from a truncated floating-point computation. Hence ordinary simple-root Hensel theory is exhaustive on these B discs.

## Simple roots for each exact target

log13(3)=2093 mod2197, hence log13(3)/13=161 mod169, equal5 mod13. The exact targets are the published rational multiples of this logarithm. The following u residue classes modulo13 give all roots on the B disc:

δ1:
- target0: u=0,12, derivatives6,7;
- target−log(3)/6: u=2,10, derivatives4,9;
- target−2log(3)/3: none.

δ5:
- target0: u=0,11, derivatives12,1;
- targetlog(3)/2: u=2,9, derivatives10,3;
- target2log(3)/3: u=3,8, derivatives9,4;
- target−4log(3)/3: u=5,6, derivatives7,6.

Each displayed residue class contains exactly one root of the exact analytic target equation. No other class can contain a root. This yields four roots on B1's disc and eight on B5's disc, not just candidate residue classes.

## Transport to all eight discs

Write H_B(P)=h13(P−B,DΓ(B,P)), and Z_actual=−ι_*+3σ_* on divisors. Directly,
DΓ(A,P)−DΓ(B,P)=Z_actual(A−B),
DΓ(B,P)−DΓ(B,A)=Z_actual(P−A).
Local height symmetry (the complement is isotropic) and functorial self-adjointness of ι,σ give cancellation of the cross-pairings, hence
H_A(P)=H_B(P)−H_B(A).
This holds first where supports are disjoint and then by the analytic extensions furnished by the normalized kernel.

For τ=ι or σ, invariance of the splitting and commutation with Γ give H_τB(τP)=H_B(P). Choose a generic auxiliary point A, with all supports disjoint. The cocycle and equivariance imply H_B(τP)−H_B(P)=H_B(τA)−H_B(A)=Kτ on the disjoint-support open. Applying τ twice gives2Kτ=0, hence Kτ=0. This proves invariance without evaluating a still-undefined colliding-support height. Transport the analytic B-disc function to each orbit disc to obtain its removable continuation; only afterward conclude H_B(τB)=0. The centered global q is also invariant: τH is linearly equivalent to H, τ commutes with Z, and global height is invariant under τ. Therefore rho is invariant under both involutions.

Each curve's four survivor discs are a free orbit of its B disc under {1,ι,σ,ισ}. They are disjoint modulo13. Transporting the isolated roots yields exactly16+32=48 local QC roots across the eight surviving discs. Eight are the known basepoint orbits (one zero-target root per disc); the other40 still require identification or exclusion.

## Candidate-level elliptic-image sieve

The exact elliptic group is E(Q)=Z·P⊕{0,T}, P=(180,900), T=(0,0). The already proved necessary residues are n={0,1,120} mod240 for epsilon0 and {0,120,239} for epsilon1. A new exact finite-field certificate uses only primes61,101,107. The orders of P are respectively26,13,104.

For each point nP+epsilonT, recover z=Y/(2X). If X=0 or the point is infinity, retain it conservatively. Otherwise require the four values z²+5,4z²+30z+80,z²+36z+185,z²+275 to have the same nonzero quadratic character, retaining every zero-value case conservatively. These are necessary conditions for the simultaneous center squares. Full enumeration of the P-cycle, with torsion translation, gives allowed n modulo each exact P-order. Project these sets modulo gcd(order,3120) and intersect with the original mod240 restriction. The survivors modulo3120 are exactly
- epsilon0: {0,1,1560};
- epsilon1: {0,1560,3119}.

This finite congruence argument is used on the isolated QC candidates, not as an independent unbounded-search exclusion.

On the B1 disc the elliptic reduction is O, forcing epsilon0 under the original mod240 restrictions; delta1 also forces n even. Thus the surviving sieve indices have n=0 mod13. On B5 the elliptic reduction is −P+T, forcing epsilon1; delta5 forces n odd, hence n=12 mod13.

The elliptic quotient pulls dX/(2Y) back to dz/(2b). Its t-derivative at each B is11 mod13. The centered elliptic logs are0 and1339=13·103 mod2197. Since f(B5)=−P+T, log(P)/13=−103=1 mod13. An independent check uses exact rational addition to compute6P: its formal parameter−X(6P)/Y(6P)=78 mod169. The invariant differential is integral with leading coefficient1, so log(6P)=78 mod169 and log(P)/13=1 mod13. This calibration does not depend on the Ψ calculation. Therefore, with q=13u,
- B1: n=11u mod13;
- B5: n=12+11u mod13.

For each curve only u=0 meets the sieve constraint. This excludes every extra QC root on its B disc. The simple zero-target root in u=0 is exactly B, since rho(B)=0 by definition and simple-root uniqueness forbids another root in that subdisc. The original simultaneous center-square conditions are invariant under both ι and σ: σ fixes(z,b), while ι fixes z and changes b's sign. Thus symmetry transports the sieve exclusion to all eight discs.

Conclusion, using the previously certified PR inputs listed below: every sixth-mode rational simultaneous center lift is one of the known four-point orbits on X1 or X5. These map to elliptic indices0 or±1 and already fail the exact collision tests. Hence the sixth-mode distinct-entry construction is impossible. This does not classify all rational points on Xδ (the square conditions were used in sieving), does not exclude the other shift-chain types, and does not settle the original size-eleven PTE problem.

## Dependencies and reproducibility

This is a theorem conditional on the correctness of the cited, previously certified upstream arithmetic statements, in exactly the same sense as the earlier sections depend on their established lemmas. The new certificate does not re-prove all95 earlier commits.

- PROOF §§5.7–5.10: coverage of sixth-mode center lifts, exact elliptic group, mod240 necessary index classes, parity/squareclass relation, and collision tests for the known orbits.
- PROOF §§5.14–5.16: genus-three models, quotient maps, full-rank logarithm hypotheses, correspondence, compatible splitting and relevant saturation. This continuation uses the full elliptic group for the final sieve, not an assumed fully saturated genus-three subgroup.
- PROOF §§5.17,5.22,5.27: elliptic global height coefficient, Richelot transport and compatible complement, and the two genus-two global height matrices.
- PROOF §5.28: the corrected graph divisor and the curve-level quadratic-Chabauty equation including the centered global term and its finite-local-value theorem.
- INFINITY_PROJECTION.md: the two certified Ψρ vectors modulo13³ and projection constants. Their Frobenius computation is not rerun by the new root script.
- COMPLETE_AWAY_TARGETS.md: exact finite target sets, not merely their rounded p-adic residues.
- The survivor-disc filter follows directly from the mod240 classes, parity, order6 of P modulo13, and the displayed degree-eight models. At13, P has order6. The delta1 elliptic images must be O or T; the delta5 images must be P=(11,3) or −P+T=(2,10). Enumerating the smooth degree-eight models gives X1 survivors (t,W)=(1,1),(1,12) and the two infinity branches with W/t⁴=±1. X5 survivors are(3,3),(3,10) and infinity branches W/t⁴=±9. Each curve has12 total F13 points. These four retained discs contain B,ιB,σB,ισB. The quotient map at D=0 is interpreted in the proper infinity chart, mapping to O or T, not discarded.

Run from any directory:

    python scripts/check_sixth_mode_qc_roots13.py
    python scripts/check_sixth_mode_candidate_sieve.py

The first script requires the existing SymPy package and computes the complement, exact derivative reductions, cup pairing, local height and centered global series, complete simple-root residue classes, and Hensel refinements. It asserts that the root module uses the coefficients actually computed, rather than an unrelated hardcoded polynomial. The second uses the Python standard library only, exhausts the finite-field elliptic cycles, verifies the exact group orders and good reduction, and independently calibrates the elliptic formal logarithm. Neither script installs software, accesses the network, or writes generated files.

The mathematical analytic-tail and cocycle arguments above are part of the proof; the executable finite algebra checks do not substitute for those arguments. Root refinements modulo169 are reproducibility outputs, not rational-point guesses.

## Scope

The sixth-mode shift-chain construction has no valid distinct-entry PTE lift. The other17 shift-chain modes and the unrestricted size-eleven PTE existence problem are not settled by this result. The 48 local QC roots are not claimed to be48 rational points, and this result does not give a full classification of X1(Q) or X5(Q).
