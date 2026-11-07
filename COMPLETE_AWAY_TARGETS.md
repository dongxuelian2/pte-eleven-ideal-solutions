# Complete away-from-13 targets for the sixth-mode curves

This completes the local target calculation left open in PROOF.md §5.29, for the two sixth-mode genus-three curves only. It does not classify their rational points or solve the original PTE problem.

Use the models, basepoints Bδ, H, involutions ι,σ, and Z=1+3σ of §§5.28–5.29. The published correspondence is Γ=Γι−3Γσ, acting as −Z. Its divisor is
DΓ(B,P)=Rι−3Rσ−ιP+3σP−ιB+3σB.
The local function is hℓ(P−B,DΓ). With ρB=h13−global, the QC target is the NEGATIVE of the sum of the away local heights.

## 2-adic vanishing without claiming a full wild model

For X1, work over any finite unramified extension K/Q2. If t=2u is integral, then
F1(2u)/16=[4(20−u²)²+5(4u−15)²][4(20−u²)²+275(4u−15)²]≡3 (mod4).
This cannot be a square: any square root must reduce to 1 modulo2 and its square is then 1 modulo4. Thus all K-points have integral s=1/t, including s=0 at infinity.

Put N∞=80s²−1, D∞=4s−30s², S∞=N∞²+140D∞² and
Y∞=(s⁴W+S∞)/2.
Then
Y∞²−S∞Y∞=−18225D∞⁴/4.
The coefficients are integral, S∞≡1 (mod2), and the right side is 0 modulo2. Over the closed s-unit disc, Hensel lifting splits the equation into two analytic discs, exchanged by ι.

Insert their boundary points in an invariant semistable vertex set. Any skeleton refinement within either disc is a genus-zero hanging tree: it carries neither a cycle contribution nor a positive-genus component contribution to the graph measure for Z. The local graph potential is therefore constant on each disc. Since ι commutes with Z, its pullback preserves the measure; the difference of potentials is constant, and ι²=1 makes this difference zero. Both disc constants agree, and normalization at B makes them zero.

For the boundary-vertex refinement use Baker–Payne–Rabinoff, Proposition 3.13(3):
https://arxiv.org/html/1404.0279

X5 becomes X1 over the unramified quadratic extension Q2(√5), by a1=√5 a5 in the original biquadratic model. This preserves ι,σ,Z. Consequently h2 is identically zero on both curves. This argument does not assert that a complete wild regular model has been computed.

## 3-adic graph and accessible points

For X1 the branch polynomial is
(t⁴−80t²−1200t+10900)(t⁴+4240t²−66000t+253900).
Choose z0²=−5,z3²=−275, both congruent to1 modulo3. Their residues modulo3⁸ are 5731 and3733. The branch roots are −2z±√q1(z), for z=±z0,±z3.

There are two four-root clusters:
A, near −1, depth1/2, with two twins of depth5/2;
B, near4, depth3/2, with one twin of depth2 and two singletons.
Over Q3(√3), their integral depths are1,3 and ν-values4,12. A is übereven; B has genus0. The top cluster gives two unstable rational components which contract to the bridges below.

The semistable cluster-to-graph rules used here are Theorems7.1 and8.6 of Dokchitser–Dokchitser–Maistret–Morgan:
https://link.springer.com/article/10.1007/s00208-021-02319-y

In lengths normalized by v3(3)=1 the graph has:
two edges A+→A− of length4;
bridges A+→B and B→A− of length2;
a loop at B of length1.
σ exchanges the long edges, fixes the bridges, and reverses the loop.
The cycle columns are (1,−1,0,0,0), (1,1,−2,−2,0), (0,0,0,0,1).
Their weighted Gram matrix is diag(8,24,1); Z acts by (−2,4,−2).
The resulting edge densities are (−1/12,−1/12,2/3,2/3,−2).

Normalize the potential j at B1, which retracts to bridge coordinate s=1/2 measured from A. Then
j(A±)=−1/12, j(B)=−3/4,
jbridge(s)=−1/12+s/3−s²/3,
jlong(x)=x²/24−x/6−1/12,
jloop(x)=x²−x−3/4.
The outgoing slopes sum to zero at each vertex.

Integer t-valuations over Q3, or the unramified field Q3(√5) used for X5, leave only A±, bridge coordinates1/2 and3/2, and the B-loop midpoint accessible. There is no entry into an A twin (which would require valuation1/2), and no B vertex (depth3/2); deeper entry into the B twin still retracts to its midpoint.
For X1, F1(4+9u)/3⁷ reduces to 2u²+2u+1, with no F3-root. Thus its valuation is odd and that midpoint is impossible.
The X1 t-values −1,0,1 realize all three remaining potential values.
For X5, the original z-values17,0,16,331 realize all four positions; q1/5 and q0q3 are Q3-squares. The two coefficients in t1=√5a−2z belong to the integral unramified basis, so their minimum valuation gives the retraction without cancellation.
The X5 basepoint has v3(t1−4)=1 and potential −1/3 in the X1 normalization.

Hence the normalized potentials have exact value sets
j1={−1/12,0,−1/3},
j5={1/4,1/3,0,−2/3}.
These are exhaustive local sets, not a finite search of rational points.

## Independent scalar and sign calibration

The graph theorem determines a potential only up to its height normalization; we calibrate it by intersection theory rather than guess a factor.

For X1 choose P with t=−1 and W/3≡−1 modulo3. Use the opposite correspondence C′, whose divisor is −DΓ. Write
D_C′=div(g)−E, g=M³/W,
E=A_P+A_B, A_Q=3σQ−ιQ−H.
Then v3(g(P)/g(B))=2.
The horizontal intersection (P−B)·E is3/2: only 3(P·σP) contributes, since v3(σt−t)=1 and the A-chart coordinate divides by √3.

For the vertical correction, order graph vertices (A+,A−,B,b+,b−), where b± lie at bridge coordinate1/2. The divisor vectors are
d=(1,0,0,−1,0), e=(3,−1,0,1,−3).
The two parallel long edges have effective resistance2; the split bridge lengths are1/2 and3/2. The graph Laplacian solves
Lφ=d, φ=(5/24,1/24,−1/8,−1/4,0),
and φ·e=1/3.
Thus h3(P−B,D_C′)=(2−3/2−1/3)χ3(3)=−log13(3)/6.
Changing to the published DΓ reverses the sign. Since j(P)=−1/12, the published scalar is −2log13(3).

## Exact local and QC targets

The complete local sets are therefore
S3,1={0,1/6,2/3} log13(3),
S3,5={0,−1/2,−2/3,4/3} log13(3).

Together with §5.29 (vanishing at5,11,19,53 and all good primes) and the dyadic proof above, these are the only away contributions.
The QC target sets are Ωδ=−S3,δ.
Since log13(3)=2093 modulo2197:
Ω1={0,1482,1534},
Ω5={0,2145,663,871} modulo2197.

The exact sets are the rational multiples of the logarithm just displayed; the residues are only their precision-three reductions. This calculation concerns curve-level QC heights and does not erase nonzero local contributions previously used for genus-two generator pairings.

Independent review checked the dyadic analytic argument, 3-adic geometry and exhaustive accessibility, network calibration, and correspondence sign. The remaining tasks are the moving normalized 13-adic functions, certified root isolation and final rational-point sieve.
