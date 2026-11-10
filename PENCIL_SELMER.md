# A three-dimensional 2-Selmer group for a splitting pencil

## Status

The staged calculation leaves exactly eight global squareclasses. Its present global completeness uses the Magma class-group calculation under GRH. With the verified local conditions described below, this computes a three-dimensional 2-Selmer group under that assumption. Consequently the Jacobian rank satisfies 1 <= rank <= 3 under GRH. The unconditional lower bound rank >= 1 remains unchanged.

This is not a rank-one certificate, a rational-point classification, an exclusion of this construction, or a solution of the original eleven-integer problem. No conclusion about finiteness of Sha is assumed.

## Construction and theorem

Put P(x)=x⁵−110x³+2629x−2520=∏(x−a), where
a ranges over {−8,−7,1,5,9}, and Q(x)=x⁶−110x⁴+2629x²+2520x.
The even part of P(x)(Q(x)+t) is the constant −2520t. If Q+t
has six distinct rational roots, none among those five roots of P, and
t≠0, the eleven roots A of P(Q+t) satisfy
f_A(x)+f_A(−x)=−5040t≠0. Thus A and −A are disjoint; clearing
denominators constructs a distinct-entry ideal size-eleven integer pair.
The five collision parameters are {40320,35280,−5040,−25200,−45360}.

The discriminant of Q+t is 64f(t), where f is the quintic below.
Complete rational splitting therefore requires a rational point on C:y²=f(t).
This necessary discriminant condition is not sufficient for splitting.

**Theorem (conditional on GRH for the class-group computation).**
For this curve C, dim_F2 Sel²(J/Q)=3 and 1≤rank J(Q)≤3.
Unconditionally J(Q)_tors=0 and rank J(Q)≥1.

This uses the standard odd-degree hyperelliptic 2-descent identification
of the Selmer group with norm-square classes in K*/K*² satisfying all
local Kummer-image conditions. No assertion of novelty for the descent
method is made. The computational statement concerns this specific pencil.

## Curve and field

The original curve is

y² = -729t⁵ + 6403320t⁴ + 930770808176t³ + 3787455316296960t² - 27536385315056910336t + 15173754195676653158400.

Use the monic model F(X) defined in the local scripts and

K = Q(a), a⁵ - 660a³ + 70983a - 102060 = 0,
θ = (-165a⁴ + 70983a² - 255150a - 750870)/2.

Each script verifies F(θ)=0. The transformation is X=-729t/4-375435, Y=531441y/32; thus F(X)=(531441/32)²f(t). There is a rational point at infinity. The irreducible quintic gives J(Q)[2]=0. At the good primes 11 and 13, exact counts give (N1,N2,#J)=(3,123,55) and (12,206,162), respectively, with #J=(N1²+N2)/2−p. Prime-to-p torsion injects into reduction; these coprime orders also exclude the 11- and 13-primary parts using the other prime. Hence J(Q)_tors=0. The point (0,123181793280) gives a nonzero class [P−∞]: a principal divisor P−∞ would give a degree-one map from a genus-two curve to P¹. It therefore has infinite order.

## Global squareclass space

Magma's GRH class-group calculation gives Cl(K)=0. This is the conditional step. Select nine prime ideals:

- the three unramified degree-one primes above 3;
- the degree-four prime above each of 5 and 7;
- the two degree-one primes above 773;
- the two degree-one primes above 3087467 containing θ-306176.

The local-image calculations justify this support. In particular the full image at 2 has even valuations, and the full images at 61 and 35498020282877 also have even valuations. At 3 we initially keep a safe larger space than needed, then impose the full local image later.

`SUnitGroup(S:Raw:=true)` returns a rank-thirteen group with torsion {±1}, hence fourteen generators modulo squares. Reducing its raw exponent vectors modulo two BEFORE expansion preserves every squareclass, while avoiding huge full-unit expansions. The complete export, including exact field coefficients and valuation vectors, is `scripts/pencil_selmer/selected.txt`.

Completeness: every Selmer class has parity support in the selected S. After removing these odd valuations with S-unit representatives, its fractional ideal is a square. Trivial class group makes that fractional ideal principal, so another field square removes it, leaving a global unit, already included in the S-unit group.

The three valuation equations are v1+v2+v3=0, v6+v7=0, v8+v9=0. In increasing-a real order the local image is span{10100,01010}, giving three additional equations. Rational interval evaluation certifies all signs. The six equations have rank six and leave dimension eight. Allowed valuations make every rational norm exponent even; the real condition makes its sign positive, so the norm is a rational square.

## Full local images

For odd p, each field component squareclass is recorded by valuation parity and the square/nonsquare character of its normalized residue unit. The scripts use completions, with explicit unit normalization. At 2 the algorithm compares squareclasses directly, normalizing only by EVEN powers of uniformizers.

For each local factor fi of F of degree d, the rational 2-torsion Kummer representative is (-1)^d fi(θ) outside its component and (-1)^(d-1)(F/fi)(θ) inside. Products of the other factors implement F/fi. The inside formula differs from the literal product only by fi'(θ)². The scripts identify the matching component to high precision and certify finite valuations on the others.

Supplementary point/divisor witnesses, all interpreted in the monic X coordinate:

- p=2: X=5/4, 2, 885;
- p=5: u(X)=X²+85X+23;
- p=7: u(X)=X²+22X+78271;
- p=773: X=205278;
- p=3087467: X=48537585605194255596273016 and u(X)=X²+2190851906626X+19753949860961491482589978;
- p=61: X=2657;
- p=35498020282877: u(X)=X²+401107556428946897390988312X+895917242855804582442647861843282088065835134709300681.

The degree-one Kummer classes are X-θ, not θ-X. Quadratic classes are u(θ). To validate a quadratic witness one checks F modulo u is a square in EACH factor field, not merely that its norm is square.

At odd p the generated local-image dimension equals dim J(Qp)[2], namely the number of local factors minus one. At 2 it equals 2+dim J(Q2)[2]=4. Thus the displayed witnesses span full images; a point search alone is not being used as a completeness assertion.

The intersections are:

- global valuation and real constraints: 256 classes, dimension 8;
- p5: 128, dimension 7;
- p7: 64, dimension 6;
- p3: 16, dimension 4;
- p773 and p3087467: still 16;
- p2: 8, dimension 3;
- p61 and p35498020282877: all eight survive.

The final masks in the fourteen-generator ordering are
0,4,3169,3173,12298,12302,15467,15471,
with basis 4,3169,12298. At all other primes, the classes are unramified and norm-square, which is exactly the local image at odd primes of good reduction. Real conditions have already been imposed. This establishes the stated conditional Selmer dimension.

## Precision and independent verification

For odd p, a normalized unit's nonzero residue determines its squareclass. For the dyadic completions e<=2, a unit modulo π^N with N>2e has determined squareclass. The script requires unit precision>16, safely beyond this threshold. Each discarded power of π is even. Higher-precision audit scripts strengthen component-uniqueness assertions and rerun at factorization precision400 and completion precision450.

Independent Python checks reproduce all real signs with rational isolating intervals, all p5/p7 squareclasses with exact modular arithmetic, and the supplemental witnesses. Independent linear algebra reproduces each finite intersection. The dyadic implementation and precision argument were separately audited; it is not claimed that the ramified dyadic square tests have a second independent Python implementation.

The successful mathematical outputs from Magma 2.30-1 are retained as plain text beside the code. Service envelopes, failed calculations and bounded point searches are omitted. The odd and dyadic files retain the higher-precision calculation.

## Next mathematical step

Dimension three does not decide whether the remaining classes come from rational Jacobian points. A nonzero Cassels–Tate pairing on the two remaining independent directions could lower the rank upper bound to one. Alternatively, finding two further independent rational divisor classes would establish rank three. Neither has been done here. A rank-one result would still require subgroup/saturation checks and a complete rational-point computation before classifying the sextic splitting parameters.

## Reproduction

Run `python scripts/pencil_selmer/verify.py` from the repository root
(Python 3 with SymPy). This verifies construction identities, good-reduction
counts, real signs, exact p=5,7 local characters, all supplemental point
and quadratic-divisor witnesses, every intersection, and the exact identity
linking the known point to mask 4. Generated intermediate files go to
`build/pencil_selmer/`, not the source tree.

The `.m` files are independent Magma inputs. Run `selected_sunits.m` for
the conditional class-group and complete S-unit calculation; the four
local files use the fixed exported fourteen representatives. Their recorded
outputs are mathematical certificate inputs to Python, not a replacement
for rerunning Magma if independent verification of its number-field and
completion arithmetic is required. Recomputing S-units can change their
basis and ordering; in that case re-export the basis into all local files.
The original and higher-precision outputs agree, including the eight dyadic
survivors. The package does not claim Python independently implements
Magma's class group or ramified dyadic square tests.

The known point has squareclass mask 4. The five coefficients in
`known_root.txt` define r∈K satisfying r²=(−375435−θ)G3 exactly;
`verify.py` checks this polynomial identity over Q, independently of Magma.

For the standard descent conventions see M. Stoll, “Implementing 2-descent
for Jacobians of hyperelliptic curves”, Acta Arithmetica 98 (2001), 245–277,
and the [Magma Jacobian handbook](https://docs.magma-maths.org/ArithmeticGeometry/HyperellipticCurves/jacobians-number-fields.html).
