# A certified denominator obstruction for type 18

Let `P(X)=(X−1)(X−9)(X−25)(X−49)(X−81)`. A prime q>5 is called rigid here if zero is the only c in F_q for which all five roots of `P(X)−c` are square residues in F_q, including zero and counting multiplicity. The complete classifier finds 521 such primes, the largest being 26431. Their product gives the following bound for any type-18 configuration with distinct, disjoint chain sets:

    D > 2.07 × 10^153,

where D is the common reduced denominator of the five doubled B centers.

## From rigid primes to a lower bound

By [the individual-square argument](PROOF.md), all five doubled B centers of any rational tuple have a common reduced denominator D, and every prime dividing D is `1 mod 10`. Write the centers as `n_i/D` and `N=∏n_i`. Newton's identities give

    ∏(X−n_i²) = D^10 P(X/D²) − Δ,
    Δ = N² − 945² D^10.

If q is rigid and prime to D, reduce after dividing out the powers of D. The resulting root polynomial is `P(X)−Δ/D^10`, with all five roots squares in F_q. Hence q divides Δ. If q divides D then q divides `DΔ` directly. Thus for any finite certificate of rigid primes, writing M for their product and M0 for the product restricted to `q ≠ 1 mod 10`, one has

    M0 | Δ,      M | DΔ.

The nonnegative numbers `n_i²/D²` have sum 165, so the arithmetic-geometric mean inequality gives `N²≤33^5 D^10`. Since `0<945²<33^5`, this implies `|Δ|<33^5 D^10`. When Δ=0 the roots are exactly those of the seed and the two chain sets coincide. Every noncolliding tuple must therefore satisfy

    M0 < 33^5 D^10,      M < 33^5 D^11.

See [Coppersmith–Mossinghoff–Scheinerman–VanderKam, §§3–4](https://arxiv.org/html/2304.11254v1) for the established method of using local patterns to constrain the constant of a PTE solution. The certificate here specializes that method to the rational type-18 denominator.

## An absolute cutoff for rigidity

For every prime q>28,955,161 there is a nonzero c for which `P(X)−c` has five distinct nonzero roots, each a square in F_q. Equivalently `P(Y²)−c` has ten distinct roots in F_q.

To prove this, put `R(c)=Res_X(P(X)−c,P'(X))`. The exact identities below ensure that the four critical points of P are simple, their critical values are distinct, and none of those values is P(0), for every prime outside the displayed prime supports. The identities also give `P'(0)≠0` there.

| Quantity | Factorization |
|---|---|
| `disc(P')` | `2^42 · 3^5 · 5 · 11^3 · 47 · 59` |
| `disc(R)` | `2^202 · 3^27 · 5^11 · 11^15 · 47^3 · 59^3` |
| `R(P(0))` | `3^10 · 5^6 · 7^2 · 11^5 · 59^2 · 181^2 · 1231 · 2417` |
| `P'(0)` | `3^2 · 11 · 59 · 181` |

The degree-five cover `X ↦ P(X)` is connected over the algebraic closure, so its geometric monodromy is transitive on five roots and hence primitive. At each simple critical value inertia acts as a transposition. A primitive group containing a transposition is S5. Indeed the graph of all conjugate transpositions must be connected because its components are invariant blocks, and the edge transpositions of a connected graph generate the symmetric group.

The geometric group of `Y ↦ P(Y²)` projects onto this S5 and preserves the pairs `{+Y,−Y}` above each root X. At `c=P(0)` only the ramification at Y=0 occurs. It yields one transposition inside a single pair. Conjugates give all five independent pair flips, so the geometric group is the entire signed permutation group `(C2)^5 ⋊ S5` of order 3840. The arithmetic group embeds in the same signed permutation group and contains the geometric group. They are equal, so the splitting field over F_q(c) is regular and its smooth projective curve is geometrically connected.

This Galois cover of the c-line is tame because q does not divide 3840. Over each of the four critical values of P the two ramified points of `Y ↦ P(Y²)` have index two, so inertia on the ten leaves has order two. Inertia at P(0) has order two and at infinity it has order ten. The Riemann–Hurwitz formula gives

    2g−2 = 3840·(−2 + 5·(1−1/2) + (1−1/10)) = 5376,
    g = 2689.

At most `5·1920 + 384 + 3840 = 13824` geometric points lie over the branch values, infinity, or c=0. For `q≥5381²=28955161`, the Hasse–Weil lower bound for the number of F_q-points is at least `5381²+1−5378·5381=16144`. Thus there is an F_q-point over a nonzero unbranched finite c. At an unbranched place the existence of a rational point of the Galois cover makes the decomposition group trivial. Specialization there gives all ten distinct roots in F_q. This proves the cutoff.

## Reproduction

The exact coefficient and branch checks are part of `scripts/check_all_shift_chain_types.py`. The C++ classifier counts all distinct square residues once and adds any missing multiplicities at critical roots. On large fields it first tries to exhibit a nonzero split fiber. After choosing a square residue and its corresponding shift, the remaining monic quartic can be certified by divisibility into `X^((q+1)/2)−X`. The sampling order only affects speed. The exhaustive histogram handles every prime where the preliminary step gives no witness.

    g++ -O3 -march=native -fopenmp -std=c++20 scripts/scan_type18_rigid_primes.cpp -o build/scan_type18
    OMP_NUM_THREADS=8 build/scan_type18 28955161 > build/type18_rigid_primes.csv
    cmp build/type18_rigid_primes.csv scripts/type18_rigid_primes.csv
    python3 scripts/verify_type18_rigid_primes.py scripts/type18_rigid_primes.csv --denominator-bound "$(python3 -c 'print(207 * 10**151)')"

The Python verifier uses derivatives to recompute the multiplicity of every square residue root for every retained prime. It proves the denominator bound from any verified subset of rigid primes. Exhaustiveness of the 521-prime list additionally uses the C++ classification of all 1,797,076 primes from 7 through the cutoff and the global cutoff proof above. The full scan took about 81 seconds using eight threads on an i7-14650HX. The prime calculation gives the lower bound on D. The later [rank-one exclusion in Theorem 49](PROOF.md#576-complete-exclusion-of-type-18-and-rigidity-of-a-ten-term-progression) rules out a noncolliding rational tuple for every D.
