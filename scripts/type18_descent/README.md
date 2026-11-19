# The vertex marked by 25

Theorem 49 uses the vertex assignment in Proposition 47, a class-group
bound, a rank-one computation and an exact reduction of a global point.
From the repository root, use Magma V2.30-1, PARI/GP 2.19.1, Python 3
and a C++20 compiler. With PARI on the compiler's standard paths, run

```sh
mkdir -p build/type18-class
magma scripts/check_type18_product_classes.m
c++ -O3 -std=c++20 scripts/type13_prym_class/prime_counts.cpp -lpari \
    -o build/type18-class/prime_counts
python3 scripts/type13_prym_class/check_bound.py build/type18-class/prime_counts --type18 --jobs 12
gp -fq scripts/type18_descent/check_factor_base.gp
magma scripts/type18_descent/rank.m
magma scripts/type18_descent/pairing.m
gp -fq scripts/type18_descent/logarithm.gp
```

The explicit-formula margin exceeds 0.284. The sequential GP factor-base
check can also be run on disjoint rational-prime intervals with the matching
PARI source and library. For example,

```sh
python3 scripts/type13_prym_class/cover.py /path/to/pari-2.19.1/src/basemath/buch2.c \
    /usr/include/pari /usr/lib build/type18-class/ranges --type18 --jobs 12
```

The last interval ends at 85,000,000. After this check the bound 1000 in
the Magma script is unconditional. The pairing script checks both explicit quartics for everywhere local
solubility, verifies their Jacobians using the binary-quartic invariants,
and evaluates their pairing.
The class group has order two. The Selmer group has dimension three and
the pairing on the two-dimensional quotient has value one. The final GP
script computes the coordinates of 24 times the seed over the exact number
field. Its formal parameter has valuation three and leading residue
`t^3 + t` in the residue field of order sixteen.
