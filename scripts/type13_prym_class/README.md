# Reproducing the Type 13 class-group bound

Proposition 27 uses PARI/GP 2.19.1, a C++20 compiler and Python 3.
Run the following commands from the repository root with PARI's headers
and shared library on the standard compiler paths.

```sh
mkdir -p build/type13-class
c++ -std=c++20 -O3 scripts/type13_prym_class/prime_counts.cpp -lpari \
    -o build/type13-class/prime_counts
python3 scripts/type13_prym_class/check_bound.py build/type13-class/prime_counts --jobs 12
gp -fq scripts/type13_prym_class/check_factor_base.gp
```

The first check enumerates prime-ideal powers and proves the explicit-formula
inequality using rational intervals. The GP check verifies ideal relations
over the resulting generating range. For a parallel run of the same relation
test, supply the PARI source file matching the installed library, its header
directory, the library directory and an output directory, for example

```sh
python3 scripts/type13_prym_class/cover.py /path/to/pari-2.19.1/src/basemath/buch2.c \
    /usr/include/pari /usr/lib build/type13-class --jobs 12
```

The worker changes the start and end of the outer prime iterator in
`bnftestprimes` and preserves its factorization and ideal-reduction code.
It checks the norm of every factor-base ideal in each process.
After the two PARI checks, run the two Magma certificates.

```sh
magma scripts/check_type13_beta_prym_rank.m
magma scripts/check_type13_beta_prym.m
```

The first gives the class groups and the Selmer rank bound. The second
checks the independent divisor classes, the 23-adic determinants and
the simultaneous signed reductions at 1181.
