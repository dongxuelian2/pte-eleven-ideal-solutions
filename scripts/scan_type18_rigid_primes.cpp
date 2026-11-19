// g++ -O3 -march=native -fopenmp -std=c++20 scan_type18_rigid_primes.cpp -o scan_type18
#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {
using U32 = std::uint32_t;
using U64 = std::uint64_t;
constexpr std::array<int, 6> seed = {-893025, 1057221, -172810, 8778, -165, 1};

// All inputs to reduce are unsigned 64-bit integers. With a floor reciprocal
// the estimated quotient can be low by only one.
struct Modulus {
    U32 prime;
    U64 reciprocal;

    explicit Modulus(U32 q) : prime(q), reciprocal(std::numeric_limits<U64>::max() / q) {}

    U32 reduce(U64 value) const {
        const U64 quotient = (static_cast<__uint128_t>(value) * reciprocal) >> 64;
        const U64 rest = value - quotient * prime;
        return rest >= prime ? rest - prime : rest;
    }

    U32 coefficient(int value) const {
        const int q = static_cast<int>(prime);
        const int rest = value % q;
        return rest < 0 ? rest + q : rest;
    }
};

template <std::size_t N>
U32 evaluate(const std::array<U32, N>& coefficients, U32 x, const Modulus& mod) {
    U32 value = coefficients.back();
    for (std::size_t i = N - 1; i > 0; --i) {
        value = mod.reduce(static_cast<U64>(value) * x + coefficients[i - 1]);
    }
    return value;
}

// Since x is a root, synthetic division computes the quotient by X-x.
// Repeating until the quotient is nonzero at x gives the root multiplicity.
unsigned multiplicity(std::array<U32, 6> coefficients, U32 x, const Modulus& mod) {
    for (std::size_t degree = 5; degree > 0; --degree) {
        for (std::size_t i = degree - 1; i > 0; --i) {
            coefficients[i] = mod.reduce(static_cast<U64>(coefficients[i + 1]) * x
                                         + coefficients[i]);
        }
        if (mod.reduce(static_cast<U64>(coefficients[1]) * x + coefficients[0]) != 0) {
            return 5 - degree;
        }
        // The first pass checks the remainder of the original polynomial.
        // The next pass uses the coefficients of its quotient.
        for (std::size_t i = 0; i < degree; ++i) {
            coefficients[i] = coefficients[i + 1];
        }
    }
    return 5;
}

// Find a single split fiber before paying for a complete histogram at a large
// prime. The sampling order merely affects speed. Each returned witness is
// proved by polynomial divisibility, with a full histogram as the fallback.
bool all_square_roots(const std::array<U32, 6>& polynomial, U32 known_root,
                      const Modulus& mod) {
    std::array<U32, 4> neg_quotient;
    U32 coefficient = 1;
    for (unsigned i = 4; i > 0; --i) {
        coefficient = mod.reduce(static_cast<U64>(coefficient) * known_root + polynomial[i]);
        neg_quotient[i - 1] = coefficient == 0 ? 0 : mod.prime - coefficient;
    }
    // The monic quartic after dividing P(X)-P(known_root) by X-known_root
    // divides X^((q+1)/2)-X precisely when its roots are distinct squares.
    const U32 exponent = (mod.prime + 1) / 2;
    std::array<U32, 4> x_power = {0, 1, 0, 0};
    for (U32 bit = std::bit_floor(exponent) >> 1; bit != 0; bit >>= 1) {
        std::array<U64, 7> square = {};
        for (unsigned i = 0; i < 4; ++i) {
            square[2 * i] += static_cast<U64>(x_power[i]) * x_power[i];
            for (unsigned j = i + 1; j < 4; ++j) {
                square[i + j] += 2ULL * x_power[i] * x_power[j];
            }
        }
        for (unsigned degree = 6; degree >= 4; --degree) {
            const U64 lead = mod.reduce(square[degree]);
            for (unsigned j = 0; j < 4; ++j) {
                square[degree - 4 + j] += lead * neg_quotient[j];
            }
        }
        if (exponent & bit) {
            const U64 lead = mod.reduce(square[3]);
            for (unsigned i = 3; i > 0; --i) {
                x_power[i] = mod.reduce(square[i - 1] + lead * neg_quotient[i]);
            }
            x_power[0] = mod.reduce(lead * neg_quotient[0]);
        } else {
            for (unsigned i = 0; i < 4; ++i) {
                x_power[i] = mod.reduce(square[i]);
            }
        }
    }
    return x_power == std::array<U32, 4>{0, 1, 0, 0};
}

bool rigid_prime(U32 q, std::vector<std::uint8_t>& counts) {
    const Modulus mod(q);
    std::array<U32, 6> polynomial;
    std::transform(seed.begin(), seed.end(), polynomial.begin(),
                   [&](int c) { return mod.coefficient(c); });
    if (q >= 250'000) {
        U32 sample = q;
        for (unsigned trial = 0; trial < 2048; ++trial) {
            sample = 1664525U * sample + 1013904223U;
            const U32 x = mod.reduce(sample);
            const U32 square = mod.reduce(static_cast<U64>(x) * x);
            if (evaluate(polynomial, square, mod) != 0 && all_square_roots(polynomial, square, mod)) {
                return false;
            }
        }
    }
    counts.assign(q, 0);

    // z=0,...,(q-1)/2 represents every square exactly once. Maintain z²
    // modulo q by addition. The c=0 fiber always splits, so it can be skipped.
    U32 square = 0;
    for (U32 z = 0; z <= (q - 1) / 2; ++z) {
        const U32 shift = evaluate(polynomial, square, mod);
        if (shift != 0 && ++counts[shift] == 5) {
            return false;
        }
        square += 2 * z + 1;
        if (square >= q) {
            square -= q;
        }
    }

    // Distinct roots never overcount. Add missing multiplicities only at the
    // at most four critical roots, and only when the first pass was inconclusive.
    std::array<U32, 5> derivative;
    for (unsigned i = 1; i <= 5; ++i) {
        derivative[i - 1] = mod.reduce(static_cast<U64>(i) * polynomial[i]);
    }
    square = 0;
    for (U32 z = 0; z <= (q - 1) / 2; ++z) {
        if (evaluate(derivative, square, mod) == 0) {
            const U32 shift = evaluate(polynomial, square, mod);
            if (shift != 0) {
                auto shifted = polynomial;
                shifted[0] = mod.reduce(static_cast<U64>(shifted[0]) + q - shift);
                counts[shift] += multiplicity(shifted, square, mod) - 1;
                if (counts[shift] == 5) {
                    return false;
                }
            }
        }
        square += 2 * z + 1;
        if (square >= q) {
            square -= q;
        }
    }
    return true;
}

std::vector<U32> primes_in_range(U32 first, U32 limit) {
    std::vector<bool> composite(limit + 1);
    for (U32 q = 2; q <= limit / q; ++q) {
        if (!composite[q]) {
            for (U32 k = q * q; k <= limit; k += q) {
                composite[k] = true;
            }
        }
    }
    std::vector<U32> primes;
    for (U32 q = first | 1U; q <= limit; q += 2) {
        if (!composite[q]) {
            primes.push_back(q);
        }
    }
    return primes;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: scan_type18 PRIME_LIMIT [FIRST_PRIME] > rigid_primes.csv\n";
        return 2;
    }
    std::size_t parsed = 0;
    unsigned long limit, first = 7;
    try {
        limit = std::stoul(argv[1], &parsed);
        if (argv[1][parsed] != '\0' || limit < 7 || limit > 29'000'000) {
            throw std::out_of_range("prime bound");
        }
        if (argc == 3) {
            first = std::stoul(argv[2], &parsed);
            if (argv[2][parsed] != '\0' || first < 7 || first > limit) {
                throw std::out_of_range("first prime");
            }
        }
    } catch (const std::exception&) {
        std::cerr << "Require 7 <= FIRST_PRIME <= PRIME_LIMIT <= 29000000\n";
        return 2;
    }

    const auto start = std::chrono::steady_clock::now();
    const auto primes = primes_in_range(first, limit);
    std::vector<std::uint8_t> rigid(primes.size());
    #pragma omp parallel
    {
        std::vector<std::uint8_t> counts;
        #pragma omp for schedule(dynamic, 64)
        for (std::size_t i = 0; i < primes.size(); ++i) {
            rigid[i] = rigid_prime(primes[i], counts);
        }
    }

    unsigned retained = 0;
    std::cout << "q\n";
    for (std::size_t i = 0; i < primes.size(); ++i) {
        if (rigid[i]) {
            std::cout << primes[i] << '\n';
            ++retained;
        }
    }
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    std::cerr << "Checked " << primes.size() << " primes in [" << first << ',' << limit << ']'
              << ", retained " << retained << " in " << elapsed.count()
              << " s with up to " << omp_get_max_threads() << " threads\n";
    return std::cout.good() ? 0 : 1;
}
