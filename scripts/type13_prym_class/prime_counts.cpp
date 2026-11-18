// Count prime-ideal powers in rational bins for sections 5.55 and 5.60.
// Usage: prime_counts FIRST LAST [type17]. The interval is in rational primes.
#include <pari/pari.h>

#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <tuple>

constexpr unsigned long cutoff = 1UL << 31;
constexpr char field[] =
    "x^12-4*x^11-68*x^10-18698*x^9+336998*x^8-500324*x^7"
    "+7039932*x^6-194015520*x^5+1041300636*x^4-3497537936*x^3"
    "+20028904184*x^2-109223576168*x+269574001992";

using Bin = std::tuple<unsigned, unsigned, unsigned>;
using Histogram = std::map<Bin, std::uint64_t>;

void record_powers(Histogram& counts, std::uint64_t norm, unsigned multiplicity, unsigned long limit) {
    for (std::uint64_t n = norm, m = 1; n <= limit; ++m) {
        unsigned exponent = std::bit_width(n) - 1;
        unsigned bin = ((n << 8) >> exponent) - 256;
        counts[{exponent, bin, static_cast<unsigned>(m)}] += multiplicity;
        if (n > limit / norm) break;
        n *= norm;
    }
}

int main(int argc, char** argv) {
    if (argc != 3 && (argc != 4 || std::string(argv[3]) != "type17")) return 2;
    const bool type17 = argc == 4;
    const unsigned long limit = type17 ? 1UL << 20 : cutoff;
    unsigned long lo = std::stoul(argv[1]), hi = std::stoul(argv[2]);
    if (lo < 2 || lo > hi || hi > limit) return 2;
    pari_init(160000000, 1000000);
    GEN poly = gp_read_str(type17
        ? "x^12-472*x^10+75282*x^8-4450612*x^6+64570537*x^4-277288596*x^2+234027000"
        : field);
    GEN nf = nfinit(poly, DEFAULTPREC);
    GEN index = nf_get_index(nf);
    GEN bound = utoipos(limit);
    Histogram counts;
    forprime_t primes;
    if (u_forprime_init(&primes, lo, hi)) {
        pari_sp clean = avma;
        for (ulong p; (p = u_forprime_next(&primes));) {
            if (p <= limit / p || umodiu(index, p) == 0) {
                GEN ideals = idealprimedec_limit_norm(nf, utoipos(p), bound);
                for (long i = 1; i < lg(ideals); ++i)
                    record_powers(counts, itou(pr_norm(gel(ideals, i))), 1, limit);
            } else {
                // Dedekind factorization applies. Only degree one can contribute.
                long roots = Flx_nbroots(ZX_to_Flx(poly, p), p);
                if (roots) record_powers(counts, p, roots, limit);
            }
            set_avma(clean);
        }
    }
    for (const auto& [key, count] : counts) {
        auto [exponent, bin, m] = key;
        std::cout << exponent << ' ' << bin << ' ' << m << ' ' << count << '\n';
    }
    pari_close();
}
