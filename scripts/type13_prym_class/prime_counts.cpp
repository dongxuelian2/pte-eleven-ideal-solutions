// Count prime-ideal powers in rational bins for the class-group certificates.
// Usage: prime_counts FIRST LAST [type8|type17|type18]. The interval is in rational primes.
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
    if (argc < 3 || argc > 4) return 2;
    const std::string mode = argc == 4 ? argv[3] : "type13";
    if (mode != "type8" && mode != "type13" && mode != "type17" && mode != "type18") return 2;
    const unsigned long limit = mode == "type8" ? 1UL << 22 :
                                mode == "type17" ? 1UL << 20 : cutoff;
    unsigned long lo = std::stoul(argv[1]), hi = std::stoul(argv[2]);
    if (lo < 2 || lo > hi || hi > limit) return 2;
    pari_init(160000000, 1000000);
    GEN poly = gp_read_str(mode == "type8"
        ? "x^12-6*x^11-66*x^10-1415*x^9-3987*x^8+128772*x^7+1543187*x^6"
          "-30854580*x^5+671675061*x^4-3249300991*x^3+46170279879*x^2"
          "-461672478855*x+1325918174175"
        : mode == "type17"
        ? "x^12-472*x^10+75282*x^8-4450612*x^6+64570537*x^4-277288596*x^2+234027000"
        : mode == "type18"
        ? "x^12-3*x^11-5625*x^10+44481*x^9+9843099*x^8-99843252*x^7"
          "-5998855201*x^6+67870752036*x^5+854290494639*x^4-9750057898095*x^3"
          "-2947764849192*x^2+213431600705238*x-420745260940011"
        : field);
    GEN nf = nfinit(poly, DEFAULTPREC);
    // In the ordered-pair field each pair of different fixed quartic roots
    // gives one degree-one prime away from the quartic discriminant.
    const bool pairs = mode == "type18";
    GEN root_poly = pairs ? gp_read_str("x^4-x^3-991*x^2+13345*x-19455") : poly;
    GEN exceptional = pairs ? poldisc0(root_poly, -1) : nf_get_index(nf);
    GEN bound = utoipos(limit);
    Histogram counts;
    forprime_t primes;
    if (u_forprime_init(&primes, lo, hi)) {
        pari_sp clean = avma;
        for (ulong p; (p = u_forprime_next(&primes));) {
            if (p <= limit / p || umodiu(exceptional, p) == 0) {
                GEN ideals = idealprimedec_limit_norm(nf, utoipos(p), bound);
                for (long i = 1; i < lg(ideals); ++i)
                    record_powers(counts, itou(pr_norm(gel(ideals, i))), 1, limit);
            } else {
                // Dedekind factorization applies. Only degree one can contribute.
                long roots = Flx_nbroots(ZX_to_Flx(root_poly, p), p);
                if (pairs) roots *= roots - 1;
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
