// An interval worker for the two-line PARI source adaptation in cover.py.
#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>

void range_primes(GEN bnf, GEN lo, GEN hi, GEN bound);

int main(int argc, char **argv) {
    if (argc != 3) return 2;
    const unsigned long lo = strtoul(argv[1], NULL, 10);
    const unsigned long hi = strtoul(argv[2], NULL, 10);
    if (lo < 2 || lo > hi || hi > 48000001) return 2;
    pari_init(160000000, 1000000);
    GEN pol = gp_read_str("x^12-4*x^11-68*x^10-18698*x^9+336998*x^8-500324*x^7"
        "+7039932*x^6-194015520*x^5+1041300636*x^4-3497537936*x^3"
        "+20028904184*x^2-109223576168*x+269574001992");
    GEN bnf = bnfinit0(pol, 1, NULL, DEFAULTPREC);
    GEN base = gel(bnf, 5);
    for (long i = 1; i < lg(base); ++i)
        if (cmpiu(pr_norm(gel(base, i)), 1000) > 0) return 3;
    range_primes(bnf, utoipos(lo), utoipos(hi), utoipos(48000001));
    printf("PASS %lu %lu\n", lo, hi);
    pari_close();
}
