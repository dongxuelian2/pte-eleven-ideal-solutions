"""Run PARI 2.19.1's factor-base coverage check on disjoint rational-prime intervals."""
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess


def main():
    args = ArgumentParser(description=__doc__)
    args.add_argument('source', type=Path, help='PARI src/basemath/buch2.c')
    args.add_argument('include', type=Path, help='directory containing pari.h')
    args.add_argument('library', type=Path, help='directory containing libpari.so')
    args.add_argument('build', type=Path, help='directory for generated source and binary')
    args.add_argument('--jobs', type=int, default=1)
    args.add_argument('--type18', action='store_true', help='check the Type 18 vertex field')
    opts = args.parse_args()
    bound = 85_000_000 if opts.type18 else 48_000_001
    opts.build.mkdir(parents=True, exist_ok=True)
    text = opts.source.read_text()
    header = 'bnftestprimes(GEN bnf, GEN BOUND)\n{'
    start = text.index(header)
    end = text.index('set_avma(av0);\n}', start)
    routine = text[start:end]
    iterator = 'forprime_init(&S, gen_2, BOUND);'
    assert routine.count(iterator) == 1
    routine = routine.replace(header, 'range_primes(GEN bnf, GEN lo, GEN hi, GEN BOUND)\n{')
    routine = routine.replace(iterator, 'forprime_init(&S, lo, hi);')
    source = opts.build / 'range_pari.c'
    source.write_text(text[:start] + routine + text[end:])
    obj = opts.build / 'range_pari.o'
    binary = (opts.build / 'cover_range').resolve()
    flags = ['cc', '-O3', '-ffunction-sections', '-fdata-sections']
    subprocess.run([*flags, '-I', opts.include, '-c', source, '-o', obj], check=True)
    worker = Path(__file__).with_name('cover_range.c')
    if opts.type18:
        flags.append('-DTYPE18')
    subprocess.run([*flags, '-I', opts.include.parent, worker, obj, '-Wl,--gc-sections',
                    '-L', opts.library, f'-Wl,-rpath,{opts.library.resolve()}',
                    '-lpari', '-lm', '-o', binary], check=True)

    def check(lo):
        hi = min(bound, lo+999_999)
        result = subprocess.run([binary, str(lo), str(hi)], text=True,
                                capture_output=True, check=True)
        assert result.stdout.strip() == f'PASS {lo} {hi}'
        return hi

    with ThreadPoolExecutor(max_workers=opts.jobs) as pool:
        for hi in pool.map(check, range(2, bound+1, 1_000_000)):
            print(f'Certified through {hi}', flush=True)


if __name__ == '__main__':
    main()
