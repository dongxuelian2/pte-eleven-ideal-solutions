# Ten-term progressions and single shifts in the ideal Prouhet–Tarry–Escott problem

Daoyu Dong and Qichao Wang study rational and integer sets with matching power sums. A rational arithmetic progression of length ten is determined among ten-element rational multisets by its first nine power sums.

For a symmetric ideal size-eleven target, every distinct-entry precursor under a single shift has at least twelve entries per side. The proof excludes all eighteen chain types of a size-ten precursor. For any ideal size-eleven pair of distinct integer sets, pairing their entries in increasing order gives total distance at least 22, or 24 for a symmetric pair. The unrestricted integer existence question remains open.

## Read the paper

Start with the [main theorems](papers/pte-rigidity/introduction.tex). The [marked-vertex proof](papers/pte-rigidity/progression.tex) establishes progression rigidity. For the shift obstruction, the [partition table](papers/pte-rigidity/shifts.tex) leads to the [individual models](papers/pte-rigidity/cases.tex) and the [Type 10 quadratic Chabauty calculation](papers/pte-rigidity/type10.tex). The [distance proof](papers/pte-rigidity/distance.tex) applies to ideal integer pairs of any size.

Build the [complete manuscript](papers/pte-rigidity/main.tex) with LuaLaTeX, Latin Modern and DejaVu Serif. From the repository root:

~~~sh
cd papers/pte-rigidity
latexmk -lualatex -outdir=../../build/pte-rigidity main.tex
~~~

The resulting PDF is at <code>build/pte-rigidity/main.pdf</code> relative to the repository root.

## Verify a result

The [computational appendix](papers/pte-rigidity/verification.tex) gives the dependencies, commands and expected results. To check the partition enumeration and polynomial identities, run the following from the repository root with Python 3 and SymPy:

~~~sh
python3 scripts/check_all_shift_chain_types.py
~~~

The [marked-vertex certificate instructions](scripts/type18_descent/README.md) cover the class-group bound, rank and dyadic logarithm. The [Type 10 proof](proofs/type-10-quadratic-chabauty.md) explains how the local roots enter the final elliptic sieve.

## Further arithmetic

The [proof dossier](PROOF.md) contains derivations, exceptional fibers and results on common-center pairs and multisets. The [supporting proofs](proofs/) give the local height computations, inverse-shift argument and rigid-prime bounds. The independent [splitting-pencil descent](proofs/splitting-pencil.md) bounds its Jacobian rank between one and three. The [Kummer fibration calculation](kummer_obstruction/README.md) records the obstruction for its seven models.
