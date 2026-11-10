# Twisted Kummer fibration obstruction

The proof in KUMMER_FIBRATION_OBSTRUCTION.md excludes genus-one fibrations over Q, finite odd-degree extensions, and quadratic extensions for all seven nonzero Selmer twists in ../PENCIL_SELMER.md. It does not exclude rational points or compute a Cassels–Tate pairing. The rank bound remains 1 through 3.

The proof and extension cases received independent mathematical review. Run `python kummer_obstruction/check_fibration_obstruction.py` and `python kummer_obstruction/verify_model_determinants.py` with SymPy. Run direct_labels.m in Magma to reproduce the seven exact true-Selmer label identifications; direct_labels_output.json records the completed successful run. The JSON models and determinant checker are independent exact inputs/checks. These checks complement, not replace, the geometric and cohomological proof.

No pairing matrix or new point on an unknown class has been obtained. Point-search failures establish no nonexistence theorem. Other exploratory calculations are not part of this deliverable.
