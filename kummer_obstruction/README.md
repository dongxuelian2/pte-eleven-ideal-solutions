# Twisted Kummer fibration obstruction

The proof in KUMMER_FIBRATION_OBSTRUCTION.md excludes genus-one fibrations over Q, finite odd-degree extensions, and all extensions of degree less than ten for all seven nonzero Selmer twists in ../proofs/splitting-pencil.md. It does not exclude rational points or compute a Cassels–Tate pairing. The rank bound remains 1 through 3.

The degree bound is sharp. A subgroup of index ten partitions the nodes into two eight-element orbits; over its fixed field, `2H−sum(E_i)` for either orbit is a nonzero square-zero divisor. The characteristic-zero K3 divisor criterion gives a genus-one fibration there. For the seven Selmer classes, local solubility and the conic Hasse principle make its base P¹. Thus their minimum extension degree is exactly ten even with a split base required; a section is not asserted. For every nonzero true translation twist of a genus-two Jacobian with a rational Weierstrass point and quintic Galois group S5, the minimum degree remains ten when a possibly nonsplit genus-zero base is allowed.

The proof and extension cases received independent mathematical review. Run `python kummer_obstruction/check_fibration_obstruction.py` and `python kummer_obstruction/verify_model_determinants.py` with SymPy. Run direct_labels.m in Magma to reproduce the seven exact true-Selmer label identifications; direct_labels_output.json records the completed successful run. The JSON models and determinant checker are independent exact inputs/checks. These checks complement, not replace, the geometric and cohomological proof.

No pairing matrix or new point on an unknown class has been obtained. Point-search failures establish no nonexistence theorem. Other exploratory calculations are not part of this deliverable.
