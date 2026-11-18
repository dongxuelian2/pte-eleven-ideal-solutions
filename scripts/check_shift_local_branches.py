"""Rational tangent certificates for simultaneous local shift lifts (§5.57)."""

from itertools import combinations, combinations_with_replacement

from sympy import Matrix


# Type number, half-length of the central chain, and the two ordered
# tilings of the positive half-integers by the remaining reflected pairs.
tilings = [
    (1, 1, (4,), (2, 1, 1, 1)),
    (2, 1, (1, 3), (3, 1, 1)),
    (3, 1, (3, 1), (2, 1, 2)),
    (5, 1, (2, 2), (2, 2, 1)),
    (6, 1, (1, 1, 2), (4, 1)),
    (7, 1, (1, 2, 1), (3, 2)),
    (8, 1, (1, 1, 1, 1), (5,)),
    (9, 2, (3,), (1, 2, 1, 1)),
    (10, 2, (1, 2), (1, 3, 1)),
    (11, 2, (2, 1), (1, 2, 2)),
    (12, 2, (1, 1, 1), (1, 4)),
    (14, 3, (2,), (1, 1, 2, 1)),
    (15, 3, (1, 1), (1, 1, 3)),
    (17, 4, (1,), (1, 1, 1, 2)),
    (18, 5, (), (1, 1, 1, 1, 1)),
]


def block_indices(lengths, start):
    return [j for j, length in enumerate(lengths, start) for _ in range(length)]


# The numbering in the partition table of §5.32.
patterns = [(m, aa, bb) for m in range(1, 6) for count in range(6)
            for aa in combinations_with_replacement(range(1, 6), count)
            for bb in combinations_with_replacement(range(1, 6), 5-count)
            if sum(aa) == 5-m and sum(bb) == 5]
assert len(patterns) == 18
assert len({row[0] for row in tilings}) == 15

for type_id, m, a_lengths, b_lengths in tilings:
    assert patterns[type_id-1] == (m, tuple(sorted(a_lengths)), tuple(sorted(b_lengths)))
    a_blocks = [0]*m + block_indices(a_lengths, 1)
    b_blocks = block_indices(b_lengths, len(a_lengths)+1)
    assert len(a_blocks) == len(b_blocks) == 5

    # Center j has infinitesimal displacement +v_j on its positive chain
    # and -v_j on its reflected chain. Each moment row is rescaled by a
    # common nonzero factor, using doubled positive vertex positions.
    jac = Matrix.zeros(4, 5)
    for blocks, color in ((a_blocks, 1), (b_blocks, -1)):
        for k, j in enumerate(blocks):
            if j:
                for row in range(4):
                    jac[row, j-1] += color*(2*k+1)**(2*row+1)
    minor = jac[:, :4]
    assert minor.det() != 0
    velocity = list(minor.inv()*(-jac[:, 4])) + [1]
    assert jac*Matrix(velocity) == Matrix.zeros(4, 1)
    if type_id == 3:
        assert [201*v for v in velocity] == [210, 200, 126, 246, 201]

    # Signed labels distinguish all eleven prescribed chains, and zero
    # denotes the fixed central chain. Vertex positions are still doubled.
    velocity = [0] + velocity
    points = [(sign*(2*k+1), sign*j, sign*velocity[j])
              for blocks in (a_blocks, b_blocks) for k, j in enumerate(blocks)
              for sign in (-1, 1)]
    prescribed_edges = 0
    for (x, chain_x, vx), (y, chain_y, vy) in combinations(points, 2):
        if chain_x == chain_y:
            assert x != y
            prescribed_edges += abs(x-y) == 2
        elif abs(x-y) in (0, 2):
            assert vx != vy
    assert prescribed_edges == 9

print("Nonsingular local lift branches for Types", ", ".join(str(t[0]) for t in tilings))
