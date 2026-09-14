"""Exact rational controls, independent of the floating-point helper library.

Finite tests only: not a proof of the width or transport asymptotics.
Runs in bounded space using Python's standard library.
"""
from fractions import Fraction as F
from itertools import product


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


# Proposition 3.2: all 256 sign words, two fully specified architectures.
# Work with the unnormalised all-ones vector; divide scalar outputs by 8.
bits = list(product((0, 1), repeat=3))
idx = {x: i for i, x in enumerate(bits)}
add = lambda x, y: tuple(a ^ b for a, b in zip(x, y))
A = lambda x: (x[2], x[0] ^ x[2], x[1])
chars = [[(-1) ** (dot(x, a) % 2) for x in bits] for a in bits]
assert len({A(x) for x in bits}) == 8
assert len({add(x, A(x)) for x in bits}) == 8
Y = [chars[idx[A(a)]] for a in bits]
X = [chars[idx[add(a, A(a))]] for a in bits]
W1 = [[F(sum(Y[a][i] * chars[a][j] for a in range(8)), 8)
       for j in range(8)] for i in range(8)]
phases = [(F(1 - j*j, 1 + j*j), F(2*j, 1 + j*j)) for j in range(8)]
assert all(re*re + im*im == 1 for re, im in phases)
outcomes, lawful_sets = [], []
for lam in [[(F(1), F(0))] * 8, phases]:
    W2r = [[sum(lam[a][0] * chars[a][i] * X[a][j] for a in range(8))/8
            for j in range(8)] for i in range(8)]
    W2i = [[sum(lam[a][1] * chars[a][i] * X[a][j] for a in range(8))/8
            for j in range(8)] for i in range(8)]
    lawful, values = set(), set()
    for g in product((-1, 1), repeat=8):
        v = [g[i] * dot(W1[i], g) for i in range(8)]
        re = [g[i] * dot(W2r[i], v) for i in range(8)]
        im = [g[i] * dot(W2i[i], v) for i in range(8)]
        if len(set(zip(re, im))) == 1:
            lawful.add(g)
            values.add((re[0], im[0]))
    lawful_sets.append(lawful)
    outcomes.append(values)
code = {tuple(s*x for x in c) for c in chars for s in (-1, 1)}
assert lawful_sets[0] == lawful_sets[1] == code
assert [len(v) for v in outcomes] == [2, 16]
print('PASS: exact rational three-insertion control, all 256 words, 2 vs 16 scalar matrices')

# A two-query packet can have a degree-zero observable: g^2=1.
# The vanishing-linear-functional lemma cannot imply a degree-4 lower bound.
for g in product((-1, 1), repeat=2):
    packet = [x*x for x in g]
    assert packet == [1, 1]
print('PASS: query cancellation negative control for the withdrawn degree barrier')

# Exact N=2, J=2 Gaussian-grid directions, up to the common normalisation.
# Multiplication of either coordinate by i and coordinate swap are grid permutations.
# Verify P E = E U entrywise without numerical tolerance or an input sample.
grid = list(product((-1, 1), repeat=4))
grididx = {x: i for i, x in enumerate(grid)}
for mode in ('phase0', 'swap'):
    perm = []
    for x in grid:
        y = (-x[2], x[1], x[0], x[3]) if mode == 'phase0' else (x[1], x[0], x[3], x[2])
        perm.append(grididx[y])
        row_a = [(x[0], -x[2]), (x[1], -x[3])]
        row_b = [(y[0], -y[2]), (y[1], -y[3])]
        row_b_U = [(-row_b[0][1], row_b[0][0]), row_b[1]] if mode == 'phase0' else row_b[::-1]
        assert row_a == row_b_U
    assert sorted(perm) == list(range(16))
    inv = [perm.index(z) for z in range(16)]
    for z in range(16):
        # The Boolean inner-product oracle plus Hadamards is BV XOR lookup.
        for out in range(16):
            amp = sum((-1)**(((perm[z] ^ out) & y).bit_count() % 2) for y in range(16))
            assert amp == (16 if out == perm[z] else 0)
        assert z ^ inv[perm[z]] == 0
print('PASS: exact non-scalar N=2 transport and two Boolean lookups, all basis inputs')
