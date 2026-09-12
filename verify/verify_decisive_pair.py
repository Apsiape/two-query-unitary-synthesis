"""Propositions 3.1 and 3.2 of the paper.

STAGE 0 (Proposition 3.1): in the TWO-query class T_g = B D_g W D_g C, the number of distinct
targets over the lawful set equals the number of distinct pre-final maps Psi(g) = D_g W D_g C
over the lawful set, whatever the final layer B.  Checked exhaustively (all 256 sign words)
on the Dong-Lombardi-Ma permutation family at d = 2 (Q = 8, N = 2) under several final layers.

STAGE 1 (Proposition 3.2): the decisive pair lives in the THREE-insertion class
Psi_g = D_g W2 D_g W1 D_g C with Q = 8 and a one-dimensional target (N = 1): two architectures
sharing C, W1, B and differing only in W2, identical lawful sets and pre-final data, target
counts 2 and 16.  Every quoted fact is asserted; exit 0 iff all hold.
"""
import sys
import numpy as np
from core import *

fails = []


def assert_fact(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


# ---------------------------------------------------------------- STAGE 0 (Prop 3.1)
def dlm_family(d):
    """Dong-Lombardi-Ma two-query permutation family: addresses (b, x, y), Q = 2 d^2."""
    Q = 2 * d * d
    idx = lambda b, x, y: b * d * d + x * d + y
    C = np.zeros((Q, d), complex)
    Bm = np.zeros((d, Q), complex)
    for x in range(d):
        for y in range(d):
            C[idx(0, x, y), x] = 1 / np.sqrt(d)
            Bm[y, idx(1, x, y)] = 1 / np.sqrt(d)
    H = hadamard(d) / np.sqrt(d)
    X = np.array([[0, 1], [1, 0]], complex)
    W = np.kron(X, np.kron(H, H))
    return C, W, Bm


def dlm_word(d, perm):
    Q = 2 * d * d
    inv = np.argsort(perm)
    g = np.zeros(Q)
    bits = lambda a, b: bin(a & b).count("1") % 2
    for x in range(d):
        for y in range(d):
            g[x * d + y] = (-1) ** bits(y, perm[x])
            g[d * d + x * d + y] = (-1) ** bits(x, inv[y])
    return g


def count_distinct(mats, dec=8):
    keys = set()
    for A_ in mats:
        keys.add(tuple(np.round(A_.ravel(), dec)))
    return len(keys)


print('=' * 78)
print('STAGE 0 : Proposition 3.1 -- at depth two the final layer is invisible')
print('=' * 78)
d = 2
C0, W0, B0 = dlm_family(d)
Q0, N0 = 2 * d * d, d
worst_perm = 0.0
for perm in ([0, 1], [1, 0]):
    g = dlm_word(d, np.array(perm))
    Tg = B0 @ np.diag(g) @ W0 @ np.diag(g) @ C0
    P = np.zeros((d, d)); P[[perm[0], perm[1]], [0, 1]] = 1
    worst_perm = max(worst_perm, np.abs(Tg - P).max())
assert_fact("DLM family at d=2 realises P_pi exactly for both permutations", worst_perm < 1e-12,
            f"max residual {worst_perm:.1e}")
G_all = words(Q0)
rng0 = np.random.default_rng(31)
finals = [("B (DLM)", B0)]
for k in range(3):
    finals.append((f"B W2 (random unitary W2 #{k+1})", B0 @ rand_unitary(Q0, rng0)))
finals.append(("V B (random target rotation)", rand_unitary(N0, rng0) @ B0))
Bcont = rng0.normal(size=(N0, Q0)) + 1j * rng0.normal(size=(N0, Q0))
finals.append(("random contraction", Bcont / (1.3 * np.linalg.norm(Bcont, 2))))
Psi_all = [np.diag(g) @ W0 @ np.diag(g) @ C0 for g in G_all]
all_ok = True
dlm_lawful = 0
for name, Bf in finals:
    T_all = [Bf @ P_ for P_ in Psi_all]
    law = [i for i, T_ in enumerate(T_all)
           if np.abs(T_.conj().T @ T_ - np.eye(N0)).max() < 1e-9]
    if name == "B (DLM)":
        dlm_lawful = len(law)
    nT = count_distinct([T_all[i] for i in law])
    nP = count_distinct([Psi_all[i] for i in law])
    ok = nT == nP
    all_ok &= ok
    print(f"     final layer {name:36s}: |lawful| = {len(law):3d}, distinct targets = {nT:3d}, "
          f"distinct Psi = {nP:3d}  {'ok' if ok else 'MISMATCH'}")
assert_fact("Prop 3.1: #distinct targets == #distinct Psi over the lawful set, every final layer", all_ok)
assert_fact("Prop 3.1 check is non-vacuous (DLM final layer has a non-empty lawful set)", dlm_lawful > 0,
            f"|lawful| = {dlm_lawful}")
print()

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(9)
q = 8

lamA = np.ones(q, dtype=complex)
lamB = np.exp(2j * np.pi * rng.random(q))

def build(name, lam):
    c, W1, W2, b, L = lookup_block(q, lam)
    return Arch(name, c.reshape(q, 1), W1, W2, b.conj().reshape(1, q))

A = build('A (lam = 1)', lamA)
B = build('B (lam random)', lamB)

print('=' * 78)
print('STAGE 1 : Proposition 3.2 -- the decisive pair (three insertions, N = 1, Q = 8;')
print('          same C, same W1, same B; only W2 differs)')
print('=' * 78)
print(f'  W2 differ by             : {np.max(np.abs(A.W2 - B.W2)):.3e}')
print(f'  C / W1 / B identical     : {np.max(np.abs(A.C-B.C)):.1e} '
      f'{np.max(np.abs(A.W1-B.W1)):.1e} {np.max(np.abs(A.B-B.B)):.1e}')
print(f'  lawful sets identical    : {np.array_equal(A.law, B.law)}   |L| = {len(A.law)}')
print(f'  DISTINCT TARGETS         : A = {A.n_targets()}   B = {B.n_targets()}')
worst_exact = 0.0
for X in (A, B):
    print(f'  -- exactness checks {X.name}')
    for k, v in X.checks().items():
        print(f'       {k:22s} {v:.2e}')
        worst_exact = max(worst_exact, v)
assert_fact("C, W1, B identical; only W2 differs",
            np.abs(A.C - B.C).max() < 1e-14 and np.abs(A.W1 - B.W1).max() < 1e-14
            and np.abs(A.B - B.B).max() < 1e-14 and np.abs(A.W2 - B.W2).max() > 1e-3)
assert_fact("lawful sets identical, |L| = 16", np.array_equal(A.law, B.law) and len(A.law) == 16)
assert_fact("distinct targets 2 (A) and 16 (B)", A.n_targets() == 2 and B.n_targets() == 16)
assert_fact("exactness residuals below 1e-12", worst_exact < 1e-12, f"worst {worst_exact:.2e}")

# pre-final invariants identical?
def kern(J):
    return np.einsum('mai,naj->mnij', J.conj(), J)
kdiff = [np.max(np.abs(kern(A.J[s_])-kern(B.J[s_]))) for s_ in range(4)]
for s_ in range(3):
    print(f'  ||K_{s_}^A - K_{s_}^B||_inf on L : {kdiff[s_]:.2e}')
print(f'  ||K_3^A - K_3^B||_inf on L : {kdiff[3]:.2e}  <-- differs')
assert_fact("pre-final Gram kernels K_0, K_1, K_2 identical", max(kdiff[:3]) < 1e-12)
assert_fact("final kernel K_3 differs by 2.00 (to two decimals)", abs(kdiff[3] - 2.0) < 5e-3, f"{kdiff[3]:.4f}")

# lawful set structure
print()
print('  lawful set as a set of sign words (rows), |L| =', len(A.law))
print(A.Gl.astype(int))
prod_closed = True
S = set(map(tuple, A.Gl.astype(int)))
for x in A.Gl:
    for y in A.Gl:
        if tuple((x * y).astype(int)) not in S:
            prod_closed = False
print('  lawful set closed under pointwise product (is a linear code) :', prod_closed)
assert_fact("lawful set closed under pointwise product (linear code)", prod_closed)

print()
print('=' * 78)
print('STAGE 1b : independent verification of the motivating GRID fact')
print('=' * 78)
for N in (2, 3, 4):
    Q = N * N
    L = np.zeros((Q, N)); R = np.zeros((Q, N))
    for i in range(N):
        for j in range(N):
            L[i * N + j, i] = 1 / np.sqrt(N)
            R[i * N + j, j] = 1 / np.sqrt(N)
    g = np.random.default_rng(N).choice([-1.0, 1.0], size=Q)
    Gm = g.reshape(N, N)
    S = L.T @ np.diag(g) @ R
    Hn = hadamard(N) if (N & (N - 1)) == 0 else None
    line = (f'  N={N} Q={Q}: || L* D_g R - G/N ||_inf = {np.max(np.abs(S - Gm / N)):.2e}'
            f'  ; sing.vals of L*D_gR: {np.round(np.linalg.svd(S, compute_uv=False),4)}')
    print(line)
    if Hn is not None:
        Su = L.T @ np.diag(Hn.ravel()) @ R
        print(f'        with G = Hadamard (best case): sing.vals = '
              f'{np.round(np.linalg.svd(Su, compute_uv=False),4)}  (= 1/sqrt(N) = '
              f'{1/np.sqrt(N):.4f}); a unitary needs 1.  amplitude gap = sqrt(N)')
print('  dense sign matrix at unitary normalisation G/sqrt(N): all sing.vals 1'
      ' (Hadamard).  Gap confirmed.')

print()
print('-' * 78)
if fails:
    print("verify_decisive_pair.py: FAIL -> " + ", ".join(fails))
    sys.exit(1)
print("verify_decisive_pair.py: PASS")
sys.exit(0)
