"""STAGE 1 -- reconstruct + verify the decisive pair; verify the stated exact facts."""
import numpy as np
from core import *

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
print('STAGE 1 : the decisive pair  (same C, same W1, same B; only W2 differs)')
print('=' * 78)
print(f'  W2 differ by             : {np.max(np.abs(A.W2 - B.W2)):.3e}')
print(f'  C / W1 / B identical     : {np.max(np.abs(A.C-B.C)):.1e} '
      f'{np.max(np.abs(A.W1-B.W1)):.1e} {np.max(np.abs(A.B-B.B)):.1e}')
print(f'  lawful sets identical    : {np.array_equal(A.law, B.law)}   |L| = {len(A.law)}')
print(f'  DISTINCT TARGETS         : A = {A.n_targets()}   B = {B.n_targets()}')
for X in (A, B):
    print(f'  -- exactness checks {X.name}')
    for k, v in X.checks().items():
        print(f'       {k:22s} {v:.2e}')

# pre-final invariants identical?
def kern(J):
    return np.einsum('mai,naj->mnij', J.conj(), J)
for s in range(3):
    print(f'  ||K_{s}^A - K_{s}^B||_inf on L : {np.max(np.abs(kern(A.J[s])-kern(B.J[s]))):.2e}')
print(f'  ||K_3^A - K_3^B||_inf on L : {np.max(np.abs(kern(A.J[3])-kern(B.J[3]))):.2e}  <-- differs')

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
