"""Corollary 3 (address length) -- the compression to address-respecting subspaces.

Setting: a general two-query architecture whose oracle acts through address projectors
of ARBITRARY rank,  D_x = sum_{j <= M} x_j P_j  (P_j mutually orthogonal projectors on
the Q-dimensional register, x in T^M), so that

    T_{x,y} = B D_x W D_y C.

Claim verified: with  L = sum_j ran(P_j B^*)  and  R = sum_j ran(P_j C)  (each of dimension
at most M N) and isometries V_L, V_R onto them, the target-independent factorization

    T_{x,y} = (B V_L) Dt_x (V_L^* W V_R) Dt_y (V_R^* C),
    Dt_x = sum_j x_j V_L^* P_j V_L  (a diagonal phase matrix in an address-respecting basis),

holds exactly, with  ||B V_L||_F = ||B||_F,  ||V_R^* C||_F = ||C||_F,  and the middle
V_L^* W V_R a contraction.  Theorem 1 then applies with Q replaced by max(dim L, dim R) <= M N,
which is how the address-length bound is obtained.

Exit 0 iff every check passes.  numpy only.
"""
import sys
import numpy as np

rng = np.random.default_rng(20260911)


def rand_iso(m, n):
    Qm, _ = np.linalg.qr(rng.normal(size=(m, n)) + 1j * rng.normal(size=(m, n)))
    return Qm


def address_basis(mats):
    """Orthonormal basis of  sum_j ran(A_j)  respecting the address decomposition."""
    cols, labels = [], []
    for j, A in enumerate(mats):
        u, sv, _ = np.linalg.svd(A, full_matrices=False)
        k = int((sv > 1e-10).sum())
        cols.append(u[:, :k])
        labels += [j] * k
    return np.concatenate(cols, axis=1), np.array(labels)


maxdev = 0.0
dims_ok = frob_ok = contr_ok = orth_ok = True
trials = 80
for trial in range(trials):
    N = int(rng.integers(2, 5))
    M = int(rng.integers(2, 6))
    ranks = rng.integers(1, 7, size=M)
    Q = int(ranks.sum())
    Uq = rand_iso(Q, Q)
    P, s = [], 0
    for r in ranks:                              # address projectors of arbitrary rank
        Vb = Uq[:, s:s + r]
        P.append(Vb @ Vb.conj().T)
        s += r
    C = rand_iso(Q, N)
    B = rand_iso(Q, N).conj().T
    if trial % 4 == 0:
        B[:, :ranks[0]] = 0                      # a dead address block on the output side
    W = rng.normal(size=(Q, Q)) + 1j * rng.normal(size=(Q, Q))
    W /= np.linalg.norm(W, 2) * 1.2              # strictly contractive middle

    VL, labL = address_basis([Pj @ B.conj().T for Pj in P])
    VR, labR = address_basis([Pj @ C for Pj in P])
    dims_ok &= VL.shape[1] <= M * N and VR.shape[1] <= M * N
    orth_ok &= np.allclose(VL.conj().T @ VL, np.eye(VL.shape[1]), atol=1e-10)
    orth_ok &= np.allclose(VR.conj().T @ VR, np.eye(VR.shape[1]), atol=1e-10)
    Wt = VL.conj().T @ W @ VR
    contr_ok &= np.linalg.norm(Wt, 2) <= 1 + 1e-12
    frob_ok &= abs(np.linalg.norm(B @ VL, "fro") - np.linalg.norm(B, "fro")) < 1e-9
    frob_ok &= abs(np.linalg.norm(VR.conj().T @ C, "fro") - np.linalg.norm(C, "fro")) < 1e-9
    for _ in range(8):
        x = np.exp(1j * rng.uniform(0, 2 * np.pi, M))
        y = np.exp(1j * rng.uniform(0, 2 * np.pi, M))
        Dx = sum(x[j] * P[j] for j in range(M))
        Dy = sum(y[j] * P[j] for j in range(M))
        T = B @ Dx @ W @ Dy @ C
        Tc = (B @ VL) @ np.diag(x[labL]) @ Wt @ np.diag(y[labR]) @ (VR.conj().T @ C)
        maxdev = max(maxdev, np.abs(T - Tc).max())

print("=" * 74)
print("verify_address_compression.py -- Corollary 3 factorization")
print("=" * 74)
ok = True
for name, cond, detail in [
    ("exact factorization T = (B V_L) Dt_x (V_L^* W V_R) Dt_y (V_R^* C)", maxdev < 1e-12,
     f"max |T - T_compressed| = {maxdev:.1e} over {trials} architectures"),
    ("address bases orthonormal", orth_ok, ""),
    ("dim L, dim R <= M N", dims_ok, ""),
    ("outer Frobenius norms preserved", frob_ok, ""),
    ("compressed middle is a contraction", contr_ok, ""),
]:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    ok &= cond
print("-" * 74)
print("verify_address_compression.py: " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
