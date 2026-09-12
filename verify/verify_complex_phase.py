"""Theorem 4.1 (independent unimodular selectors) -- machine checks.

Architecture: T_{x,y} = B D_x W D_y C with B (N x Q), W (Q x Q), C (Q x N) arbitrary
contractions, D_x = diag(x), D_y = diag(y), x, y in the torus T^Q (|x_p| = |y_q| = 1).
The repeated real sign word is the special case x = y = g in {+-1}^Q.

Notation (paper, Section 4): b_p = B e_p (p-th column of B), c_q = e_q^T C (q-th row of C),
leverages lam_p = ||b_p||^2, mu_q = ||c_q||^2, Lam = diag(lam), Mu = diag(mu),
P_lam, P_mu the projections onto their supports.

Checks (all deterministic, seed fixed):
 (a) exact expansion   Tr(G^* T_{x,y}) = x^T M y  with  M_{pq} = W_{pq} (c_q G^* b_p),
     M independent of x and y -- for arbitrary complex x, y;
 (b) pointwise absorption  |x^T M y| <= ||B||_F ||C||_F ||K||_op  for unimodular x, y,
     where K = Lam^{-1/2} M Mu^{-1/2} (zero-leverage coordinates dropped); tested on random
     torus points AND at the torus supremum found by alternating maximisation;
 (c) the two variance operators of the matrix Gaussian series K have the closed forms
        E K K^* = Lam^{-1/2} [ conj(B^* B) o (W P_mu W^*) ] Lam^{-1/2},
        E K^* K = Mu^{-1/2}  [ conj(C C^*) o (W^* P_lam W) ] Mu^{-1/2},
     agree with Monte Carlo, and are each dominated by the corresponding support projection
     (exact eigenvalue check on the closed forms);
 (d) the real-selector case: |Re Tr(G^* T_g)| <= |g^T M g| for every sign word g;
 (e) the three-term unitarity identity (Proposition 2.2) for contractive C and W:
        T^* T = I - (I - C^* C) - C^* D (I - W^* W) D C - Psi^* Pi Psi,
        Psi = D W D C,  Pi = I - B^* B,
     and the failure of the shorter forms when C^* C != I or W^* W != I;
 (f) the empirical E||K|| / sqrt(log 2Q) against the proved cap sqrt(2) (Tropp's bound for
     the pooled real series of 2N^2 terms with variance parameter <= 1; the paper quotes
     the looser 2 obtained by the triangle inequality).

Exit 0 iff every check passes.  numpy only.
"""
import sys
import numpy as np

rng = np.random.default_rng(20260912)
TOL = 1e-9


def contraction(m, n, scale=1.3):
    A = rng.normal(size=(m, n)) + 1j * rng.normal(size=(m, n))
    return A / (np.linalg.norm(A, 2) * scale)


def gaussian(N):
    """Standard complex Gaussian: E|G_ij|^2 = 1."""
    return (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)


def inv_sqrt_diag(v):
    return np.diag([1 / np.sqrt(t) if t > 1e-14 else 0.0 for t in v])


def support_proj(v):
    return np.diag((v > 1e-14).astype(float))


def unimodular(Q):
    return np.exp(1j * rng.uniform(0, 2 * np.pi, Q))


fails = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


print("=" * 74)
print("verify_complex_phase.py -- Theorem 4.1 with independent unimodular selectors")
print("=" * 74)

# ---------------------------------------------------------------- (a), (b), (c), (d)
max_expansion_dev = 0.0
max_mc_dev = 0.0
absorption_violations = 0
absorption_checks = 0
worst_absorption_ratio = 0.0
sup_violations = 0
worst_sup_ratio = 0.0


def torus_sup(M, iters=40):
    """Alternating maximisation of |x^T M y| over the torus: for fixed y the best x is the
    conjugate phase of M y (value sum_p |(M y)_p|), and symmetrically for y."""
    Q1, Q2 = M.shape
    y = unimodular(Q2)
    best = 0.0
    for _ in range(iters):
        v = M @ y
        x = np.exp(-1j * np.angle(v))
        u = M.T @ x
        y = np.exp(-1j * np.angle(u))
        best = max(best, abs(x @ M @ y))
    return best
real_case_violations = 0
domination_ok = True

trials = 40
mc_samples = 6000
for trial in range(trials):
    N = int(rng.integers(2, 6))
    Q = int(rng.integers(3, 9))
    B, W, C = contraction(N, Q), contraction(Q, Q), contraction(Q, N)
    if trial % 3 == 0:                      # zero-leverage coordinates on both sides
        B[:, 0] = 0
        C[1, :] = 0
    lam = np.sum(np.abs(B) ** 2, axis=0)
    mu = np.sum(np.abs(C) ** 2, axis=1)
    Li, Mi = inv_sqrt_diag(lam), inv_sqrt_diag(mu)
    Pl, Pm = support_proj(lam), support_proj(mu)
    fro = np.linalg.norm(B, "fro") * np.linalg.norm(C, "fro")

    S1 = np.zeros((Q, Q), complex)
    S2 = np.zeros((Q, Q), complex)
    for s in range(mc_samples):
        G = gaussian(N)
        M = W * (C @ G.conj().T @ B).T          # M_{pq} = W_{pq} (c_q G^* b_p)
        K = Li @ M @ Mi
        S1 += K @ K.conj().T / mc_samples
        S2 += K.conj().T @ K / mc_samples
        if s < 4:
            # (a) exact expansion for arbitrary complex x, y
            x = rng.normal(size=Q) + 1j * rng.normal(size=Q)
            y = rng.normal(size=Q) + 1j * rng.normal(size=Q)
            lhs = np.trace(G.conj().T @ B @ np.diag(x) @ W @ np.diag(y) @ C)
            max_expansion_dev = max(max_expansion_dev, abs(lhs - x @ M @ y))
            # (b) pointwise absorption on the torus
            Kop = np.linalg.norm(K, 2)
            for _ in range(60):
                xu, yu = unimodular(Q), unimodular(Q)
                val = abs(xu @ M @ yu)
                bound = fro * Kop
                absorption_checks += 1
                if val > bound + 1e-9:
                    absorption_violations += 1
                if bound > 1e-12:
                    worst_absorption_ratio = max(worst_absorption_ratio, val / bound)
            sup_val = max(torus_sup(M) for _ in range(3))
            if sup_val > fro * Kop + 1e-9:
                sup_violations += 1
            if fro * Kop > 1e-12:
                worst_sup_ratio = max(worst_sup_ratio, sup_val / (fro * Kop))
            # (d) real selectors: |Re Tr(G^* T_g)| <= |g^T M g|
            g = rng.choice([1.0, -1.0], Q)
            Tg = B @ np.diag(g) @ W @ np.diag(g) @ C
            if abs(np.real(np.trace(G.conj().T @ Tg))) > abs(g @ M @ g) + 1e-9:
                real_case_violations += 1

    # (c) closed forms and domination
    S1f = Li @ (np.conj(B.conj().T @ B) * (W @ Pm @ W.conj().T)) @ Li
    S2f = Mi @ (np.conj(C @ C.conj().T) * (W.conj().T @ Pl @ W)) @ Mi
    max_mc_dev = max(max_mc_dev, np.abs(S1 - S1f).max(), np.abs(S2 - S2f).max())
    dom1 = np.linalg.eigvalsh((Pl - S1f + (Pl - S1f).conj().T) / 2).min()
    dom2 = np.linalg.eigvalsh((Pm - S2f + (Pm - S2f).conj().T) / 2).min()
    domination_ok &= dom1 > -TOL and dom2 > -TOL

check("(a) exact expansion Tr(G^* T_{x,y}) = x^T M y, complex x,y",
      max_expansion_dev < 1e-9, f"max deviation {max_expansion_dev:.2e}")
check("(b) pointwise absorption at the torus supremum (alternating maximisation)",
      sup_violations == 0, f"0 violations expected, worst sup/bound {worst_sup_ratio:.4f}")
check("(b) pointwise absorption |x^T M y| <= ||B||_F ||C||_F ||K|| on the torus",
      absorption_violations == 0,
      f"{absorption_checks} checks, 0 violations expected, worst ratio {worst_absorption_ratio:.3f}")
check("(c) closed-form variance operators agree with Monte Carlo",
      max_mc_dev < 3e-2, f"max deviation {max_mc_dev:.2e} at MC noise ~1e-2")
check("(c) both closed-form variance operators <= support projections (exact)",
      domination_ok)
check("(d) real selectors are the special case x = y = g",
      real_case_violations == 0)

# ---------------------------------------------------------------- (e)
N, Q = 3, 6
Cc = contraction(Q, N)                       # strictly contractive: C^* C != I
Bc = contraction(N, Q)
Wc = contraction(Q, Q)                       # strictly contractive: W^* W != I
g = rng.choice([1.0, -1.0], Q)
D = np.diag(g)
T = Bc @ D @ Wc @ D @ Cc
Psi = D @ Wc @ D @ Cc
Pi = np.eye(Q) - Bc.conj().T @ Bc
lhs = T.conj().T @ T
term_C = np.eye(N) - Cc.conj().T @ Cc
term_W = Cc.conj().T @ D @ (np.eye(Q) - Wc.conj().T @ Wc) @ D @ Cc
term_B = Psi.conj().T @ Pi @ Psi
rhs_three = np.eye(N) - term_C - term_W - term_B
rhs_two = np.eye(N) - term_W - term_B         # the earlier draft's form (needs C^* C = I)
rhs_one = np.eye(N) - term_B                  # needs C^* C = I and W^* W = I
psd_ok = all(np.linalg.eigvalsh((X + X.conj().T) / 2).min() > -1e-10
             for X in (term_C, term_W, term_B))
check("(e) three-term T^*T identity holds for contractive C, W (Proposition 2.2)",
      np.allclose(lhs, rhs_three, atol=1e-10))
check("(e) all three subtracted terms are positive semidefinite", psd_ok)
check("(e) the two-term form fails when C^*C != I (as it must)",
      not np.allclose(lhs, rhs_two, atol=1e-6),
      f"max discrepancy {np.abs(lhs - rhs_two).max():.2e}")
check("(e) the one-term form fails when W^*W != I (as it must)",
      not np.allclose(lhs, rhs_one, atol=1e-6),
      f"max discrepancy {np.abs(lhs - rhs_one).max():.2e}")
# the isometric special case reduces to the one-term criterion
Ciso, _ = np.linalg.qr(rng.normal(size=(Q, N)) + 1j * rng.normal(size=(Q, N)))
Wu, _ = np.linalg.qr(rng.normal(size=(Q, Q)) + 1j * rng.normal(size=(Q, Q)))
Psi_i = D @ Wu @ D @ Ciso
T_i = Bc @ Psi_i
check("(e) isometric C, W: T^*T = I - Psi^* Pi Psi exactly",
      np.allclose(T_i.conj().T @ T_i, np.eye(N) - Psi_i.conj().T @ Pi @ Psi_i, atol=1e-10))

# ---------------------------------------------------------------- (f)
print("  (f) empirical E||K|| / sqrt(log 2Q), proved cap sqrt(2) (pooled real series, v <= 1):")
ratio_ok = True
for (N, Q) in [(4, 8), (8, 16), (8, 32)]:
    Ciso, _ = np.linalg.qr(rng.normal(size=(Q, N)) + 1j * rng.normal(size=(Q, N)))
    Bco = Ciso.conj().T
    Wu, _ = np.linalg.qr(rng.normal(size=(Q, Q)) + 1j * rng.normal(size=(Q, Q)))
    lam = np.sum(np.abs(Bco) ** 2, axis=0)
    mu = np.sum(np.abs(Ciso) ** 2, axis=1)
    Li, Mi = inv_sqrt_diag(lam), inv_sqrt_diag(mu)
    norms = []
    for _ in range(300):
        G = gaussian(N)
        K = Li @ (Wu * (Ciso @ G.conj().T @ Bco).T) @ Mi
        norms.append(np.linalg.norm(K, 2))
    ratio = np.mean(norms) / np.sqrt(np.log(2 * Q))
    ratio_ok &= ratio <= np.sqrt(2.0)
    print(f"      N={N:2d} Q={Q:2d}: E||K|| = {np.mean(norms):.3f}, ratio = {ratio:.3f}")
check("(f) empirical ratio within the proved cap sqrt(2) (hence within the quoted 2)", ratio_ok)

print("-" * 74)
if fails:
    print("verify_complex_phase.py: FAIL -> " + ", ".join(fails))
    sys.exit(1)
print("verify_complex_phase.py: PASS")
sys.exit(0)
