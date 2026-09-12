"""
Verifier for the two-query width bound in its audited real-selector form (paper, Corollary 4.2).

Checks the four lines the theorem actually is:

  (i)   Re Tr(G^* T_g)  =  g^T Herm(M) g          [exact, M independent of g]
  (ii)  =  v^* Herm(K) v,   v = D^{1/2} g,  ||v||^2 = Tr D <= 2 min(N,Q)
  (iii) =>  |Re Tr(G^* T_g)| <= Tr(D) ||K||_op    POINTWISE IN G, every g
  (iv)  both square functions of the Gaussian series K are <= I

with the SYMMETRIC weighting  d_p = lambda_p + mu_p,  lambda_p = ||b_p||^2,
mu_p = ||c_p||^2.  That weighting is the theorem: either half alone satisfies
one square function and breaks the other, which is checked explicitly below.

Also checks Proposition 2.2 (lawfulness, isometric case), the zero-leverage lemma, and
Remark 4.3 (the one-sided weightings fail on a substantial fraction of instances).

Deterministic (seed 20260821), numpy only, exit 0 iff all pass.
"""
import numpy as np
import itertools

SEED = 20260821
rng = np.random.default_rng(SEED)
fails = []


def check(name, residual, tol=1e-11):
    ok = residual <= tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:54s} {residual:.3e}")
    if not ok:
        fails.append(name)


def haar(n):
    z = (rng.standard_normal((n, n)) + 1j*rng.standard_normal((n, n)))/np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q*(np.diag(r)/abs(np.diag(r)))


def contraction(m, n, scale=0.9):
    A = (rng.standard_normal((m, n)) + 1j*rng.standard_normal((m, n)))/np.sqrt(2)
    return scale*A/np.linalg.norm(A, 2)


def herm(X):
    return (X + X.conj().T)/2


def build(B, W, C, G):
    """M with Tr(G^* T_g) = g^T M g.

    Tr(G^* T_g) = sum_{pq} g_p g_q W_pq * (b_p^T conj(G) c_q),
    with b_p = B[:,p] and c_q = C[q,:].  Hence
        M = W .* (B^T conj(G) C^T).
    """
    return W * (B.T @ G.conj() @ C.T)


def leverages(B, C):
    lam = np.einsum('ip,ip->p', B.conj(), B).real      # ||b_p||^2
    mu = np.einsum('qj,qj->q', C.conj(), C).real       # ||c_q||^2
    return lam, mu


print(__doc__.strip())
print(f"\nseed = {SEED}\n")

# ---------------------------------------------------------- Prop 2.2
print("Proposition 2.2 -- lawfulness criterion (isometric C and W)")
worst = 0.0
for _ in range(150):
    Q, N = int(rng.integers(3, 10)), int(rng.integers(1, 4))
    N = min(N, Q-1)
    C, B, W = haar(Q)[:, :N], haar(Q)[:, :N].conj().T, haar(Q)
    g = rng.choice([-1.0, 1.0], size=Q)
    Dg = np.diag(g.astype(complex))
    Psi = Dg @ W @ Dg @ C
    Pi = np.eye(Q) - B.conj().T @ B
    T = B @ Psi
    worst = max(worst, np.linalg.norm(
        T.conj().T @ T - (np.eye(N) - Psi.conj().T @ Pi @ Psi), 'fro'))
check("T*T = I - Psi^* Pi Psi", worst)

# ---------------------------------------------------------- (i)
print("\n(i) the exact expansion -- an identity, M independent of g")
worst = 0.0
for _ in range(150):
    Q, N = int(rng.integers(2, 8)), int(rng.integers(1, 5))
    B, W, C = contraction(N, Q), contraction(Q, Q), contraction(Q, N)
    G = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2)
    M = build(B, W, C, G)
    for _ in range(4):
        g = rng.choice([-1.0, 1.0], size=Q)
        Dg = np.diag(g.astype(complex))
        lhs = np.real(np.trace(G.conj().T @ (B @ Dg @ W @ Dg @ C)))
        worst = max(worst, abs(lhs - np.real(g @ herm(M) @ g)))
check("Re Tr(G^* T_g) = g^T Herm(M) g", worst)

# ---------------------------------------------------------- (ii)+(iii)
print("\n(ii)-(iii) pointwise absorption of the 2^Q sign supremum")
print("           v = D^{1/2} g,  Tr D <= 2 min(N,Q);  exhaustive over the cube")
worst_viol, worst_trD, tight = 0.0, 0.0, np.inf
for _ in range(40):
    Q, N = int(rng.integers(2, 9)), int(rng.integers(1, 5))
    B, W, C = contraction(N, Q), contraction(Q, Q), contraction(Q, N)
    G = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2)
    M = build(B, W, C, G)
    lam, mu = leverages(B, C)
    d = lam + mu
    supp = d > 1e-13
    if supp.sum() == 0:
        continue
    Dh = np.diag(np.where(supp, d, 1.0)**-0.5)
    K = Dh @ herm(M) @ Dh
    bound = d.sum()*np.linalg.norm(herm(K), 2)
    sup = max(abs(np.real(np.array(g) @ herm(M) @ np.array(g)))
              for g in itertools.product([-1.0, 1.0], repeat=Q))
    worst_viol = max(worst_viol, max(0.0, sup - bound))
    worst_trD = max(worst_trD, max(0.0, d.sum() - 2*min(N, Q)))
    tight = min(tight, bound/max(sup, 1e-300))
check("|Re Tr(G^* T_g)| <= Tr(D) ||Herm(K)||_op, all g", worst_viol)
check("Tr D <= 2 min(N,Q)", worst_trD)
print(f"           tightest bound/sup ratio: {tight:.3f}   (>= 1 required)")

# ---------------------------------------------------------- (iv)
print("\n(iv) BOTH square functions <= I  -- and the symmetric weighting is why")
print("     S1 = D^{-1/2}[conj(B^*B) o (W diag(mu/d) W^*)]D^{-1/2}")
print("     S2 = D^{-1/2}[conj(CC^*) o (W^* diag(lam/d) W)]D^{-1/2}")
w_sym, w_lam_only, w_mu_only = 0.0, 0.0, 0.0
n_inst = n_lam_break = n_mu_break = 0
for _ in range(200):
    Q, N = int(rng.integers(2, 10)), int(rng.integers(1, 5))
    B, W, C = contraction(N, Q), contraction(Q, Q), contraction(Q, N)
    lam, mu = leverages(B, C)

    def squares(d):
        supp = d > 1e-13
        if supp.sum() == 0:
            return None
        dd = np.where(supp, d, 1.0)
        Dh = np.diag(dd**-0.5)
        GB = np.conj(B.conj().T @ B)          # conj(B^*B), PSD
        GC = np.conj(C @ C.conj().T)          # conj(CC^*), PSD
        S1 = Dh @ (GB * (W @ np.diag(mu/dd) @ W.conj().T)) @ Dh
        S2 = Dh @ (GC * (W.conj().T @ np.diag(lam/dd) @ W)) @ Dh
        e1 = np.linalg.eigvalsh(herm(S1)).max()
        e2 = np.linalg.eigvalsh(herm(S2)).max()
        return e1, e2

    r = squares(lam + mu)                      # the symmetric weighting
    if r:
        w_sym = max(w_sym, max(0.0, r[0]-1.0), max(0.0, r[1]-1.0))
    n_inst += 1
    r = squares(lam.copy())                    # lambda alone
    if r:
        w_lam_only = max(w_lam_only, max(r[0]-1.0, r[1]-1.0))
        if max(r[0]-1.0, r[1]-1.0) > 1e-6:
            n_lam_break += 1
    r = squares(mu.copy())                     # mu alone
    if r:
        w_mu_only = max(w_mu_only, max(r[0]-1.0, r[1]-1.0))
        if max(r[0]-1.0, r[1]-1.0) > 1e-6:
            n_mu_break += 1
check("d = lambda + mu : BOTH square functions <= I", w_sym, tol=1e-9)
print(f"     d = lambda alone : worst excess over I = {w_lam_only:+.3f}"
      f"   {'(BREAKS)' if w_lam_only > 1e-6 else ''}")
print(f"     d = mu     alone : worst excess over I = {w_mu_only:+.3f}"
      f"   {'(BREAKS)' if w_mu_only > 1e-6 else ''}")
frac_lam, frac_mu = n_lam_break / n_inst, n_mu_break / n_inst
print(f"     Remark 4.3: lambda alone breaks on {frac_lam:.1%}, mu alone on {frac_mu:.1%} "
      f"of {n_inst} instances (paper: 'a substantial fraction')")
check("Remark 4.3: lambda alone fails domination on >= 10% of instances", max(0.0, 0.10 - frac_lam), tol=0.0)
check("Remark 4.3: mu alone fails domination on >= 10% of instances", max(0.0, 0.10 - frac_mu), tol=0.0)

# ---------------------------------------------------------- zero leverage
print("\nZero-leverage lemma -- d_p = 0 kills BOTH row p and column p of M")
worst = 0.0
for _ in range(100):
    Q, N = int(rng.integers(3, 9)), int(rng.integers(1, 4))
    B, W, C = contraction(N, Q), contraction(Q, Q), contraction(Q, N)
    p0 = int(rng.integers(0, Q))
    B[:, p0] = 0.0                              # forces lambda_{p0} = 0
    C[p0, :] = 0.0                              # forces mu_{p0} = 0
    G = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2)
    M = build(B, W, C, G)
    lam, mu = leverages(B, C)
    assert lam[p0] + mu[p0] < 1e-24
    worst = max(worst, np.linalg.norm(M[p0, :]) + np.linalg.norm(M[:, p0]))
check("d_p = 0  =>  row p and column p of M vanish", worst)

print("\n" + "="*74)
if fails:
    print(f"FAILED: {len(fails)} check(s): " + "; ".join(fails))
    raise SystemExit(1)
print("ALL CHECKS PASSED")
raise SystemExit(0)
