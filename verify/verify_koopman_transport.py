"""Theorem 4 (matching upper bound) -- the permutation-transport construction, small cases.

The construction (paper, Section 7; notes/permutation-transport.md):
  * quantize N(0,1) into K equiprobable intervals; q_k = conditional mean on interval k;
    s^2 = mean_k q_k^2;  e^2 = 1 - s^2 = mean-square quantization error;
  * addresses a in [K]^{2N}, M = K^{2N};  w_a = (q_{a_k} + i q_{a_{N+k}})_{k} / (s sqrt 2);
  * fixed encoding  (E psi)_a = w_a^* psi / sqrt M,  an exact isometry;
  * for a target U with realification O_U, the coupling (A, B) = (Q(X), Q(O_U X)), X ~ N(0, I_{2N}),
    is doubly stochastic; any Birkhoff decomposition into permutations pi_l with weights p_l gives
        sum_l p_l Delta_l^* Delta_l  <=  (4 e^2 / s^2) I_N,   Delta_l = P_{pi_l} E - E U;
  * a uniform coherent label register over R copies turns the mixture into ONE permutation Pi_U with
        ||Pi_U Et - Et U||^2 = || sum_l (k_l/R) Delta_l^* Delta_l ||.

Checks (numpy only, deterministic):
 (1) quantizer: e^2 <= 32 K^{-1/3} for K = 4, 8, 16, 64;
 (2) E^* E = I_N exactly (N = 1, 2; K = 4);
 (3) coupling error, N = 2, K = 4, three random targets, Monte Carlo:
        E (w_B - U w_A)(w_B - U w_A)^*  <=  (4 e^2 / s^2) I_N   (max eigenvalue vs bound),
     and the packet reading  sum_{a,b} Pr(a,b) |(w_a - U^* w_b)^* psi|^2 = E|(w_A - U^* w_B)^* psi|^2;
 (4) the coherent-label identity  ||Pi Et - Et U||_op^2 = || sum_l (k_l/R) Delta_l^* Delta_l ||_op
     on a synthetic Birkhoff mixture (N = 2, K = 4, M = 256, R = 8), and
     the mixture-to-coupling identity  sum_l p_l ||Delta_l psi||^2 = sum_{a,b} D_{ba}/M |(w_a - U^* w_b)^* psi|^2.

Exit 0 iff every check passes.
"""
import sys
import math
import itertools
import numpy as np

rng = np.random.default_rng(20260912)


# ---------------------------------------------------------------- standard normal helpers
def Phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def phi(x):
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def Phi_inv(p):
    lo, hi = -12.0, 12.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if Phi(mid) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def quantizer(K):
    """Equiprobable K-level quantizer of N(0,1): boundaries t_0..t_K and conditional means q."""
    t = [-math.inf] + [Phi_inv(k / K) for k in range(1, K)] + [math.inf]
    q = np.array([K * ((phi(t[k]) if np.isfinite(t[k]) else 0.0)
                       - (phi(t[k + 1]) if np.isfinite(t[k + 1]) else 0.0)) for k in range(K)])
    s2 = float(np.mean(q ** 2))
    return np.array(t), q, s2, 1.0 - s2


def quantize(X, t):
    """Interval index of each entry of X under boundaries t (vectorized)."""
    return np.searchsorted(t[1:-1], X, side="right")


def realify(U):
    return np.block([[U.real, -U.imag], [U.imag, U.real]])


def rand_unitary(n):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Qm, R = np.linalg.qr(z)
    return Qm * (np.diag(R) / np.abs(np.diag(R)))


def w_vectors(idx, q, s2, N):
    """idx: (..., 2N) integer array of addresses -> (..., N) complex w vectors."""
    qq = q[idx]
    return (qq[..., :N] + 1j * qq[..., N:]) / math.sqrt(2.0 * s2)


fails = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


print("=" * 74)
print("verify_koopman_transport.py -- permutation transport, small cases")
print("=" * 74)

# ---------------------------------------------------------------- (1)
ok = True
for K in (4, 8, 16, 64):
    _, q, s2, e2 = quantizer(K)
    ok &= e2 <= 32 * K ** (-1 / 3) and abs(np.mean(q)) < 1e-12
    print(f"      K={K:3d}: s^2 = {s2:.5f}, e^2 = {e2:.5f}, bound 32 K^(-1/3) = {32 * K ** (-1/3):.4f}")
check("(1) quantizer error bound e^2 <= 32 K^(-1/3), zero mean", ok)

# ---------------------------------------------------------------- (2)
ok = True
for N in (1, 2):
    K = 4
    t, q, s2, e2 = quantizer(K)
    grid = np.array(list(itertools.product(range(K), repeat=2 * N)))     # M x 2N
    W = w_vectors(grid, q, s2, N)                                         # M x N
    M = W.shape[0]
    E = np.conj(W) / math.sqrt(M)                                         # rows w_a^* / sqrt M
    ok &= np.allclose(E.conj().T @ E, np.eye(N), atol=1e-12)
check("(2) E^*E = I_N exactly (N = 1, 2; K = 4)", ok)

# ---------------------------------------------------------------- (3)
N, K = 2, 4
t, q, s2, e2 = quantizer(K)
bound = 4 * e2 / s2
ok_op = ok_packet = True
samples = 400_000
for trial in range(3):
    U = rand_unitary(N)
    O = realify(U)
    X = rng.normal(size=(samples, 2 * N))
    A = quantize(X, t)
    Bq = quantize(X @ O.T, t)
    wA, wB = w_vectors(A, q, s2, N), w_vectors(Bq, q, s2, N)
    diff = wB - (U @ wA.T).T                                              # w_B - U w_A
    Sigma = (diff[:, :, None] * np.conj(diff[:, None, :])).mean(axis=0)  # E diff diff^*
    top = np.linalg.eigvalsh((Sigma + Sigma.conj().T) / 2).max()
    ok_op &= top <= bound * (1 + 0.03)                                    # 3% Monte Carlo slack
    # packet reading: E |(w_A - U^* w_B)^* psi|^2 for a random unit psi equals the same
    # quadratic form  psi^* E[(w_A - U^*w_B)(...)^*] psi  (consistency of the two readings)
    psi = rng.normal(size=N) + 1j * rng.normal(size=N)
    psi /= np.linalg.norm(psi)
    d2 = wA - (U.conj().T @ wB.T).T
    lhs = np.mean(np.abs(np.conj(d2) @ psi) ** 2)
    rhs = np.real(psi.conj() @ (U.conj().T @ Sigma @ U) @ psi)
    ok_packet &= abs(lhs - rhs) < 1e-9 + 1e-6 * abs(rhs)
    print(f"      target {trial}: max eig E diff diff^* = {top:.5f}  vs  4e^2/s^2 = {bound:.5f}")
check("(3) coupling error operator bound  E (w_B - U w_A)(.)^* <= (4e^2/s^2) I", ok_op)
check("(3) packet reading consistent with the coupling reading", ok_packet)

# ---------------------------------------------------------------- (4)
N, K, R = 2, 4, 8
t, q, s2, e2 = quantizer(K)
grid = np.array(list(itertools.product(range(K), repeat=2 * N)))
Wv = w_vectors(grid, q, s2, N)
M = Wv.shape[0]
E = np.conj(Wv) / math.sqrt(M)
U = rand_unitary(N)
# synthetic Birkhoff mixture: L permutations with integer multiplicities summing to R
L = 5
perms = [rng.permutation(M) for _ in range(L)]
k = np.array([2, 2, 2, 1, 1])
assert k.sum() == R
# Convention: P_pi maps coordinate a to pi(a), i.e. (P_pi)_{pi(a), a} = 1, so (P_pi E)_{pi(a)} = E_a
# and D_{ba} = sum_l p_l [pi_l(a) = b].
def perm_matrix(perm):
    Ppi = np.zeros((M, M))
    Ppi[perm, np.arange(M)] = 1.0
    return Ppi


D = sum((k[l] / R) * perm_matrix(perms[l]) for l in range(L))
Deltas = [perm_matrix(perms[l]) @ E - E @ U for l in range(L)]
S = sum((k[l] / R) * (Deltas[l].conj().T @ Deltas[l]) for l in range(L))
# coherent label register: Et = |+_R> (x) E, Pi = sum_j |j><j| (x) P_{pi_{l(j)}}
labels = np.repeat(np.arange(L), k)                                        # length R
Et = np.kron(np.ones((R, 1)) / math.sqrt(R), E)
Pi = np.zeros((R * M, R * M))
for j in range(R):
    Ppi = np.zeros((M, M))
    Ppi[perms[labels[j]], np.arange(M)] = 1.0
    Pi[j * M:(j + 1) * M, j * M:(j + 1) * M] = Ppi
lhs = np.linalg.norm(Pi @ Et - Et @ U, 2) ** 2
rhs = np.linalg.norm(S, 2)
check("(4) coherent-label identity ||Pi Et - Et U||^2 = ||sum (k_l/R) Delta_l^* Delta_l||",
      abs(lhs - rhs) < 1e-10, f"lhs {lhs:.6f} rhs {rhs:.6f}")
# mixture-to-coupling identity on the synthetic D
psi = rng.normal(size=N) + 1j * rng.normal(size=N)
psi /= np.linalg.norm(psi)
mix = sum((k[l] / R) * np.linalg.norm(Deltas[l] @ psi) ** 2 for l in range(L))
Upsi = U @ psi
coupling = 0.0
nz = np.argwhere(D > 0)
for b, a in nz:
    coupling += (D[b, a] / M) * abs(np.conj(Wv[a]) @ psi - np.conj(Wv[b]) @ Upsi) ** 2
check("(4) mixture-to-coupling identity sum p_l ||Delta_l psi||^2 = sum D_{ba}/M |...|^2",
      abs(mix - coupling) < 1e-10, f"{mix:.6e} vs {coupling:.6e}")

print("-" * 74)
if fails:
    print("verify_koopman_transport.py: FAIL -> " + ", ".join(fails))
    sys.exit(1)
print("verify_koopman_transport.py: PASS")
sys.exit(0)
