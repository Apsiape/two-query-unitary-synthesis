"""Section 8 (garbage model) -- the junk-blindness lemma and the conjugation transfer.

Model: a t-query circuit V(g) = A_t (D_g (x) I_K) ... (D_g (x) I_K) A_0 on C^N (x) C^K, with A_i
fixed unitaries, restricted to the clean input:  S_g = V(g) iota,  iota psi = psi (x) |0>.
For a contraction X on the target register set  T_g = S_g^* (X (x) I_K) S_g.

Checks (numpy only, deterministic):
 (1) reversal identity: T_g equals the explicit 2t-insertion product
        iota^* A_0^* (D (x) I) A_1^* ... (D (x) I) A_t^* (X (x) I) A_t (D (x) I) ... A_0 iota
     for t = 2 (four insertions of the SAME sign word), and every slot is a contraction;
 (2) junk cancellation: if  S_g psi = (U_g psi) (x) |j_g>  (keep form) then  T_g = U_g^* X U_g,
     with no dependence on j_g;
 (3) balanced reflections: for R = R^*, R^2 = I, Tr R = 0 and unitary U, the image U^* R U is again
     a balanced reflection with ||U^* R U||_F = sqrt N  (the conjugation orbit carries the code);
 (4) junk-blindness lemma witness: a linear functional L with L(S) != 0 is not invariant under
     S -> e^{i theta} S, so no nonzero linear functional of the packet is junk-blind.

Exit 0 iff every check passes.
"""
import sys
import numpy as np

rng = np.random.default_rng(20260912)


def rand_unitary(n):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Qm, R = np.linalg.qr(z)
    return Qm * (np.diag(R) / np.abs(np.diag(R)))


def contraction(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return A / (np.linalg.norm(A, 2) * 1.2)


fails = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


print("=" * 74)
print("verify_conjugation_transfer.py -- junk-blindness and the conjugation transfer")
print("=" * 74)

# ---------------------------------------------------------------- (1)
N, K = 3, 2
dim = N * K
A0, A1, A2 = rand_unitary(dim), rand_unitary(dim), rand_unitary(dim)
iota = np.kron(np.eye(N), np.array([[1.0], [0.0]]))        # psi -> psi (x) |0>, dim x N
X = contraction(N)
ok_rev = ok_contr = True
for _ in range(20):
    g = rng.choice([1.0, -1.0], N)
    DI = np.kron(np.diag(g), np.eye(K))
    V = A2 @ DI @ A1 @ DI @ A0
    S = V @ iota
    T = S.conj().T @ np.kron(X, np.eye(K)) @ S
    slots = [iota.conj().T @ A0.conj().T, A1.conj().T, A2.conj().T @ np.kron(X, np.eye(K)) @ A2, A1, A0 @ iota]
    explicit = slots[0] @ DI @ slots[1] @ DI @ slots[2] @ DI @ slots[3] @ DI @ slots[4]
    ok_rev &= np.allclose(T, explicit, atol=1e-12)
    ok_contr &= all(np.linalg.norm(s, 2) <= 1 + 1e-12 for s in slots)
check("(1) reversal identity: T_g is an exact four-insertion product of the same sign word", ok_rev)
check("(1) every slot of the reversed product is a contraction", ok_contr)

# ---------------------------------------------------------------- (2)
ok = True
for _ in range(20):
    U = rand_unitary(N)
    j = rng.normal(size=K) + 1j * rng.normal(size=K)
    j /= np.linalg.norm(j)
    S = np.kron(U, j.reshape(-1, 1))                        # psi -> (U psi) (x) |j>
    T = S.conj().T @ np.kron(X, np.eye(K)) @ S
    ok &= np.allclose(T, U.conj().T @ X @ U, atol=1e-12)
check("(2) keep form: S^*(X (x) I)S = U^* X U, junk cancels identically", ok)

# ---------------------------------------------------------------- (3)
N = 4
R = np.diag([1.0, 1.0, -1.0, -1.0])
ok = True
for _ in range(20):
    U = rand_unitary(N)
    Y = U.conj().T @ R @ U
    ok &= np.allclose(Y, Y.conj().T) and np.allclose(Y @ Y, np.eye(N)) and abs(np.trace(Y)) < 1e-10
    ok &= abs(np.linalg.norm(Y, "fro") - np.sqrt(N)) < 1e-10
check("(3) U^* R U is a balanced reflection with Frobenius norm sqrt N", ok)

# ---------------------------------------------------------------- (4)
S = np.kron(rand_unitary(3), np.array([[1.0], [0.0]]))
Lmat = rng.normal(size=S.shape) + 1j * rng.normal(size=S.shape)
L = lambda Z: np.sum(np.conj(Lmat) * Z)                     # a generic linear functional
val = L(S)
invariant = all(abs(L(np.exp(1j * th) * S) - val) < 1e-12 for th in (np.pi / 3, np.pi, 1.0))
check("(4) a nonzero linear functional is not junk-blind (phase rotation changes it)",
      abs(val) > 1e-8 and not invariant)

print("-" * 74)
if fails:
    print("verify_conjugation_transfer.py: FAIL -> " + ", ".join(fails))
    sys.exit(1)
print("verify_conjugation_transfer.py: PASS")
sys.exit(0)
