"""Core helpers for the Lookup-or-Amplify dichotomy test.

Architecture: (B, W1, W2, C).  W1,W2 in U(Q); C: Q x N isometry (C*C=I_N);
B: N x Q coisometry (B B* = I_N).  P_B = B* B rank-N projection.
Psi(g) = D_g W2 D_g W1 D_g C .   Lawful iff (I-P_B)Psi(g) = 0 ; U_g = B Psi(g).

HARD CAPS honoured everywhere: N<=6, Q<=16, matrices <=24x24, <=2^16 words.
"""
import numpy as np

# ---------------------------------------------------------------- words
def words(Q, cap=1 << 16):
    assert (1 << Q) <= cap, f"cube 2^{Q} exceeds cap"
    idx = np.arange(1 << Q)
    bits = ((idx[:, None] >> np.arange(Q)[None, :]) & 1).astype(np.int8)
    return (1 - 2 * bits).astype(np.float64)


def rand_unitary(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r).copy()
    return q * (d / np.abs(d))


# ---------------------------------------------------------------- frames
def frames(C, W1, W2, G):
    """J[0]=C (broadcast), J[1]=D_g C, J[2]=D_g W1 J[1], J[3]=D_g W2 J[2]."""
    M, Q = G.shape
    N = C.shape[1]
    J1 = G[:, :, None] * C[None, :, :]
    J2 = G[:, :, None] * np.einsum('ab,mbn->man', W1, J1)
    J3 = G[:, :, None] * np.einsum('ab,mbn->man', W2, J2)
    return [np.broadcast_to(C, (M, Q, N)), J1, J2, J3]


class Arch:
    def __init__(self, name, C, W1, W2, B, G=None, tol=1e-9):
        self.name = name
        self.C, self.W1, self.W2, self.B = C, W1, W2, B
        self.Q, self.N = C.shape
        assert W1.shape == (self.Q, self.Q) and W2.shape == (self.Q, self.Q)
        assert B.shape == (self.N, self.Q)
        self.G = words(self.Q) if G is None else G
        J = frames(C, W1, W2, self.G)
        PB = B.conj().T @ B
        resid = np.linalg.norm(J[3] - np.einsum('ab,mbn->man', PB, J[3]), axis=(1, 2))
        self.law = np.where(resid < tol)[0]
        self.J = [j[self.law].copy() for j in J]        # keep only lawful (memory)
        self.Gl = self.G[self.law]
        self.U = np.einsum('aq,mqn->man', B, self.J[3])
        del J, resid

    # ------- sanity ---------------------------------------------------
    def checks(self):
        out = {}
        out['C isometry']  = np.max(np.abs(self.C.conj().T @ self.C - np.eye(self.N)))
        out['B coisometry']= np.max(np.abs(self.B @ self.B.conj().T - np.eye(self.N)))
        out['W1 unitary']  = np.max(np.abs(self.W1.conj().T @ self.W1 - np.eye(self.Q)))
        out['W2 unitary']  = np.max(np.abs(self.W2.conj().T @ self.W2 - np.eye(self.Q)))
        out['U unitary on L'] = np.max(np.abs(
            np.einsum('mai,maj->mij', self.U.conj(), self.U) - np.eye(self.N)))
        # T_x* T_x = I - Psi*(I-P_B)Psi  (the stated exact fact), on a sample of ALL words
        PB = self.B.conj().T @ self.B
        Gs = self.G[np.linspace(0, len(self.G) - 1, min(64, len(self.G))).astype(int)]
        Js = frames(self.C, self.W1, self.W2, Gs)[3]
        T = np.einsum('aq,mqn->man', self.B, Js)
        lhs = np.einsum('mai,maj->mij', T.conj(), T)
        rhs = np.eye(self.N) - np.einsum('mai,mab,mbj->mij', Js.conj(),
                                         np.broadcast_to(np.eye(self.Q) - PB, (len(Gs), self.Q, self.Q)), Js)
        out['T*T identity'] = np.max(np.abs(lhs - rhs))
        # J_s(x) is an isometry at every stage?
        out['J_s isometry'] = max(np.max(np.abs(
            np.einsum('mai,maj->mij', self.J[s].conj(), self.J[s]) - np.eye(self.N)))
            for s in range(4))
        # J_3 = B* U  on lawful
        out['J3 = B* U'] = np.max(np.abs(self.J[3] - np.einsum('qa,man->mqn',
                                                               self.B.conj().T, self.U)))
        return out

    def n_targets(self, dec=6):
        keys = set()
        for uu in self.U:
            keys.add(tuple(np.round(uu.ravel(), dec)))
        return len(keys)


# ---------------------------------------------------------------- families
def hadamard(n):
    H = np.array([[1.0]])
    while H.shape[0] < n:
        H = np.block([[H, H], [H, -H]])
    return H


def lookup_block(q, lam):
    """Exact t=3 lookup compiler on one block of size q=2^m (N=1).
    Returns c, W1, W2, b, L(codewords).  U_g = lam_g on codeword g."""
    m = int(np.log2(q)); assert 1 << m == q
    polys = {2: [1, 1], 3: [1, 1, 0], 4: [1, 1, 0, 0], 5: [1, 0, 1, 0, 0]}
    A = np.zeros((m, m), dtype=np.uint8)
    A[1:, :-1] = np.eye(m - 1, dtype=np.uint8)
    A[:, -1] = np.array(polys[m], dtype=np.uint8)
    xs = ((np.arange(q)[:, None] >> np.arange(m)[None, :]) & 1).astype(np.uint8)
    chi = lambda v: (-1.0) ** (xs @ v % 2)
    f = np.ones(q, complex) / np.sqrt(q)
    c = f.copy(); b = f.copy()
    L = [chi(xs[i]) for i in range(q)]
    Y = [chi((A @ xs[i]) % 2) * f for i in range(q)]
    W1 = sum(np.outer(Y[i], (L[i] * c).conj()) for i in range(q))
    X = [L[i] * Y[i] for i in range(q)]
    W2 = sum(lam[i] * np.outer(L[i] * b, X[i].conj()) for i in range(q))
    return c, W1, W2, b, np.array(L)


def balanced_promise_block(q):
    """u = q^{-1/2}(1..1); C=e_0; W1 real orth with W1 e_0 = u; W2 = 2uu*-I; B=u*."""
    u1 = np.ones(q) / np.sqrt(q)
    e0 = np.zeros(q); e0[0] = 1.0
    w = u1 - e0
    W1b = np.eye(q) - 2 * np.outer(w, w) / (w @ w)
    if np.linalg.norm(W1b @ e0 - u1) > 1e-9:
        W1b = -W1b
    assert np.linalg.norm(W1b @ e0 - u1) < 1e-9
    W2b = 2 * np.outer(u1, u1) - np.eye(q)
    return e0.astype(complex), W1b.astype(complex), W2b.astype(complex), u1.astype(complex)


def blocks_to_arch(name, blocks):
    """blocks = list of (c,W1,W2,b) each of size q_i -> block-diagonal arch with N=len(blocks)."""
    qs = [len(b[0]) for b in blocks]
    Q = sum(qs); N = len(blocks)
    C = np.zeros((Q, N), complex); B = np.zeros((N, Q), complex)
    W1 = np.zeros((Q, Q), complex); W2 = np.zeros((Q, Q), complex)
    o = 0
    for i, (c, w1, w2, b) in enumerate(blocks):
        q = qs[i]; s = slice(o, o + q)
        C[s, i] = c; B[i, s] = b.conj(); W1[s, s] = w1; W2[s, s] = w2
        o += q
    return Arch(name, C, W1, W2, B)


def diagonal_corner(N, q):
    Q = N * q
    C = np.zeros((Q, N), complex); B = np.zeros((N, Q), complex)
    for i in range(N):
        C[i * q, i] = 1.0; B[i, i * q] = 1.0
    return Arch(f'diagonal corner N={N} q={q}', C, np.eye(Q, dtype=complex),
                np.eye(Q, dtype=complex), B)
