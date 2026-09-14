"""Targeted finite controls for the 2026-09-13 release corrections.

Not a replacement for the analytic proofs or the nine-script numerical suite.
Run under the campaign's memory/time wrapper.
"""
from itertools import product
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20260913)
worst = 0.0
checks = 0


def near(a, b):
    global worst, checks
    err = float(np.linalg.norm(a - b))
    worst = max(worst, err)
    checks += 1
    assert err < 1e-11, err


def unitary(d):
    a = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return np.linalg.qr(a)[0]


# Address, answer, spectator all retained: (a,b,r), each of size two.
H = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
S = np.kron(np.kron(np.eye(2), H), np.eye(2))
Jin = np.eye(8)[:, [0, 4]]
Jout = np.eye(8)[:, [3, 7]]
A0, A1, A2 = (unitary(8) for _ in range(3))
for f in product((0, 1), repeat=2):
    oracle = np.zeros((8, 8))
    diag = np.empty(8)
    for a, b, r in product((0, 1), repeat=3):
        j = 4*a + 2*b + r
        oracle[4*a + 2*(b ^ f[a]) + r, j] = 1
        diag[j] = (-1)**(b*f[a])
    D = np.diag(diag)
    near(S @ oracle @ S, D)
    full = A2 @ oracle @ A1 @ oracle @ A0
    B, W, C = Jout.conj().T @ A2 @ S, S @ A1 @ S, S @ A0 @ Jin
    near(Jout.conj().T @ full @ Jin, B @ D @ W @ D @ C)

# An explicit almost-clean isometry; full-output leakage is not channel error.
U = unitary(2)
Jperp = np.eye(8)[:, [2, 6]]
theta = 0.13
VJ = np.cos(theta)*Jout @ U + np.sin(theta)*Jperp @ U
near(VJ.conj().T @ VJ, np.eye(2))
near(np.linalg.norm(VJ - Jout @ U, 2), 2*np.sin(theta/2))
projected = np.linalg.norm(Jout.conj().T @ VJ - U, 2)
assert projected <= np.linalg.norm(VJ - Jout @ U, 2) + 1e-12
checks += 1

# Real-linear norm and complex covariance pushforward have different constants.
for n, s in ((1, .5), (2, .8), (4, 1.0)):
    L = np.hstack((np.eye(n), 1j*np.eye(n))) / (s*np.sqrt(2))
    near(L @ L.conj().T, np.eye(n)/s**2)
    for _ in range(8):
        x = rng.normal(size=2*n)
        near(np.linalg.norm(L @ x)**2, np.linalg.norm(x)**2/(2*s**2))
        a = rng.normal(size=(2*n, 2*n))
        Sigma = a @ a.T
        Sigma /= np.linalg.norm(Sigma, 2)
        assert np.linalg.eigvalsh(np.eye(n)/s**2 - L @ Sigma @ L.conj().T)[0] > -1e-11
        checks += 1

for path in (ROOT / 'paper/paper.tex', ROOT / 'paper/paper.md'):
    text = path.read_text(encoding='utf-8')
    assert 'From oracle circuits to Definition' in text
    assert 'negatively in the clean model' in text
    assert 'full-output approximate-clean guarantee' in text
    assert 'whereas covering all' in text
    assert 'deterministic spine' not in text
print(f'PASS: {checks} finite circuit/covariance controls; maximum residual {worst:.3g}')
print('PASS: corrected TeX/Markdown scope and bridge markers')
