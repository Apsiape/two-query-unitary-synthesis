"""
PIN C_rho against the Dong-Lombardi-Ma permutation family.

The chain (paper, Corollary 5.1):
    log Pack_{rho sqrt N} <= 2 C_S^2 * w_2^2 / (rho^2 N),   w_2 <= C1 * min(N,Q) * sqrt(log 2Q)
    =>  log Pack <= (kappa / rho^2) * N * log(2Q),              kappa := 2 C_S^2 C1^2
where C_S is the Sudakov constant (eps sqrt(log Pack) <= C_S E sup) and C1 the width constant
(2 for Theorem 4.1; 8 for the audited real-selector form).  Proposition 5.3 of the paper.

DLM give a CLEAN TWO-QUERY algorithm for all d! permutation unitaries.
That family is inside our class, so it must fit under the ceiling.
Requirement:   K  >=  rho^2 * logPack(rho) / (d * log(2Q))   for every rho, d.

If K < that maximum, the theorem is FALSE -- refuted by a published construction.
"""
import math
import sys
from functools import lru_cache

fails = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


results = {}

@lru_cache(maxsize=65537)
def log_factorial(n):
    return math.lgamma(n + 1)

_log_derangements = [0.0, -math.inf]  # !0 = 1 and !1 = 0


def log_derangement(n):
    """Evaluate !n=(n-1)(!(n-1)+!(n-2)) in log space, not n!/e.

    The recurrence is exact; its floating-point evaluation is not an interval
    certificate. In particular, the impossible one-point derangement has zero
    weight rather than the spurious 1/e used by the former approximation.
    """
    while len(_log_derangements) <= n:
        j = len(_log_derangements)
        a, b = _log_derangements[-1], _log_derangements[-2]
        hi, lo = max(a, b), min(a, b)
        _log_derangements.append(math.log(j - 1) + hi + math.log1p(math.exp(lo - hi)))
    return _log_derangements[n]

def log_ball(d, r):
    """log |{pi : hamming_dist(pi, id) <= r}|  = log sum_{k>=d-r} C(d,k) * !(d-k)
       (agreements >= d-r).  Computed in log-space."""
    terms = []
    for k in range(d - r, d + 1):          # k = number of agreements
        if k < 0:
            continue
        terms.append(log_factorial(d) - log_factorial(k) - log_factorial(d - k)
                     + log_derangement(d - k))
    if not terms:
        return -math.inf
    m = max(terms)
    return m + math.log(sum(math.exp(t - m) for t in terms))

def log_pack_permutations(d, alpha):
    """Numerical evaluation of the exact GV lower-bound expression for a code with minimum
       Hamming distance D = alpha*d, via greedy deletion:
           M >= d! / |ball of radius D-1|      (Gilbert-Varshamov)
       Returns natural log; floating-point values are not certified enclosures."""
    D = max(1, int(math.ceil(alpha * d)))
    lb = log_ball(d, D - 1)
    return log_factorial(d) - lb

print(__doc__.strip())
print()
check("!0=1, !1=0, !2=1, !3=2, !4=9",
      all(abs(math.exp(log_derangement(n)) - value) < 1e-12
          for n, value in enumerate([1, 0, 1, 2, 9])))
check("exact small deletion ball: d=4, radius=2 has 7 permutations",
      abs(math.exp(log_ball(4, 2)) - 7) < 1e-12)

CONVENTIONS = {
    "Q = 2d^2  (address count, per the DLM packet reconstruction)": lambda d: 2*d*d,
    "Q = d     (the alternative convention)":                        lambda d: d,
}

for label, Qof in CONVENTIONS.items():
    print("=" * 78)
    print(label)
    print("=" * 78)
    print(f"  {'d':>8} {'best alpha':>11} {'rho^2':>8} {'logPack':>12} "
          f"{'d*log(2Q)':>12} {'required kappa':>12}")
    worst_K = 0.0
    for n in [4, 6, 8, 10, 12, 14, 16]:
        d = 2**n
        Q = Qof(d)
        denom = d * math.log(2*Q)
        best = (0.0, None, None, None)
        for i in range(1, 200):
            alpha = i/200.0
            if alpha >= 1.0:
                continue
            lp = log_pack_permutations(d, alpha)
            if lp <= 0:
                continue
            rho2 = 2.0*alpha                     # ||P-P'||_F^2 = 2*D ; normalized: 2*alpha
            K = rho2 * lp / denom
            if K > best[0]:
                best = (K, alpha, rho2, lp)
        K, alpha, rho2, lp = best
        worst_K = max(worst_K, K)
        results.setdefault(label, {})[d] = K
        print(f"  {d:>8} {alpha:>11.3f} {rho2:>8.3f} {lp:>12.1f} "
              f"{denom:>12.1f} {K:>12.4f}")
    print(f"\n  -> REQUIRED  kappa = C * C1^2  >=  {worst_K:.4f}   (largest over d tested)")
    # Asymptotics, done correctly.
    #   GV:  log M ~ (1-a) d (ln((1-a)d) - 1)   with a = alpha = D/d
    #   K   = rho^2 * logM / (d ln 2Q),  rho^2 = 2a
    #       ~ 2a(1-a) ln d / ln(2Q)
    #   Q = 2d^2 : ln 2Q ~ 2 ln d  =>  K -> a(1-a),      max 1/4 at a = 1/2
    #   Q = d    : ln 2Q ~   ln d  =>  K -> 2a(1-a),     max 1/2 at a = 1/2
    #
    # NOTE: one must NOT assume the packing tends to all of S_d as rho^2 -> 2.
    # It does not.  For minimum distance alpha*d each permutation has about
    # d!/((1-alpha)d)! close neighbours, so greedy deletion leaves only about
    # ((1-alpha)d)! codewords -- consistent with the Singleton-type bound
    # M <= d!/(D-1)!.  Assuming otherwise inflates the requirement ~4x.
    lim = 0.25 if "2d^2" in label else 0.5
    print(f"  -> asymptotic requirement (GV, max at alpha = 1/2): K >= {lim:.4f}")
    print(f"     finite-d column above is converging up to this value\n")

print("=" * 78)
print("ASSERTIONS (paper, Proposition 5.3)")
print("=" * 78)
main = [k for k in results if "2d^2" in k][0]
expected = {16: 0.1256, 256: 0.1591, 4096: 0.1849, 65536: 0.1999}
for d, K_exp in expected.items():
    check(f"required kappa at d = {d} matches the paper's table ({K_exp:.4f})",
          abs(results[main][d] - K_exp) < 2e-3, f"computed {results[main][d]:.4f}")
seq = [results[main][d] for d in sorted(results[main])]
check("required kappa increases with d and stays below the asymptotic 1/4",
      all(a < b for a, b in zip(seq, seq[1:])) and seq[-1] < 0.25)


# No finite list of quadratures establishes a universal Sudakov constant.
# The paper uses the standard theorem with an unspecified absolute constant.
print("INFO: no numerical upper bound on the universal Sudakov constant is certified")
print("-" * 78)
if fails:
    print("verify_crho_requirement.py: FAIL -> " + ", ".join(fails))
    sys.exit(1)
print("verify_crho_requirement.py: PASS")
sys.exit(0)
