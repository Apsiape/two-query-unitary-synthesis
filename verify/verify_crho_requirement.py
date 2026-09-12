"""
PIN C_rho against the Dong-Lombardi-Ma permutation family.

The chain (paper, Corollary 5.1):
    log Pack_{rho sqrt N} <= 2 C_S^2 * w_2^2 / (rho^2 N),   w_2 <= C1 * min(N,Q) * sqrt(log 2Q)
    =>  log Pack <= (K / rho^2) * N * log(2Q),              K := 2 C_S^2 C1^2
where C_S is the Sudakov constant (eps sqrt(log Pack) <= C_S E sup) and C1 the width constant
(2 for Theorem 4.1; 8 for the audited real-selector form).  Proposition 5.3 of the paper.

DLM give a CLEAN TWO-QUERY algorithm for all d! permutation unitaries.
That family is inside our class, so it must fit under the ceiling.
Requirement:   K  >=  rho^2 * logPack(rho) / (d * log(2Q))   for every rho, d.

If K < that maximum, the theorem is FALSE -- refuted by a published construction.
"""
import math
import sys

fails = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


results = {}

def log_factorial(n):
    return math.lgamma(n + 1)

def log_derangement(n):
    # !n = round(n!/e); log form
    return log_factorial(n) - 1.0 if n >= 1 else 0.0

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
    """Rigorous lower bound on the size of a permutation code with minimum
       Hamming distance D = alpha*d, via greedy deletion:
           M >= d! / |ball of radius D-1|      (Gilbert-Varshamov)
       Returns natural log."""
    D = max(1, int(math.ceil(alpha * d)))
    lb = log_ball(d, D - 1)
    return log_factorial(d) - lb

print(__doc__.strip())
print()

CONVENTIONS = {
    "Q = 2d^2  (address count, per the DLM packet reconstruction)": lambda d: 2*d*d,
    "Q = d     (the alternative convention)":                        lambda d: d,
}

for label, Qof in CONVENTIONS.items():
    print("=" * 78)
    print(label)
    print("=" * 78)
    print(f"  {'d':>8} {'best alpha':>11} {'rho^2':>8} {'logPack':>12} "
          f"{'d*log(2Q)':>12} {'required K':>12}")
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
    print(f"\n  -> REQUIRED  K = C * C1^2  >=  {worst_K:.4f}   (largest over d tested)")
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
    check(f"required K at d = {d} matches the paper's table ({K_exp:.4f})",
          abs(results[main][d] - K_exp) < 2e-3, f"computed {results[main][d]:.4f}")
seq = [results[main][d] for d in sorted(results[main])]
check("required K increases with d and stays below the asymptotic 1/4",
      all(a < b for a, b in zip(seq, seq[1:])) and seq[-1] < 0.25)

# The Sudakov constant, bounded through the Sudakov-Fernique comparison with independent
# N(0, eps^2/2) variables:  E sup X  >=  (eps/sqrt 2) * E max_m g_i,  and
# E max of m standard normals  =  int x * m * phi(x) * Phi(x)^(m-1) dx  >=  0.677 sqrt(log m)
# for every m >= 2 (the minimum of the ratio is at m = 2).  Hence eps sqrt(log m) <= C_S E sup
# with  C_S = sqrt(2) / 0.677 = 2.09.


def sudakov_fernique_ratio(m, n=400000, x_max=12.0):
    xs = [-x_max + 2 * x_max * (i + 0.5) / n for i in range(n)]
    h = 2 * x_max / n
    tot = 0.0
    for x in xs:
        ph = math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi)
        Ph = 0.5 * (1.0 + math.erf(x / math.sqrt(2)))
        tot += x * m * ph * Ph ** (m - 1) * h
    return tot / math.sqrt(math.log(m))


ratios = {m: sudakov_fernique_ratio(m, n=60000) for m in (2, 3, 4, 8, 16, 64, 256)}
min_ratio = min(ratios.values())
print("  E max_m / sqrt(log m):", ", ".join(f"m={m}: {r:.4f}" for m, r in ratios.items()))
check("E max of m standard normals >= 0.677 sqrt(log m), minimum at m = 2 (= 1/sqrt(pi log 2))",
      min_ratio >= 0.677 and abs(ratios[2] - 1 / math.sqrt(math.pi * math.log(2))) < 1e-3,
      f"min ratio {min_ratio:.4f}")
C_S = math.sqrt(2) / min_ratio
K_thm41 = 2 * C_S ** 2 * 2 ** 2
K_audited = 2 * C_S ** 2 * 8 ** 2
print(f"  C_S <= {C_S:.3f};  K = 2 C_S^2 C1^2 = {K_thm41:.1f} (C1 = 2, Theorem 4.1), "
      f"{K_audited:.0f} (C1 = 8, audited form)")
check("C_S <= 2.09 and K <= 35 for the constant-2 route", C_S <= 2.09 and K_thm41 <= 35.0)
check("both instantiations exceed the required K (no contradiction with the DLM family)",
      K_thm41 >= worst_K and K_audited >= worst_K,
      f"required {max(results[main].values()):.4f} (finite d), 0.25 (asymptotic)")

print("-" * 78)
if fails:
    print("verify_crho_requirement.py: FAIL -> " + ", ".join(fails))
    sys.exit(1)
print("verify_crho_requirement.py: PASS")
sys.exit(0)
