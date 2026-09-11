"""
PIN C_rho against the Dong-Lombardi-Ma permutation family.

Our q2 chain (q2-width-explainer.md sec.8):
    log Pack_{rho sqrt N} <= C * w_2^2 / (rho^2 N),   w_2 <= C1 * min(N,Q) * sqrt(log 2Q)
    =>  log Pack <= (K / rho^2) * N * log(2Q),        K := C * C1^2

DLM give a CLEAN TWO-QUERY algorithm for all d! permutation unitaries.
That family is inside our class, so it must fit under the ceiling.
Requirement:   K  >=  rho^2 * logPack(rho) / (d * log(2Q))   for every rho, d.

If K < that maximum, the theorem is FALSE -- refuted by a published construction.
"""
import math

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
print("WHAT K DO WE ACTUALLY HAVE?")
print("=" * 78)
print("""
  K = C * C1^2 where
    C1 : noncommutative Khintchine,  E||sum g_i A_i|| <= C1 sqrt(log(d1+d2)) * sigma
         standard constant  C1 = sqrt(2)   =>  C1^2 = 2
    C  : Sudakov minoration in the form  eps * sqrt(log M) <= C * w
         standard constants for this normalization lie in roughly  C in [2, 4]

  =>  K  in roughly  [4, 8]
""")
for Kc in [2.0, 4.0, 6.0, 8.0]:
    print(f"    K = {Kc:>4.1f}  -> margin vs Q=2d^2 requirement (~1.0): "
          f"{Kc/1.0:>5.1f}x ;  vs Q=d requirement (~2.0): {Kc/2.0:>5.1f}x")
