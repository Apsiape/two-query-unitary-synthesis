# Changelog

## 1.0.0 — 2026-09-11

First standalone release.

- Paper: the two-query Gaussian mean width theorem extended from real sign words to
  independent unimodular selectors (Theorem 4.1, constant 2), with the audited real form
  retained (Corollary 4.2, constant 8); corrected unitarity criterion (Proposition 2.2);
  decisive pair (Proposition 3.1); capacity and qubit lower bound with the metric erratum
  stated in place (Corollaries 5.1–5.2); the constant requirement (Proposition 5.3);
  tightness at polynomial workspace (Proposition 6.1, Corollary 6.2); address-length lower
  bound independent of workspace (Theorem 6.3); matching permutation-transport upper bound
  (Theorem 6.4); thresholds (Corollary 6.5); junk-blindness and conjugation transfer for the
  garbage model (Section 7); scope (Section 8).
- Notes: the full permutation-transport construction with all estimates and its
  fixed-encoding converse; the full garbage-to-clean transfer.
- Verification: eight deterministic numpy scripts under `verify/`, run by `run_all.py`.
- Lean: five kernel-checked theorems in `lean/TwoQuery/Q2Width.lean`, built and axiom-checked
  on this date; recorded output in `lean/RECORDED-AXIOMS.txt`.
- Ledger: `STATUS.md` labels every claim as audited, proved here, or not claimed.
