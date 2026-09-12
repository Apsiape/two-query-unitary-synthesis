# Changelog

## 1.0.1 — 2026-09-11

Fixes from the five-lane refutation-first audit of 1.0.0 (details in `STATUS.md`, "Audit
round of 2026-09-11"). No theorem was found false.

- Proposition 2.2: third term for a contractive $`C`$; verifier now uses contractive $`C`$
  and $`W`$.
- Section 3 rewritten: new Proposition 3.1 (at depth two the final layer cannot change the
  target count); the decisive pair is stated in the three-insertion, one-dimensional-target
  model in which it is verified (Proposition 3.2).
- Corollary 5.2: quantitative packing lemma for $`U(N)`$ with a proof; Proposition 5.3:
  derived Sudakov constant $`C_S\le2.09`$ replacing an unsupported interval.
- Corollary 6.2: correct statement at exponential $`Q`$ (the logarithm is necessary, by
  Theorem 6.4). Theorem 6.3: two different address decompositions allowed. Theorem 6.4:
  repeated-sign-word form, operator-norm error, all symbols defined.
- Section 7: Theorem 7.2 restated in address form with a rate function; Stinespring
  continuity cited as [KSW08] plus [vE23] with the enlargement argument; palindromic form
  noted; Corollary 7.3 scoped to constant error.
- Literature corrections: parallel-query bound attributed to [LMW24]; [DLM26] indexing;
  [Ros26] ancilla count; [AK07] footnote 13; model paragraph of Section 1.2.
- Verifiers: assertions added to the decisive-pair and constant scripts; torus supremum by
  alternating maximisation; end-to-end permutation transport at $`N=1`$; real-linear
  junk-blindness; two-decomposition compression; linter rejects tab bytes.
- Wording: "exactly", "nothing was known", "for every oracle" and similar overclaims
  removed from the abstract, introduction and README.

## 1.0.0 — 2026-09-11

First standalone release.

- Paper: the two-query Gaussian mean width theorem extended from real sign words to
  independent unimodular selectors (Theorem 4.1, constant 2), with the audited real form
  retained (Corollary 4.2, constant 8); corrected unitarity criterion (Proposition 2.2);
  decisive pair (Section 3); capacity and qubit lower bound with the metric erratum
  stated in place (Corollaries 5.1–5.2); the constant requirement (Proposition 5.3);
  tightness at $`Q=\Theta(N^2)`$ (Proposition 6.1, Corollary 6.2); address-length lower
  bound independent of workspace (Theorem 6.3); matching permutation-transport upper bound
  (Theorem 6.4); thresholds (Corollary 6.5); junk-blindness and conjugation transfer for the
  garbage model (Section 7); scope (Section 8).
- Notes: the full permutation-transport construction with all estimates and its
  fixed-encoding converse; the full garbage-to-clean transfer.
- Verification: eight deterministic numpy scripts under `verify/`, run by `run_all.py`.
- Lean: five kernel-checked theorems in `lean/TwoQuery/Q2Width.lean`, built and axiom-checked
  on this date; recorded output in `lean/RECORDED-AXIOMS.txt`.
- Ledger: `STATUS.md` labels every claim as audited, proved here, or not claimed.
