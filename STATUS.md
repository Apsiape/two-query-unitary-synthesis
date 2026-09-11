# STATUS — claim-by-claim ledger

*Updated 2026-09-11 for version 1.0.0. This file is the authority on what this repository
claims and how each claim is supported. When it disagrees with a sentence elsewhere, this
file wins.*

## Labels

| label | meaning |
|---|---|
| **AUDITED** | proved in the paper, machine-checked here, and independently re-derived by an audit external to the author on the date given |
| **PROVED HERE** | proved in the paper and machine-checked by a script in `verify/`; no independent re-derivation yet |
| **CERTIFIED (Lean)** | the deterministic part is kernel-checked in `lean/`, with the analytic inputs as explicit hypotheses; see `lean/README.md` |
| **NOT CLAIMED** | discussed for context only; no theorem asserted |

Every PROVED HERE item is a candidate for independent audit. The author regards the
verifier scripts as necessary, not sufficient: they catch algebraic and numerical
mistakes at small size and cannot certify an asymptotic argument.

## Ledger

| # | claim | label | machine check |
|---|---|---|---|
| Prop 2.2 | unitarity criterion with the corrected middle term | PROVED HERE | `verify_proof_steps.py`, `verify_complex_phase.py` (e) |
| Prop 3.1 | decisive pair: identical pre-final data, 2 versus 16 targets | PROVED HERE (exhaustive at $`Q=8`$, $`N=2`$) | `verify_decisive_pair.py` |
| Thm 4.1 | two-query width with independent unimodular selectors, constant 2 | PROVED HERE | `verify_complex_phase.py` (a)–(d), (f) |
| Cor 4.2 | real sign words, constant 2 (as a corollary of Thm 4.1) | PROVED HERE | same |
| Cor 4.2, audited form | real sign words, symmetric weighting, constant 8 | **AUDITED 2026-08-16**; the asymmetric $`\sqrt{LR}`$ route was verified in the same audit for real selectors | `verify_proof_steps.py`; **CERTIFIED (Lean)** for the deterministic spine |
| Cor 5.1 | capacity in the normalized Frobenius metric only | AUDITED 2026-08-16 (as a corollary of the audited width bound) | metric erratum recorded in the paper |
| Cor 5.2 | universal clean two-query synthesis needs $`\Omega(2^n)`$ qubits | AUDITED 2026-08-16 (as above) | — |
| Prop 5.3 | the permutation family forces $`K\ge1/4`$ | PROVED HERE | `verify_crho_requirement.py` |
| Prop 6.1 | the Dong–Lombardi–Ma permutation family lies in the class and fills a constant fraction of the ceiling | PROVED HERE | `verify_permutation_family.py` |
| Cor 6.2 | both factors tight at polynomial $`Q`$ | PROVED HERE | — (Sudakov backwards on Prop 6.1) |
| Thm 6.3 | address-length lower bound $`\Omega(2^n)`$ bits, any workspace | PROVED HERE (new in this version) | `verify_address_compression.py` |
| Thm 6.4 | permutation transport: $`O(2^n\log(1/\varepsilon))`$ qubits and address bits, two Boolean queries, clean | PROVED HERE (new in this version; full construction in `notes/permutation-transport.md`) | `verify_koopman_transport.py` (ingredients at small $`N`$) |
| Cor 6.5 | both thresholds $`\Theta(2^n)`$ up to $`\log(1/\varepsilon)`$ | PROVED HERE | follows from Cor 5.2, Thm 6.3, Thm 6.4 |
| Note §7.1 | monomial-action converse $`(N-1)\log_2(1/\varepsilon)`$ for fixed encodings | PROVED HERE | — |
| Lemma 7.1 | junk-blind linear functionals vanish | PROVED HERE | `verify_conjugation_transfer.py` |
| Thm 7.2 | conjugation transfer: clean-$`2t`$ capacity ⇒ garbage-$`t`$ impossibility | PROVED HERE | `verify_conjugation_transfer.py` (reversal identity, junk cancellation, orbit) |
| Cor 7.3 | one-query garbage-model bound via Thm 7.2 + Thm 6.3 | PROVED HERE | — |
| clean four-insertion width bound | — | **NOT CLAIMED** (open) | — |
| unrestricted garbage-model two-query bound | — | **NOT CLAIMED** (open; equals the row above by Thm 7.2) | — |
| restricted garbage-model two-query bound (junk register below $`2n-O(1)`$ qubits) | — | **NOT CLAIMED** (in working records, not audited) | — |
| three or more adaptive queries | — | **NOT CLAIMED** | — |

## The audit of record

The two-query width theorem in its real-selector, symmetric-weighting form (constant 8)
passed an independent gate audit on 2026-08-16 in the author's private research
repository. That audit re-derived the six proof steps, checked the trace budget, the
pointwise absorption, and the Tropp application, and additionally verified the asymmetric
$`\sqrt{LR}`$ normalization for real selectors. The extension to independent unimodular
selectors (Theorem 4.1) reuses the audited asymmetric route and adds only the observation
that the leverage-weighted selector norms are constant on the torus; it has been checked
numerically here but not yet audited.

## What changed relative to the earlier draft

1. The unitarity criterion assumed an isometric middle layer; corrected (Prop 2.2).
2. The earlier draft's summary claimed coverage of diagonal-unitary oracles while its
   remark disclaimed it. Resolved in favour of coverage, by the unimodular route
   (Theorem 4.1, Remark 4.6).
3. The audited constant 8 is retained for the audited form; the unimodular route gives 2.
4. New: address-length bound (Thm 6.3), matching upper bound (Thm 6.4), thresholds
   (Cor 6.5), garbage-model transfer (Section 7).

## Lean

Five statements in `lean/TwoQuery/Q2Width.lean` are kernel-checked with axioms
`[propext, Classical.choice, Quot.sound]` only. The analytic inputs (trace budget,
pointwise absorption, matrix Khintchine) are hypotheses of the statements, not axioms and
not proved. `lean/RECORDED-AXIOMS.txt` holds the recorded `#print axioms` output. See
`lean/README.md`.

## Reproducing the checks

```
python -u verify/run_all.py
python verify/md_math_lint.py paper/paper.md notes/*.md README.md STATUS.md
cd lean && lake build && lake env lean axiom_check.lean
```

The verifier suite needs only Python 3 and numpy and finishes in a few minutes. The Lean
check needs the pinned toolchain in `lean/lean-toolchain` and a Mathlib build cache.
