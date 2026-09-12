# STATUS — claim-by-claim ledger

*Updated 2026-09-11 for version 1.0.1. This file is the authority on what this repository
claims and how each claim is supported. When it disagrees with a sentence elsewhere, this
file wins.*

## Labels

| label | meaning |
|---|---|
| **AUDITED** | proved in the paper, machine-checked here, and independently re-derived by an audit external to the author on the date given |
| **PROVED HERE** | proved in the paper and machine-checked by a script in `verify/` where a finite check exists; re-derived in the 2026-09-11 refutation-first audit round (see below) but not by an auditor external to the project |
| **CERTIFIED (Lean)** | the deterministic part is kernel-checked in `lean/`, with the analytic inputs as explicit hypotheses; see `lean/README.md` |
| **NOT CLAIMED** | discussed for context only; no theorem asserted |

The verifier scripts are necessary, not sufficient: they catch algebraic and numerical
mistakes at small size and cannot certify an asymptotic argument. A "—" in the machine-check
column means the claim is asymptotic or a direct consequence of checked claims and admits
no finite check of its own.

## Ledger

| # | claim | label | machine check |
|---|---|---|---|
| Prop 2.2 | three-term unitarity criterion for contractive $`C`$ and $`W`$ | PROVED HERE | `verify_complex_phase.py` (e); isometric case in `verify_proof_steps.py` |
| Prop 3.1 | at depth two the final layer does not change the number of distinct targets | PROVED HERE | `verify_decisive_pair.py` stage 0 (exhaustive, $`Q=8`$, $`N=2`$, six final layers) |
| Prop 3.2 | decisive pair in the **three-insertion** class, $`Q=8`$, $`N=1`$: identical pre-final data, 2 versus 16 targets | PROVED HERE (exhaustive) | `verify_decisive_pair.py` stage 1 |
| Thm 4.1 | two-query width with independent unimodular selectors, constant 2 | PROVED HERE | `verify_complex_phase.py` (a)–(d), (f) |
| Cor 4.2 | real sign words, constant 2 (as a corollary of Thm 4.1) | PROVED HERE | same |
| Cor 4.2, audited form | real sign words, symmetric weighting, constant 8 | **AUDITED 2026-08-16**; the asymmetric $`\sqrt{LR}`$ route was verified in the same audit for real selectors | `verify_proof_steps.py`; **CERTIFIED (Lean)** for the deterministic spine |
| Cor 5.1 | capacity in the normalized Frobenius metric only | AUDITED 2026-08-16 (as a corollary of the audited width bound) | metric erratum recorded in the paper |
| Cor 5.2 | quantitative packing of $`U(N)`$ and the $`\Omega(2^n)`$ qubit bound | AUDITED 2026-08-16 for the bound; the quantitative packing lemma is PROVED HERE | — |
| Prop 5.3 | the permutation family forces $`K\ge1/4`$; $`C_S\le2.09`$, $`K\le35`$ | PROVED HERE | `verify_crho_requirement.py` (table values, monotonicity, the Sudakov–Fernique constant) |
| Prop 6.1 | the Dong–Lombardi–Ma permutation family lies in the class and attains the ceiling up to the constant | PROVED HERE | `verify_permutation_family.py` |
| Cor 6.2 | both factors tight at $`Q=\Theta(N^2)`$; the logarithm is necessary at exponential $`Q`$ | PROVED HERE | — (Sudakov backwards on Prop 6.1; Thm 6.4) |
| Thm 6.3 | address-length lower bound $`\Omega(2^n)`$ bits, any workspace, two decompositions allowed | PROVED HERE (new in 1.0.0) | `verify_address_compression.py` |
| Thm 6.4 | permutation transport: $`O(2^n\log(1/\varepsilon))`$ qubits and address bits, two Boolean queries, clean, repeated-sign-word form | PROVED HERE (new in 1.0.0; full construction in `notes/permutation-transport.md`) | `verify_koopman_transport.py` (1)–(5), including an end-to-end run at $`N=1`$ |
| Cor 6.5 | both thresholds $`\Theta(2^n)`$ up to $`\log(1/\varepsilon)`$ | PROVED HERE | follows from Cor 5.2, Thm 6.3, Thm 6.4 |
| Note §7.1 | monomial-action converse $`(N-1)\log_2(1/\varepsilon)`$ for fixed encodings | PROVED HERE | — |
| Lemma 7.1 | junk-blind linear functionals vanish (complex- and real-linear) | PROVED HERE | `verify_conjugation_transfer.py` (4) |
| Thm 7.2 | conjugation transfer: clean-$`2t`$ capacity in address form ⇒ garbage-$`t`$ impossibility | PROVED HERE | `verify_conjugation_transfer.py` (1)–(3) |
| Cor 7.3 | one-query garbage-model address-length bound via Thm 7.2 + Thm 6.3, constant error | PROVED HERE | — |
| clean four-insertion width bound | — | **NOT CLAIMED** (open) | — |
| unrestricted garbage-model two-query bound | — | **NOT CLAIMED** (open; reduces to the row above by Thm 7.2, one direction only) | — |
| restricted garbage-model two-query bound (junk register below $`2n-O(1)`$ qubits) | — | **NOT CLAIMED** (in working records, not audited) | — |
| three or more adaptive queries | — | **NOT CLAIMED** | — |

## The audit of record

The two-query width theorem in its real-selector, symmetric-weighting form (constant 8)
passed an independent gate audit on 2026-08-16 in the author's private research
repository. That audit re-derived the six proof steps, checked the trace budget, the
pointwise absorption, and the Tropp application, and additionally verified the asymmetric
$`\sqrt{LR}`$ normalization for real selectors. The extension to independent unimodular
selectors (Theorem 4.1) reuses the audited asymmetric route and adds only the observation
that the leverage-weighted selector norms are constant on the torus.

## Audit round of 2026-09-11 (version 1.0.1)

Five independent refutation-first audits were run on version 1.0.0, one per theorem group
plus one on repository consistency and literature. No claim was found false. The defects
found and fixed in 1.0.1:

1. Proposition 2.2 stated the unitarity criterion with two terms; for a contractive $`C`$ a
   third term $`I-C^{*}C`$ is needed (counterexample $`C=0`$). Fixed; verifier (e) now uses
   contractive $`C`$ and $`W`$.
2. The decisive pair of Section 3 was verified in a three-insertion model with a
   one-dimensional target, not the two-query model the text implied; and the phenomenon
   provably cannot occur at two queries (new Proposition 3.1). Section 3 rewritten.
3. Corollary 6.2's closing sentence left open a refinement that Theorem 6.4 excludes.
   Replaced by the correct statement.
4. Theorem 7.2's hypothesis was stated with $`\mathrm{polylog}(Q)`$ and on the register
   dimension rather than the address count; the conclusion did not follow as written, and
   "for every junk dimension" needed the address form of Theorem 6.3. Restated.
5. The Stinespring-continuity citation: [KSW08] gives the junk-register unitary on an
   enlarged environment; the version on the junk register itself is [vE23]. Both now cited
   and the enlargement argued.
6. Literature: the poly(n)-parallel-query bound is due to [LMW24], not [DLM26]; the
   translation to [DLM26]'s indexing was stated backwards; [Ros26]'s query result uses
   poly(n) ancillas, not $`2^{O(n)}`$; [AK07]'s footnote 13 concerns the every-oracle
   unitarity requirement, not cleanliness. All corrected, and the model paragraph of
   Section 1.2 now states that no requirement is imposed on non-lawful oracles.
7. The Sudakov constant: the earlier draft's "constants in $`[2,4]`$" was unsupported;
   replaced by a derived bound $`C_S\le2.09`$ and a quantitative packing lemma for $`U(N)`$.
8. Two verifiers had no assertions; two others tested Lemma 2.1 of the transport note only
   where it is vacuous or formed the coupling bound only synthetically. All four now assert,
   and the transport verifier runs the construction end to end at $`N=1`$.
9. A tab byte had replaced `\times` in one formula; the markdown linter now rejects tabs.
10. Overclaims in the abstract, introduction and README ("exactly", "nothing was known",
    "for every oracle") were reworded to what the theorems support.

## What changed relative to the earlier private draft (version 1.0.0)

1. The unitarity criterion assumed an isometric middle layer; corrected.
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

The verifier suite needs only Python 3 and numpy and finishes in about a minute. The Lean
check needs the pinned toolchain in `lean/lean-toolchain` and a Mathlib build cache.
