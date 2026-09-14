# STATUS — claim-by-claim ledger

*Updated 2026-09-13 for version 1.0.4. This file records support and scope.
The self-contained mathematical manuscript is `paper/paper.tex`; its Markdown
companion and this ledger must agree with it. See `VALIDATION.md` for the current
local review. Historical audit labels below are preserved, not newly certified.*

The fresh-context independent AI audit of 2026-09-13 found no headline proof blocker
after the quantifier and finite-count corrections. See
`audit/2026-09-13-independent-audit.md`. It is an additional mathematical review,
not external peer review or complete formal verification. Current release checks
and publication status are recorded in `publication/RELEASE-CHECKLIST.md`.

## Labels

| label | meaning |
|---|---|
| **AUDITED** | historical independently re-derived argument on the date given, supplemented by finite checks; not a full machine formalization |
| **PROVED HERE** | analytic proof supplied and locally reviewed; finite instances checked where applicable, not an external peer-review certification |
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
| Prop 3.2 | three-insertion pair, $`Q=8,N=1`$: 2 versus 16 scalar matrices, not distinct channels | PROVED HERE (exact finite certificate) | Appendix B; `verify_exact_controls.py` exhausts all 256 words in rational arithmetic |
| Thm 4.1 | two-query width with independent unimodular selectors, constant 2 | PROVED HERE | `verify_complex_phase.py` (a)–(d), (f) |
| Cor 4.2 | real sign words, constant 2 (as a corollary of Thm 4.1) | PROVED HERE | same |
| Cor 4.2, audited form | real sign words, symmetric weighting, constant 8 | **AUDITED 2026-08-16**; the asymmetric $`\sqrt{LR}`$ route was verified in the same audit for real selectors | `verify_proof_steps.py`; **CERTIFIED (Lean)** for the deterministic spine |
| Cor 5.1 | capacity in the normalized Frobenius metric only | AUDITED 2026-08-16 (as a corollary of the audited width bound) | metric erratum recorded in the paper |
| Cor 5.2 | quantitative packing of $`U(N)`$ and the $`\Omega(2^n)`$ qubit bound | AUDITED 2026-08-16 for the bound; the quantitative packing lemma is PROVED HERE | — |
| Prop 5.3 | the permutation family forces $`\kappa\ge1/4`$; no numerical upper bound on the universal Sudakov constant is claimed | PROVED HERE | `verify_crho_requirement.py` (finite table values and calibration only) |
| Prop 6.1 | the Dong–Lombardi–Ma permutation family lies in the class and attains the ceiling up to the constant | PROVED HERE | `verify_permutation_family.py` |
| Cor 6.2 | both factors tight at $`Q=\Theta(N^2)`$; the logarithm is necessary at exponential $`Q`$ | PROVED HERE | — (Sudakov backwards on Prop 6.1; Thm 6.4) |
| Thm 6.3 | address-length lower bound $`\Omega(2^n)`$ bits, any workspace, two decompositions allowed | PROVED HERE (new in 1.0.0) | `verify_address_compression.py` |
| Thm 6.4 | permutation transport: $`O(2^n\log(1/\varepsilon))`$ qubits and address bits, two Boolean queries, clean | PROVED HERE; complete proof now in PDF Appendix A | `verify_koopman_transport.py` (1)–(5); exact non-scalar $`N=2`$ controls in `verify_exact_controls.py` |
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

Version 1.0.3 additionally received a local proof/dependency review of the full
manuscript, not a blind or external gate. The headline clean two-query thresholds
survive that review. The unsupported minimal-degree claim was withdrawn; the
claimed numerical Sudakov constant was removed; attribution/model wording was
corrected. The detailed result and reproducibility limits are in `VALIDATION.md`.

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
   replaced by a claimed bound $`C_S\le2.09`$ and a quantitative packing lemma for $`U(N)`$.
   **Superseded in 1.0.3:** the universal numerical-constant claim lacked a proof
   beyond finitely many quadratures and has been removed. The packing lemma remains.
8. Two verifiers had no assertions; two others tested Lemma 2.1 of the transport note only
   where it is vacuous or formed the coupling bound only synthetically. All four now assert,
   and the transport verifier runs the construction end to end at $`N=1`$.
9. A tab byte had replaced `\times` in one formula; the markdown linter now rejects tabs.
10. Overclaims in the abstract, introduction and README ("exactly", "nothing was known",
    "for every oracle") were reworded to what the theorems support.

## Referee pass on the compiled paper (version 1.0.2)

A sixth audit refereed the LaTeX rendition and the PDF against the markdown and the five
earlier reports. It found no new mathematical error and confirmed every round-one fix.
Changes made in 1.0.2: the capacity constant is now $`\kappa`$, the quantizer level count
$`J`$, the workspace dimension of Section 7 $`k`$, and the transport construction uses
$`\theta`$, $`\alpha,\beta`$, $`\Xi`$ for the cell map, labels and coupling, so that no
symbol in the paper has two meanings; the abstract's junk-blindness sentence is restricted
to linear width arguments, as Lemma 7.1 proves; Theorem 6.3 states its error tolerance;
Corollary 6.5's error is "below an absolute constant"; the Meckes citation gives the
theorem number; Section 8 records that any $`N\cdot\mathrm{poly}(k,\log Q)`$ width bound
would already answer the Aaronson–Kuperberg question negatively, and that Theorem 7.2
carries this to the garbage model; the paper names its repository URL; the condensed
construction no longer ends with a proof mark.

## What changed relative to the earlier private draft (version 1.0.0)

1. The unitarity criterion assumed an isometric middle layer; corrected.
2. The earlier draft's summary claimed coverage of diagonal-unitary oracles while its
   remark disclaimed it. Resolved in favour of coverage, by the unimodular route
   (Theorem 4.1, Remark 4.6).
3. The audited constant 8 is retained for the audited form; the unimodular route gives 2.
4. New: address-length bound (Thm 6.3), matching upper bound (Thm 6.4), thresholds
   (Cor 6.5), garbage-model transfer (Section 7).

## Lean

Seven statements in `lean/TwoQuery/Q2Width.lean` are kernel-checked with axioms
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

`paper/paper.tex` is the self-contained submission source, with embedded bibliography
and complete proof appendices. The Markdown paper and notes are reading companions.
The main mathematical claims must agree; appendix proofs need not be duplicated in
both formats. `tools/build_paper.py` builds in a clean temporary directory and creates
`release/arxiv-source.zip` containing just `paper.tex`. Generated bundles and render
images are ignored; no release or upload is performed by this script.
