# Validation and release scope — version 1.0.4

Version 1.0.4 adds a fresh-context independent AI audit, documented in
`audit/2026-09-13-independent-audit.md`. This is not external peer review or a
formalization of the full paper. It found no headline proof blocker and identified
two corrections: the transfer theorem must fix the architecture before packing
instructions, and the finite permutation-count check must not substitute n!/e
for exact derangements. Both are repaired; the displayed finite table is unchanged.
The audit also sharpened the rate-function and numerical-verification wording.
See `publication/RELEASE-CHECKLIST.md` for the current build and publication state.

The following records the prior version 1.0.3 local proof/dependency review,
not a fresh certification of every historical claim. The reviewer read the main manuscript and
both notes, re-derived the load-bearing arguments, inspected the formalization,
checked the cited primary statements, and added exact finite controls.

## Mathematical assessment

No unresolved defect was found in the headline clean two-query threshold after
the corrections below. At fixed sufficiently small error, the lower and upper
bounds match in both total qubits and Boolean oracle input length. The latter
lower bound is independent of retained workspace. The result does not settle
two-query synthesis with unrestricted target-dependent garbage or full AK.

| Dependency | What was checked |
|---|---|
| Width theorem | Selector-independent trace expansion; complex Gaussian normalization; both leverage supports; constant weighted selector norms; both Schur square functions; compression by correlation matrices; Tropp's rectangular Gaussian bound. |
| Address compression | Orthogonal address blocks preserve each compressed endpoint range; dimensions at most MN; Frobenius norms are unchanged; the middle compression is contractive; allowing independent compressed phases only enlarges the supremum. |
| Capacity and lower bound | Real Gaussian-process metric is Frobenius distance divided by square root of two; Sudakov is used at normalized scale; Haar packing follows from an explicit Lipschitz/concentration argument. |
| Transport | Dimension-independent quantizer error; exact encoding isometry; uniform coupling marginals; finite Birkhoff decomposition; covariance bound without residual independence; correct conjugation by U; fixed label-size rounding; operator-norm rather than average-input error. |
| Physical upper-bound implementation | Fixed encoder and inverse, BV forward/inverse lookups of one tagged Boolean table, complete address count, restoration of direction/answer/scratch, full-isometry cleanup and reference-uniform error. Offline table construction is unrestricted, not efficiently compiled. |
| Garbage transfer | The doubled product is an insertion-class representation, not necessarily a physically clean unitary. Junk cancels in observables. Direct channel duality controls error; the balanced-reflection orbit packing is now justified explicitly. |
| Finite scalar example | An explicit rational-phase construction exhausts all 256 sign words exactly. It counts matrices at N=1, not distinct physical channels, and is not used in the threshold. |

Zero endpoints make the width identically zero; otherwise each supported matrix
dimension is positive. The log-normalization and support restrictions remain
unchanged. The matching upper bound is approximate, not an exact finite-table
implementation of a continuum of unitaries.

## Corrections made

1. Junk-invariant linear functionals vanish, but that alone does not imply a
   minimal oracle degree or optimal transfer depth. The stronger claim was
   withdrawn everywhere it was used as an active conclusion. A cancellation
   control illustrates why packet degree cannot be inferred that way.
2. The claimed universal Sudakov constant 2.09 was not justified by its finite
   quadratures. Removed the numerical upper claim and its misleading assertions.
   The standard theorem with an unspecified absolute constant is sufficient.
3. ABP's cited characterization is for scalar expected outputs of t-query
   algorithms and degree-2t forms. It was not a proved equivalence with the
   matrix-valued amplitude class asserted in the earlier prose.
4. AK's relevant Boolean table entries were not oracle address bits. Corrected
   that terminology, and added the primary Banerjee report and Huang preprint.
5. Added complete proof appendices to the submission PDF and narrowed the
   blanket machine-checking statement to the actual scope.

These are substantive proof-scope and exposition corrections; they should not be
described as proof of the previously stronger claims. None changes the headline
clean two-query thresholds.

## Primary-source checks

- [Tropp, Theorem 4.1.1](https://tropp.caltech.edu/books/Tro14-Introduction-Matrix-preprint.pdf):
  rectangular Gaussian expectation bound with both variance operators and a
  logarithmic dimension factor. The conservative complex-series constant is valid.
- [Meckes, Theorem 5.17](https://case.edu/artsci/math/mwmeckes/elizabeth/Haar_book.pdf):
  concentration on U(N) in Hilbert–Schmidt distance. The stated N-minus-two
  coefficient supplies the paper's absolute constant for N at least three.
- [LMW](https://arxiv.org/html/2310.08870v1):
  one-query/discard and parallel-query context; the present two-query-clean
  theorem is not an unrestricted two-query-garbage result.
- [Dong–Lombardi–Ma](https://arxiv.org/html/2607.26478v1):
  permutation implementation and explicit one-query separations; Question 1.4
  has its own alternating-phase indexing.
- [Rosenthal, Theorems 2.2 and 3.1–3.2](https://arxiv.org/html/2111.07992v5):
  distinguish the quantum column-oracle theorem from its Boolean-oracle compiler.
- [ABP, Theorem 1.3](https://arxiv.org/html/1711.07285v2):
  scalar expectation characterization, not the former unqualified matrix-valued claim.
- [AK, Theorem 6.7](https://arxiv.org/html/quant-ph/0604056):
  the every-oracle unitarity condition and Boolean-variable counting.
- [Banerjee](https://www.scottaaronson.com/showcase5/synthesis.pdf):
  exact every-oracle-clean model, improved one-query count and restricted
  two-query analysis; not the present approximate address theorem.
- [Huang](https://arxiv.org/html/2508.13215v2):
  related 1.5-query study, not a dependency here. No correctness endorsement of
  that paper's separate estimates is made.
- [Brakerski–Yuen](https://arxiv.org/html/2605.09957v1):
  address-length results with different query scope.

This was a targeted comparison, not a proof that no equivalent result exists
anywhere. Broad priority superlatives have been removed from the abstract.

## Reproducibility and release boundary

The finite suite, exact controls, source checks, scoped Lean declarations, isolated
PDF compilation and visual rendering are separate checks. Passing one does not
stand in for another or certify the asymptotic analytic proof.

The seven Lean declarations retain explicit analytic hypotheses where indicated
in `lean/README.md`. The full Gaussian concentration theorem, Schur variance
derivation, packing theory, and transport construction are not Lean-formalized.
The local build initially timed out while importing all Mathlib; targeted imports
resolved that problem. No uncapped retry was used.

`tools/build_paper.py` compiles in a fresh temporary directory and creates a
local `release/arxiv-source.zip` containing just `paper.tex`. Bibliography and
appendices are embedded. The archive is generated output, not a release.
Nothing in this pass commits, pushes, mints a DOI, or submits to arXiv.

Before publication, the author should read the final PDF and confirm the title,
attribution and disclosures. This document does not promise mathematical infallibility.

## Executed local checks

- `python verify/run_all.py`: **9/9 passed**, including the new rational-arithmetic
  controls. Numerical checks remain finite-instance evidence.
- `lake build` and `lake env lean axiom_check.lean` in `lean/`: passed;
  all seven declarations use only a subset of `propext`, `Classical.choice`,
  and `Quot.sound`. This does not discharge explicit analytic hypotheses.
- `python tools/check_release.py`: passed source/version/label/scope consistency
  and Markdown mathematics checks.
- `python tools/build_paper.py`: passed isolated three-pass compilation of only
  `paper.tex`, with embedded bibliography and appendices, no undefined references
  and no overfull boxes. The PDF is **17 pages**.
- All pages were rendered and visually inspected, with individual-page checks of
  the title page and dense proof appendices. No clipping or layout defect remained.
- `git diff --check`: passed. Generated release and rendering directories are ignored.

Final reviewed source SHA-256:
`E39591FDCFB5F31A95A65DD9D06C914797F4E144F1D8649E64F509F3C293C9BB`.

Final locally built PDF SHA-256:
`BBDEF52AA80FA533BEF33D4A3507F1C02DC430D6A18839FF9BD711E940B20016`.
Rebuilding may change PDF metadata and its hash without changing the source.
