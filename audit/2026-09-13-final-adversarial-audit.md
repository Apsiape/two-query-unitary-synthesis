# Final adversarial AI audit — 2026-09-13

## Verdict and scope

PASS within the inspected scope: no remaining mathematical or artifact-integrity
release blocker found. Critical/high findings: none. No required repair identified.
This was an independent, fresh-context, read-only AI review of the corrected
manuscript, companion notes, release documentation, and artifacts. Prior audits
and test success were not accepted as proof. This report records the auditor's
returned findings; it is not external human peer review, a complete formalization,
or an exhaustive priority certification.

The auditor made no paper/code changes and did not commit, push, tag, or publish.
Public publication is a separate author decision.

## Independently inspected arguments

1. **Circuit reduction:** answer-bit diagonalization, fixed input/output
   embeddings, spectator registers, fixed contractions, and the projection of
   full-output approximate-clean error. Reduced-channel closeness alone is not
   asserted to imply the same full-isometry error.
2. **Gaussian width:** exact conjugations, zero-leverage supports, pointwise
   torus absorption, both variance operators, Schur domination for arbitrary
   middle contractions, and the conservative complex-Gaussian constant 2.
3. **Packing and address compression:** invariant endpoint spaces of dimensions
   at most MN, preserved Frobenius budgets, arbitrary address-projector ranks,
   the Boolean M+1 eigenspace treatment, complex-Gaussian metric normalization,
   and the reference concentration statement.
4. **Finite transport:** dimension-independent quantization, fixed encoding,
   finite Birkhoff decomposition, coherent rounding with target-independent
   label dimension, covariance orientation, and the actual tagged two-query
   forward/inverse permutation implementation. The real-linear squared norm
   1/(2s^2) is distinguished from LL*=I/s^2 for the complex covariance map.
   Resource counts and full-output cleanup bounds were checked.
5. **Garbage transfer:** packing within a fixed architecture, workspace-uniform
   hypothesis, same-word 2t-insertion representation, normalized junk cancellation,
   channel-duality error control, and a packing on the balanced-reflection orbit.
   Unrestricted garbage-model two-query synthesis remains open. The all-depth
   clean conjecture is distinguished from the address-uniform garbage premise.
6. **Tightness:** permutation packing logarithm Theta(N log N) is compatible
   with the Theta(N^2) universal-target packing requirement. The kappa calibration
   is justified asymptotically; finite floating-point checks are not claimed as
   certified analytic constants.

Primary-source comparisons confirmed the uses of
[Tropp, Theorem 4.1.1](https://arxiv.org/pdf/1501.01571) and
[Meckes, Theorem 5.17](https://case.edu/artsci/math/mwmeckes/elizabeth/Haar_book.pdf).
The latter's N-2 coefficient yields an absolute multiple of N for N>=3.

## Executed checks and exact snapshot

- 62 targeted circuit/covariance controls passed; maximum residual 7.11e-15.
- Release/source consistency and Markdown lint passed.
- All 44 then-current source-manifest entries matched disk and source ZIP.
- Every then-current SHA256SUMS artifact matched.
- The isolated arXiv-source ZIP contained only paper.tex, byte-identical to source.
- `git diff --check` passed. Diagnostic jobs were capped and all exited.

Verified TeX SHA-256:
`B9BF7D21198A032E270A1322F4A590268C523975D03CC7CEFFFFE44F7A819BB8`.

Verified PDF SHA-256:
`1529C00E91C5153B8F03FBAEF6283F2C149CE9294574170FAA5F45BED79E508D`.

The 44-file manifest count refers to the audited snapshot before adding this
report and updating public-publication bookkeeping. Subsequent source packages
must be regenerated; the unchanged manuscript hashes identify the audited proof.

## Limits

The auditor did not rebuild/render the PDF, rerun the full nine-script suite or
Lean, or conduct an exhaustive novelty search. Earlier visual QA and kernel
checks remain historical evidence, not new checks by this auditor. No absolute
guarantee of correctness or publication acceptance follows from this pass.
