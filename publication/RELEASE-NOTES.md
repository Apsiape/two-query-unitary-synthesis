# Version 1.0.4 — release candidate

Private repository preparation only; public release and Zenodo publication are on hold.

The Clean Two-Query Threshold for Unitary Synthesis, Seth Douglas.

This preprint establishes the clean two-query width bound, workspace-independent
address-length lower bound, and matching fixed-error resource thresholds. It does
not settle the unrestricted garbage-model two-query problem or full AK.

This release includes the self-contained manuscript, its complete proof appendices,
two companion notes, nine finite-check scripts and seven scoped Lean declarations.
The full analytic proof is not formalized. An independent, fresh-context AI audit
is documented in `audit/2026-09-13-independent-audit.md`; this is not external peer review.

Changes since 1.0.3:

- Expanded AI-use disclosure in line with the I3322 paper, preserving Q2-specific
  verification limits and explicit author responsibility.

- Made the fixed-architecture quantifier explicit in the conjugation-transfer
  theorem; its capacity bound packs instructions, not architectures.
- Allowed the rate function to depend on target dimension, but not workspace.
- Replaced an approximate derangement count with the exact recurrence, evaluated
  in floating point and clearly distinguished from an interval certificate.
- Distinguished clean-output compression from a reduced quantum channel and
  numerical local optima from a certified torus supremum.
- Added manuscript-preferred citation and Zenodo preprint metadata; clarified the
  separate manuscript (CC BY 4.0) and source-code (MIT) licenses.

See `VALIDATION.md` for the checks actually completed and their limits.
