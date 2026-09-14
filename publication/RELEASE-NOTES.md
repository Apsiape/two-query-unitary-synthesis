# Version 1.0.4 — public source publication

Public GitHub publication and versioned archival release were authorized on
2026-09-13; the author confirmed that the Zenodo integration is enabled.
The public Zenodo record and archived files have been verified:
https://doi.org/10.5281/zenodo.22740995 (version 1.0.4).
The all-versions DOI is https://doi.org/10.5281/zenodo.22740994.

The Clean Two-Query Threshold for Unitary Synthesis, Seth Douglas.

This preprint establishes the clean two-query width bound, workspace-independent
address-length lower bound, and matching fixed-error resource thresholds. It does
not settle the unrestricted garbage-model two-query problem or full AK.

This release includes the self-contained manuscript, its complete proof appendices,
two companion notes, nine finite-check scripts and seven scoped Lean declarations.
The full analytic proof is not formalized. An independent, fresh-context AI audit
is documented in `audit/2026-09-13-independent-audit.md`; this is not external peer review.

Changes since 1.0.3:

- Completed the explicit circuit-to-matrix reduction and clean-output error
  implication; clarified the clean-only scope of the conjectural conclusion.
- Distinguished the real-linear norm from the complex covariance-map norm in
  transport, without changing the valid covariance bound.
- Removed duplicated prose and internal shorthand, explained why easy permutation
  families do not imply universality, and refreshed current snapshot provenance.
- Expanded AI-use disclosure, preserving Q2-specific
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

Licensing is file-specific: manuscript and prose are CC BY 4.0; software is MIT.
The root `.zenodo.json` supplies preprint metadata and explicitly records this
distinction; it overrides Zenodo's limited CITATION.cff importer. The manuscript
PDF and TeX retain the hashes recorded in `VALIDATION.md`.

The final adversarial AI audit found no remaining mathematical or artifact-integrity
release blocker within its inspected scope; its record is
`audit/2026-09-13-final-adversarial-audit.md`. It is not external peer review.
