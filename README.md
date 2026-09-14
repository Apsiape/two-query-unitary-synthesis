# The Clean Two-Query Threshold for Unitary Synthesis

**Seth Douglas** · version 1.0.4 · 2026-09-13

A focused research repository: one self-contained paper, two companion notes, nine
finite-check scripts, and seven scoped Lean declarations. The complete analytic proof
is not formalized. See [VALIDATION.md](VALIDATION.md) for the current readiness assessment.

Read the [manuscript PDF](paper/paper.pdf), the
[independent AI audit](audit/2026-09-13-independent-audit.md), and the
[release checklist](publication/RELEASE-CHECKLIST.md). This audit is not external peer review.

## The result in three lines

Write $`N=2^n`$. For clean two-query unitary synthesis of every $`n`$-qubit unitary:

- **lower bound** — any fixed two-query architecture, with any number of ancillas and any
  diagonal-unitary oracle, realizes a family of Gaussian mean width at most
  $`2\min(N,Q)\sqrt{\log 2Q}`$; hence universality needs $`\Omega(2^n)`$ qubits and, by an
  address-compression argument that ignores workspace entirely, $`\Omega(2^n)`$ oracle
  input bits;
- **upper bound** — a permutation-transport construction implements every unitary cleanly
  with two Boolean queries to one table, $`O(2^n\log(1/\varepsilon))`$ qubits and
  $`O(2^n\log(1/\varepsilon))`$ address bits;
- so both thresholds are $`\Theta(2^n)`$ up to the factor $`\log(1/\varepsilon)`$.

The paper also shows the width bound is tight in both factors at register dimension $`Q=\Theta(N^2)`$,
explains why the one-query argument does not extend, and reduces the garbage-model
two-query problem to a clean four-insertion width bound, which remains open.

## Layout

```
paper/paper.md                  the paper (GitHub-rendered mathematics)
paper/paper.tex, paper.pdf      self-contained paper, including complete proof appendices
notes/permutation-transport.md  the matching upper bound, in full
notes/garbage-to-clean.md       the garbage-model transfer, in full
verify/run_all.py               runs all nine finite-check scripts; exit 0 iff all pass
verify/verify_*.py              one deterministic numpy script per group of claims
verify/md_math_lint.py          lints the markdown mathematics
lean/                           Lean 4 + Mathlib: kernel-checked deterministic spine
STATUS.md                       claim-by-claim ledger: audited / proved here / not claimed
CHANGELOG.md                    version history
VALIDATION.md                   current proof/dependency review and release scope
tools/build_paper.py            isolated build and local arXiv-source bundle
tools/check_release.py          non-mutating source consistency checks
CITATION.cff                    how to cite
```

## Reproduce

```
python -m pip install -r requirements.txt
python -u verify/run_all.py
```

Expected last line: `SUITE PASSED -- 9/9 verifiers green`. Each script prints a PASS/FAIL
line per check and exits nonzero on any failure. All seeds are fixed. The whole suite runs
in about a minute on a laptop and uses well under a gigabyte of memory.

For the Lean check:

```
cd lean
lake build
lake env lean axiom_check.lean
```

Each axiom list must be a subset of `[propext, Classical.choice, Quot.sound]`;
in particular, no `sorryAx` or custom axiom is allowed.
The toolchain is pinned in `lean/lean-toolchain` and the Mathlib revision in
`lean/lakefile.toml`; a first build downloads the Mathlib cache. What is and is not
certified is stated in `lean/README.md`.

To build in a clean temporary directory, check references/layout warnings, and
produce the local submission source bundle:

```
python tools/check_release.py
python tools/build_paper.py
```

## Reading order

1. `paper/paper.md` Sections 1–2 for the model, then Section 4 for the width theorem.
2. `STATUS.md` before quoting any statement: it says which claims have passed an
   independent audit, which are proved here and self-verified, and what the audit round
   of 2026-09-11 found and fixed.
3. Appendix A of the PDF for the complete matching construction, or its Markdown note.
   Appendix B specifies the exact finite control. The notes are not required to compile
   or read the complete submitted proof.

## Relation to prior work

At one query the problem was settled by Aaronson–Kuperberg (exact model), Lombardi–Ma–Wright
(approximate model with garbage allowed, also against polynomially many parallel queries),
Dong–Lombardi–Ma (explicit separations) and Brakerski–Yuen (input length at any query
count); Rosenthal gave the $`O(2^{n/2})`$-query upper bound. Banerjee's every-oracle-clean
analysis and Huang's 1.5-query study are also acknowledged in the manuscript. This
repository establishes matching fixed-error resource thresholds for two sequential clean
queries. The clean model is strictly smaller than the garbage model, and
the results here neither improve nor contradict the one-query garbage-model bounds;
Section 7 of the paper prices the distance between the two models.

## AI-use disclosure

Language models were used extensively for proof discovery, implementation,
verification engineering, refutation-first auditing, and editing. The author
directed the work and assumes sole responsibility. AI audits are not external
peer review, and the complete analytic proof is not formalized. The full
[disclosure](paper/paper.md#acknowledgements-and-disclosure) states the scope.

## License

Text under CC BY 4.0 (`LICENSE`); code under MIT (`LICENSE-CODE`).
