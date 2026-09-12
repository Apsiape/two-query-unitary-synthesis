# The Clean Two-Query Threshold for Unitary Synthesis

**Seth Douglas** · version 1.0.2 · 2026-09-11

A self-contained repository: one paper, two companion notes, eight machine verifiers,
and a kernel-checked Lean spine.

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

The paper also shows the bound is tight in both factors at workspace $`\Theta(N^2)`$,
explains why the one-query argument does not extend, and reduces the garbage-model
two-query problem to a clean four-insertion width bound, which remains open.

## Layout

```
paper/paper.md                  the paper (GitHub-rendered mathematics)
paper/paper.tex, paper.pdf      the same paper as LaTeX source and compiled PDF
notes/permutation-transport.md  the matching upper bound, in full
notes/garbage-to-clean.md       the garbage-model transfer, in full
verify/run_all.py               runs all eight verifiers; exit 0 iff all pass
verify/verify_*.py              one deterministic numpy script per group of claims
verify/md_math_lint.py          lints the markdown mathematics
lean/                           Lean 4 + Mathlib: kernel-checked deterministic spine
STATUS.md                       claim-by-claim ledger: audited / proved here / not claimed
CHANGELOG.md                    version history
CITATION.cff                    how to cite
```

## Reproduce

```
pip install numpy
python -u verify/run_all.py
```

Expected last line: `SUITE PASSED -- 8/8 verifiers green`. Each script prints a PASS/FAIL
line per check and exits nonzero on any failure. All seeds are fixed. The whole suite runs
in about a minute on a laptop and uses well under a gigabyte of memory.

For the Lean check:

```
cd lean
lake build
lake env lean axiom_check.lean
```

Every line printed must read `depends on axioms: [propext, Classical.choice, Quot.sound]`.
The toolchain is pinned in `lean/lean-toolchain` and the Mathlib revision in
`lean/lakefile.toml`; a first build downloads the Mathlib cache. What is and is not
certified is stated in `lean/README.md`.

To rebuild the PDF from source (standard packages only; three passes resolve the
cross-references):

```
cd paper
pdflatex paper.tex && pdflatex paper.tex && pdflatex paper.tex
```

## Reading order

1. `paper/paper.md` Sections 1–2 for the model, then Section 4 for the width theorem.
2. `STATUS.md` before quoting any statement: it says which claims have passed an
   independent audit, which are proved here and self-verified, and what the audit round
   of 2026-09-11 found and fixed.
3. The notes for the two constructions the paper only condenses.

## Relation to prior work

At one query the problem was settled by Aaronson–Kuperberg (exact model), Lombardi–Ma–Wright
(approximate model with garbage allowed, also against polynomially many parallel queries),
Dong–Lombardi–Ma (explicit separations) and Brakerski–Yuen (input length at any query
count); Rosenthal gave the $`O(2^{n/2})`$-query upper bound. Against two *adaptive* queries
no capacity bound was known. This repository gives, to our knowledge, the first, and shows
it is essentially attained. The clean model is strictly smaller than the garbage model, and
the results here neither improve nor contradict the one-query garbage-model bounds;
Section 7 of the paper prices the distance between the two models.

## License

Text under CC BY 4.0 (`LICENSE`); code under MIT (`LICENSE-CODE`).
