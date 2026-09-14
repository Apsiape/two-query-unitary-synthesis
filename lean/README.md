# Scoped Lean formalization

A kernel-checked formalization of the deterministic part of the two-query width bound,
in its audited form (paper Corollary 4.2, real sign selectors, symmetric weighting,
constant 8).

## What is certified

Seven declarations in `TwoQuery/Q2Width.lean`, namespace `TwoQuery.Q2`:

| theorem | content |
|---|---|
| `trace_expansion` | $`\mathrm{Tr}(G^{*}T_g)=\sum_{p,q}g_pg_qM_{pq}`$ with $`M`$ independent of the selector |
| `conj_quadratic` | conjugating the quadratic form transposes its matrix when the selector is real |
| `two_re_quadratic` | $`2\,\mathrm{Re}(g^{\top}Mg)=g^{\top}(M+M^{*})g`$ for real $`g`$ |
| `cube_norm_const` | $`\|D^{1/2}g\|^2=\mathrm{Tr}D`$ for every $`g\in\lbrace\pm1\rbrace^Q`$: the weighted selector norm is constant on the cube |
| `width_of_khintchine` | assembly: trace budget + pointwise absorption + Khintchine bound ⇒ $`w\le8\min(N,Q)\sqrt{\log 2Q}`$ |
| `torus_norm_const` | independent unimodular selectors have constant leverage-weighted squared norm |
| `asymmetric_width_assembly` | scalar assembly for the asymmetric route, conditional on absorption and concentration |

`cube_norm_const` is the conceptual heart. It is why the supremum over $`2^Q`$ sign words
costs nothing, hence why the proof needs no chaining and never exchanges a supremum with
an expectation.

## What is not certified

The three analytic inputs enter `width_of_khintchine` as explicit hypotheses on the
statement, never as axioms:

- the trace budget $`\mathrm{Tr}D\le2\min(N,Q)`$ for contractions;
- the pointwise absorption $`w\le\mathrm{Tr}D\cdot\mathbb E\|K\|`$;
- Tropp's matrix Gaussian series bound $`\mathbb E\|K\|\le4\sqrt{\log 2Q}`$.

Mathlib has no matrix concentration inequality; formalizing Tropp's bound is a separate
project. A hypothesis is visible in the theorem's type and leaves `#print axioms` clean,
which is the discipline followed here: nothing unproved is hidden.

The complete unimodular-selector route of Theorem 4.1 (constant 2, complex phases,
two tables) is **not** formalized. Version 1.0.3 adds its constant-norm identity and
conditional scalar assembly, not the Schur variance or matrix-concentration proof.

## Reproduce

```
lake build
lake env lean axiom_check.lean
```

Every axiom list must be a subset of `[propext, Classical.choice, Quot.sound]`. The
recorded output is in `RECORDED-AXIOMS.txt`. Toolchain: `leanprover/lean4:v4.30.0`
(`lean-toolchain`); Mathlib `v4.30.0` (`lakefile.toml`, `lake-manifest.json`). A first
build fetches the Mathlib cache, which is several gigabytes. `.lake/` is build output and is
not tracked.
