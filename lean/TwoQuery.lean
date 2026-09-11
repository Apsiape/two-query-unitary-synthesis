import TwoQuery.Q2Width

/-!
# TwoQuery — the deterministic spine of the two-query width theorem

One module, `TwoQuery.Q2Width`, formalises the algebraic steps of Corollary 4.2 of
`paper/paper.md` in its audited symmetric form (real sign selectors, constant 8): the
exact trace expansion with a selector-independent coefficient matrix, the real-part
identity that the symmetric route relies on, the constancy of the weighted selector
norm over the sign cube (the fact that makes the supremum free), and the final
assembly.  The analytic inputs (trace budget, pointwise absorption, matrix
Khintchine) enter as explicit hypotheses on the statements, never as axioms.  See
`lean/README.md` for exactly what is and is not certified.
-/
