import TwoQuery

/-!
Prints the axiom dependencies of every certified theorem.  Each line must show exactly
`[propext, Classical.choice, Quot.sound]` (or a subset): no `sorryAx`, no custom axiom.

    lake env lean axiom_check.lean
-/

#print axioms TwoQuery.Q2.trace_expansion
#print axioms TwoQuery.Q2.conj_quadratic
#print axioms TwoQuery.Q2.two_re_quadratic
#print axioms TwoQuery.Q2.cube_norm_const
#print axioms TwoQuery.Q2.width_of_khintchine
