/-
# Q2 — the two-query Gaussian mean width bound (deterministic core)

INFORMAL SOURCE
  `paper/paper.md`, Corollary 4.2 in its audited form (real sign selectors, the
  symmetric weighting, constant 8):

    w_2 := E_G sup_{g ∈ {±1}^Q} Re Tr(G* T_g)  ≤  8 · min(N,Q) · sqrt(log 2Q)

  for a two-query architecture T_g = B · D_g · W · D_g · C with B, W, C
  arbitrary contractions.

WHAT IS FORMALISED, AND WHAT IS NOT
  The symmetric route has four algebraic steps and two analytic inputs.  The
  algebraic steps are deterministic linear algebra: the exact trace expansion
  with a selector-independent coefficient matrix; the passage from
  Re(gᵀMg) to a Hermitian quadratic form, which uses that g is real; and the
  constancy of the weighted selector norm ‖D^{1/2}g‖² = Tr D over the whole
  sign cube, which is why the supremum is free.  The analytic inputs are the
  trace budget with the pointwise absorption, and Tropp's matrix Gaussian
  series bound (noncommutative Khintchine).  Mathlib has no matrix
  concentration inequality, and formalising Tropp's bound is an independent
  project of comparable scale to this one.

  So the analytic inputs are taken as **explicit hypotheses on the
  statements**, never as axioms.  A hypothesis is visible in the type and
  leaves `#print axioms` clean (no `sorryAx`, no custom axiom), preserving
  this repository's kernel-cleanliness discipline while keeping the imported
  inputs legible.

  Formalised here:
    * `trace_expansion`      the exact expansion, M independent of g
    * `conj_quadratic`       conjugation is transposition when g is real
    * `two_re_quadratic`     2·Re(gᵀMg) = gᵀ(M+Mᴴ)g, **because g is real**;
                             this is the step the unimodular route of
                             Theorem 4.1 avoids (paper, Remark 4.6)
    * `cube_norm_const`      ‖D^{1/2}g‖² = Tr D, CONSTANT on {±1}^Q
    * `width_of_khintchine`  assembly from the two hypotheses

  `cube_norm_const` is the conceptual heart: it is *why* the supremum over
  the 2^Q sign words costs nothing, and hence why the proof needs no chaining
  and never exchanges sup with E.

STATUS: no `sorry`, no custom axioms.
-/
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.LinearAlgebra.Matrix.ConjTranspose
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Push

namespace TwoQuery.Q2

open Matrix BigOperators Finset

variable {N Q : ℕ}

/-- The two-query architecture map `T_g = B · D_g · W · D_g · C`. -/
noncomputable def Tg (B : Matrix (Fin N) (Fin Q) ℂ) (W : Matrix (Fin Q) (Fin Q) ℂ)
    (C : Matrix (Fin Q) (Fin N) ℂ) (g : Fin Q → ℂ) : Matrix (Fin N) (Fin N) ℂ :=
  B * Matrix.diagonal g * W * Matrix.diagonal g * C

/-- The `g`-independent coefficient matrix `A = C · Gᴴ · B`. -/
noncomputable def Amat (B : Matrix (Fin N) (Fin Q) ℂ) (C : Matrix (Fin Q) (Fin N) ℂ)
    (G : Matrix (Fin N) (Fin N) ℂ) : Matrix (Fin Q) (Fin Q) ℂ :=
  C * Gᴴ * B

/-- The paper's `M`, with `M p q = W p q * (C · Gᴴ · B) q p`. -/
noncomputable def Mmat (B : Matrix (Fin N) (Fin Q) ℂ) (W : Matrix (Fin Q) (Fin Q) ℂ)
    (C : Matrix (Fin Q) (Fin N) ℂ) (G : Matrix (Fin N) (Fin N) ℂ) :
    Matrix (Fin Q) (Fin Q) ℂ :=
  fun p q => W p q * Amat B C G q p

/-- The quadratic form `gᵀ M g` for a real selector `g`. -/
noncomputable def quad (M : Matrix (Fin Q) (Fin Q) ℂ) (g : Fin Q → ℝ) : ℂ :=
  ∑ p : Fin Q, ∑ q : Fin Q, (g p : ℂ) * (g q : ℂ) * M p q

/-! ### Step 1 — the exact expansion -/

/-- **Step 1.**  `Tr(Gᴴ T_g) = ∑_{p,q} g p · g q · M p q`, with `M`
    **independent of `g`**. -/
theorem trace_expansion (B : Matrix (Fin N) (Fin Q) ℂ) (W : Matrix (Fin Q) (Fin Q) ℂ)
    (C : Matrix (Fin Q) (Fin N) ℂ) (G : Matrix (Fin N) (Fin N) ℂ) (g : Fin Q → ℂ) :
    Matrix.trace (Gᴴ * Tg B W C g)
      = ∑ p : Fin Q, ∑ q : Fin Q, g p * g q * Mmat B W C G p q := by
  classical
  have hcyc :
      Matrix.trace (Gᴴ * Tg B W C g)
        = Matrix.trace (Amat B C G * Matrix.diagonal g * W * Matrix.diagonal g) := by
    unfold Tg Amat
    have h : Gᴴ * (B * Matrix.diagonal g * W * Matrix.diagonal g * C)
        = (Gᴴ * B * Matrix.diagonal g * W * Matrix.diagonal g) * C := by
      simp [Matrix.mul_assoc]
    rw [h, Matrix.trace_mul_comm]
    simp [Matrix.mul_assoc]
  rw [hcyc]
  -- Entrywise: (A · D · W · D) p p = ∑_q (A p q · g q) · (W q p · g p).
  have hdiag : ∀ p : Fin Q,
      (Amat B C G * Matrix.diagonal g * W * Matrix.diagonal g) p p
        = ∑ q : Fin Q, (Amat B C G p q * g q) * (W q p * g p) := by
    intro p
    -- Re-associate to (A·D)·(W·D) first, then read off entries.
    rw [Matrix.mul_assoc (Amat B C G * Matrix.diagonal g) W (Matrix.diagonal g),
      Matrix.mul_apply]
    refine Finset.sum_congr rfl fun q _ => ?_
    rw [Matrix.mul_diagonal, Matrix.mul_diagonal]
  simp only [Matrix.trace, Matrix.diag_apply]
  rw [Finset.sum_congr rfl fun p (_ : p ∈ Finset.univ) => hdiag p]
  -- Both are the same double sum; swap the order and match summands.
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => ?_
  simp only [Mmat]
  ring

/-! ### Step 1' — the real part, and why `g` must be real -/

/-- Conjugating the quadratic form transposes the matrix — **this uses that
    `g` is real**. -/
theorem conj_quadratic (M : Matrix (Fin Q) (Fin Q) ℂ) (g : Fin Q → ℝ) :
    (starRingEnd ℂ) (quad M g) = quad Mᴴ g := by
  classical
  unfold quad
  -- Distribute the conjugation FIRST, so `conj_ofReal` can clear the real
  -- casts before anything normalises `starRingEnd` to `star`.
  simp only [map_sum, map_mul, Complex.conj_ofReal]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => ?_
  simp only [Matrix.conjTranspose_apply, starRingEnd_apply]
  ring

/-- **Step 1'.**  `2 · Re(gᵀ M g) = gᵀ (M + Mᴴ) g` for a **real** selector.

    This is the paper's Remark 4.6: the identity holds *because* `g` is real,
    and it is the step the symmetric route relies on.  The unimodular route of
    Theorem 4.1 never forms a Hermitian form and so does not need it; that
    route is not formalised here. -/
theorem two_re_quadratic (M : Matrix (Fin Q) (Fin Q) ℂ) (g : Fin Q → ℝ) :
    ((2 : ℝ) * (quad M g).re : ℂ) = quad (M + Mᴴ) g := by
  classical
  have hsum : quad (M + Mᴴ) g = quad M g + quad Mᴴ g := by
    unfold quad
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun q _ => ?_
    simp [Matrix.add_apply, mul_add]
  rw [hsum, ← conj_quadratic M g, Complex.add_conj]
  push_cast
  ring

/-! ### Step 2 — the cube-invariance that makes the supremum free -/

/-- **Step 2.**  With `d ≥ 0` and `v = D^{1/2} g`, the norm `‖v‖² = ∑ d p`
    is **constant over the sign cube** `{±1}^Q`: it does not depend on `g`.

    This is why the supremum over the `2^Q` sign words costs nothing, and
    hence why the proof needs no chaining and never exchanges `sup` with `E`.
    The objective is degree 2 in the selector, so the adaptive supremum is
    free. -/
theorem cube_norm_const (d : Fin Q → ℝ) (g : Fin Q → ℝ)
    (hg : ∀ p, g p = 1 ∨ g p = -1) :
    ∑ p : Fin Q, d p * (g p) ^ 2 = ∑ p : Fin Q, d p := by
  refine Finset.sum_congr rfl fun p _ => ?_
  rcases hg p with h | h <;> rw [h] <;> ring

/-! ### Step 6 — assembly -/

/-- **Step 6, the assembly.**  Given

    * `hTr`  : the trace budget `Tr D ≤ 2 · min N Q`  (the Frobenius bound
               `‖B‖_F² + ‖C‖_F² ≤ 2 min(N,Q)` for contractions), and
    * `hKh`  : the Khintchine bound `E‖K‖ ≤ 4 · sqrt (log (2Q))`  (Tropp's
               matrix Gaussian series bound, paper step 6), and
    * `hAbs` : the pointwise absorption `w ≤ Tr D · E‖K‖`  (paper step 3,
               `|v*Kv| ≤ ‖K‖ ‖v‖²` with `‖v‖²` constant on the cube by
               `cube_norm_const`),

    the width bound follows:  `w ≤ 8 · min N Q · sqrt (log (2Q))`.

    The two analytic imports are hypotheses, not axioms: they are visible in
    the type, and `#print axioms` stays clean. -/
theorem width_of_khintchine
    (w trD EK : ℝ)
    (hTr : trD ≤ 2 * (min N Q : ℝ))
    (hKh : EK ≤ 4 * Real.sqrt (Real.log (2 * Q)))
    (hEK : 0 ≤ EK)
    (hAbs : w ≤ trD * EK) :
    w ≤ 8 * (min N Q : ℝ) * Real.sqrt (Real.log (2 * Q)) := by
  have h1 : trD * EK ≤ (2 * (min N Q : ℝ)) * (4 * Real.sqrt (Real.log (2 * Q))) :=
    mul_le_mul hTr hKh hEK (by positivity)
  have h2 : (2 * (min N Q : ℝ)) * (4 * Real.sqrt (Real.log (2 * Q)))
      = 8 * (min N Q : ℝ) * Real.sqrt (Real.log (2 * Q)) := by ring
  linarith [hAbs, h1, h2.le, h2.ge]

/-- The independent unimodular selector has a constant leverage-weighted norm.
    No concentration or operator-norm inequality is imported into this identity. -/
theorem torus_norm_const {ι : Type*} [Fintype ι]
    (d : ι → ℝ) (x : ι → ℂ) (hx : ∀ p, Complex.normSq (x p) = 1) :
    ∑ p, d p * Complex.normSq (x p) = ∑ p, d p := by
  simp only [hx, mul_one]

/-- Scalar assembly for the asymmetric route. The pointwise absorption and
    Gaussian estimate remain explicit hypotheses, not formalized conclusions. -/
theorem asymmetric_width_assembly (w left right EK h : ℝ)
    (hl : 0 ≤ left) (hr : 0 ≤ right)
    (hAbs : w ≤ left * right * EK) (hKh : EK ≤ 2 * h) :
    w ≤ 2 * left * right * h := by
  have hp : 0 ≤ left * right := mul_nonneg hl hr
  calc
    w ≤ left * right * EK := hAbs
    _ ≤ left * right * (2 * h) := mul_le_mul_of_nonneg_left hKh hp
    _ = 2 * left * right * h := by ring

end TwoQuery.Q2
