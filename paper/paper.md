# The Clean Two-Query Threshold for Unitary Synthesis

**Seth Douglas**

*Version 1.0.1 — 2026-09-11*

> Every numbered claim below that admits a finite check is exercised by a script in `verify/`,
> and `STATUS.md` records, claim by claim, which checks exist, what has been independently
> audited and what has not.
> Display and inline mathematics follow GitHub's renderer conventions.

---

## Abstract

Aaronson and Kuperberg asked whether every $`n`$-qubit unitary can be implemented by a
polynomial-size quantum circuit making few queries to a classical oracle. At one query the
answer is known to be no. Between one query and the $`O(2^{n/2})`$ queries that suffice, no
capacity bound was known against adaptive queries.

We settle the two-query case in the clean model, where the workspace must be returned to
its initial state, by pinning both of its resource thresholds up to a logarithm in the
error. Write $`N = 2^n`$.

- **Width.** For any fixed two-query architecture $`T_{x,y} = B\,D_x\,W\,D_y\,C`$ with
  arbitrary contractions $`B, W, C`$ and independent unimodular phase selectors $`x, y`$, the
  Gaussian mean width of the realizable family is at most
  $`2\min(N,Q)\sqrt{\log(2Q)}`$, where $`Q`$ is the dimension of the register the oracle acts
  on, ancillas included. The supremum over selectors is absorbed pointwise; no chaining is
  used.
- **Lower bounds.** Consequently a universal clean two-query circuit needs
  $`\Omega(2^n)`$ qubits, and, by an address-compression lemma, an oracle input length of
  $`\Omega(2^n)`$ bits regardless of workspace.
- **Matching upper bound.** A permutation-transport construction implements every unitary
  cleanly with two queries, $`O(2^n\log(1/\varepsilon))`$ qubits and
  $`O(2^n\log(1/\varepsilon))`$ address bits, at error $`\varepsilon`$.

So the clean two-query qubit threshold and address-length threshold are both
$`\Theta(2^n)`$ up to the factor $`\log(1/\varepsilon)`$. We also show that the bound is
tight at workspace $`\Theta(N^2)`$, explain why the one-query argument does not extend, and
prove a transfer theorem reducing the garbage-model two-query problem to a clean
four-insertion width bound, which remains open; no junk-blind width argument can do with
fewer insertions.

---

## 1. Introduction

### 1.1 The problem

The unitary synthesis problem of Aaronson and Kuperberg [AK07] asks whether every
$`n`$-qubit unitary can be implemented by a polynomial-time quantum algorithm with few
queries to a classical oracle. The known landmarks are:

- **[AK07, Thm 6.7]** one query, exact model: at most $`4^N`$ distinct unitaries, at most
  $`2N`$ oracle bits used;
- **[LMW24]** one query, or polynomially many parallel (non-adaptive) queries, in the
  approximate/discard model, oracle input length $`o(2^n)`$;
- **[DLM26]** explicit one-query separations, and Question 1.4 posing the general-depth
  problem;
- **[BY26]** an oracle input-length bound $`(2-o(1))\log d`$ valid at any query count;
- **[Ros26]** an upper bound: $`O(2^{n/2})`$ clean queries suffice, with $`\mathrm{poly}(n)`$
  ancillas.

Against two or more *adaptive* queries no capacity bound was known; on the other side there
were only the trivial lookup and the permutation synthesizer of [DLM26]. This paper gives,
to our knowledge, the first capacity bound at two adaptive queries and shows it is
essentially tight.

### 1.2 The model, and why "clean" is the right restriction

A fixed circuit with two oracle calls realizes a target $`U`$ if some oracle makes it
implement $`U`$ on the target register and, in the **clean** model, return every auxiliary
register to a fixed state. Nothing is required of the circuit for other oracles; in
particular, unlike the setting of [AK07, Thm 6.7], whose footnote 13 notes that the theorem
needs the circuit to implement some unitary for every oracle, we impose no such requirement,
and our lower bounds hold without it. The only known one-query universal synthesizer, the
Bernstein–Vazirani lookup attributed to Yuen in [LMW24], is not clean: its workspace ends
holding the description string. So the class we bound is strictly smaller than the discard
model of [LMW24], and **our conclusions are correspondingly stronger on a correspondingly
smaller class; they neither improve nor contradict [LMW24] or [BY26].** Section 7 prices the
distance between the two models.

### 1.3 What is new

[LMW24, Def. B.1] already name the quantity we bound: an oracle circuit is
$`S`$-small if it realizes at most $`S`$ distinct unitaries, and their Appendix B proves that
small capacity implies a lower bound. They then observe that counting arguments *"become
ineffective once the adversary has even a small number of ancilla qubits."* The implication
is theirs; the missing premise is that a two-query architecture's capacity *is* small. This
paper supplies that premise, in a metric-entropy form that degrades only as the square root
of the total qubit count: an adversary needs $`\Omega(2^n)`$ ancilla qubits to neutralize it.

Two rounds, not two queries. A phase selector $`D_x`$ is an arbitrary diagonal on the whole
register, and a product of phase oracles applied in one round is again such a diagonal.
One selector factor therefore subsumes arbitrarily many parallel queries, so Theorem 4.1
constrains *two rounds of adaptivity with arbitrarily many parallel queries per round*.

The obstruction [LMW24, §2.5] name for going past one query, bounding the operator norm of
$`2^M`$ matrices simultaneously, is met by absorbing the entire selector supremum
*pointwise* into a single operator norm (Section 4).

Beyond the width bound, this version adds: the extension from real sign words to
independent unimodular selectors, which covers diagonal-unitary oracles and two different
tables; an address-length lower bound that holds with any amount of workspace; the
matching upper bound; and the garbage-model transfer.

**Indexing.** [DLM26, Question 1.4] ask whether synthesizing $`F_tH\cdots F_1`$ requires $`t`$
sequential queries; the case $`t`$ is a lower bound against $`t-1`$ queries, so their proved
case $`t=2`$ is a one-query bound and a two-query bound has the strength of their case
$`t=3`$. Our result has that strength for Haar-scale targets; we say nothing about their
explicit family.

---

## 2. Setting

**Definition 2.1 (two-query architecture).** A triple $`(B,W,C)`$ with $`C`$ a $`Q\times N`$
contraction, $`W`$ a $`Q\times Q`$ contraction and $`B`$ an $`N\times Q`$ contraction, realizing

```math
T_{x,y} \;=\; B\,D_x\,W\,D_y\,C,
\qquad D_x=\mathrm{diag}(x),\; D_y=\mathrm{diag}(y),
\qquad x,y\in\mathbb T^Q,
```

where $`\mathbb T=\lbrace z\in\mathbb C:|z|=1\rbrace`$. The **repeated real sign word** is the
special case $`x=y=g\in\lbrace\pm1\rbrace^Q`$, written $`T_g`$; it is the shape of a two-query
circuit with one $`\pm1`$ phase oracle. Independent $`x,y`$ cover diagonal-unitary oracles
and circuits that query two different tables.

Ancillas are inside $`Q`$: an architecture on $`n`$ target qubits and $`a`$ ancillas has
$`Q=2^{n+a}`$, and we write $`w=\log_2Q`$ for the total qubit count. The bounds below are
stated in $`\log(2Q)`$, i.e. in $`w`$; Theorem 6.3 then separates the oracle input length
from the workspace.

**Notation.** $`b_p=Be_p\in\mathbb C^N`$ is the $`p`$-th column of $`B`$ and
$`c_q=e_q^{\top}C\in\mathbb C^N`$ the $`q`$-th row of $`C`$, with **leverages**
$`\lambda_p=\|b_p\|^2`$, $`\mu_q=\|c_q\|^2`$, $`\Lambda=\mathrm{diag}(\lambda)`$,
$`M_\mu=\mathrm{diag}(\mu)`$, and $`P_\lambda, P_\mu`$ the projections onto their supports.
Since $`B`$ and $`C`$ are contractions,
$`\sum_p\lambda_p=\|B\|_F^2\le\min(N,Q)`$ and $`\sum_q\mu_q=\|C\|_F^2\le\min(N,Q)`$.

**Proposition 2.2 (unitarity criterion).** Write $`\Psi(x,y)=D_xWD_yC`$ and
$`\Pi=I_Q-B^{*}B`$. Then

```math
T_{x,y}^{*}T_{x,y}
\;=\; I_N \;-\;(I_N-C^{*}C)\;-\; C^{*}D_y^{*}\,(I_Q-W^{*}W)\,D_yC \;-\; \Psi(x,y)^{*}\,\Pi\,\Psi(x,y).
```

*All three subtracted terms are positive semidefinite, so $`T_{x,y}`$ is unitary exactly when
$`C`$ is an isometry, $`(I_Q-W^{*}W)D_yC=0`$ and $`\Pi\Psi(x,y)=0`$.* When $`C`$ and $`W`$ are
isometries the criterion reduces to $`\Pi\Psi=0`$.

*Proof.* $`T^{*}T=\Psi^{*}B^{*}B\Psi=\Psi^{*}\Psi-\Psi^{*}\Pi\Psi`$ and
$`\Psi^{*}\Psi=C^{*}D_y^{*}W^{*}WD_yC=C^{*}C-C^{*}D_y^{*}(I-W^{*}W)D_yC`$, while
$`C^{*}C=I_N-(I_N-C^{*}C)`$. Each subtracted term is positive semidefinite because
$`I-C^{*}C`$, $`I-W^{*}W`$ and $`\Pi`$ are, for contractions; and a sum of positive
semidefinite operators vanishes only when each does. $`\square`$

> An earlier draft stated the criterion with only the last term, which is correct only for
> isometric $`C`$ and $`W`$; with $`C=0`$, a legal contraction, the shorter criterion is
> satisfied while $`T=0`$. `verify_complex_phase.py` (e) checks the three-term identity for
> contractive $`C`$ and $`W`$ and exhibits the failure of the shorter forms.

**Definition 2.3 (capacity).** The capacity at scale $`\rho`$ is
$`\log\mathrm{Pack}_{\rho\sqrt N}\lbrace T_{x,y}\rbrace`$, the logarithm of the largest
family of realizable operators that is pairwise $`\rho\sqrt N`$-separated in Frobenius norm.
This is the physically meaningful scale for unitary targets, since $`\|U\|_F=\sqrt N`$.
**The capacity claim is asserted only in this metric**; a raw, un-normalized reading is
false by a factor $`N`$ (Corollary 5.1).

---

## 3. Why the one-query argument does not extend

[AK07, Thm 6.7] writes $`U(X_{ij})=I+E_i+E_j`$ and applies [AK07, Lem. 6.6] to the
anticommutation relations $`E_iE_j^{\dagger}+E_jE_i^{\dagger}=0`$. This is pure affineness and
is unavailable at degree two. One might still hope that a combinatorial invariant of the
**lawful set**, the set of selector words for which the circuit's output is unitary, controls
capacity. Two facts speak against building on that hope, one at depth two and one beyond it.

**Proposition 3.1 (at depth two the final layer is invisible).** *For a two-query
architecture, write $`\Psi(x,y)=D_xWD_yC`$ and let $`\mathcal L`$ be its lawful set. Then the
number of distinct realizable targets $`\lbrace T_{x,y}:(x,y)\in\mathcal L\rbrace`$ equals the
number of distinct $`\Psi(x,y)`$, $`(x,y)\in\mathcal L`$. Given $`\Psi`$ and $`\mathcal L`$, the
count does not depend on $`B`$.*

*Proof.* For lawful $`(x,y)`$, $`T=B\Psi`$ is unitary, so by Proposition 2.2 $`\Psi`$ is an
isometry and $`B`$ is isometric on $`\mathrm{ran}\,\Psi(x,y)`$. For a contraction $`B`$ the set of
vectors on which it is isometric is the subspace $`\ker(I-B^{*}B)`$, so $`B`$ is injective on
the sum of these ranges over the lawful set, and $`B\Psi(x,y)=B\Psi(x',y')`$ forces
$`\Psi(x,y)=\Psi(x',y')`$. $`\square`$

**Proposition 3.2 (beyond depth two, pre-final data does not determine the count).** *In the
three-insertion class $`\Psi_g=D_gW_2D_gW_1D_gC`$, $`T_g=B\Psi_g`$, with $`Q=8`$ and a
one-dimensional target ($`N=1`$), there are two architectures sharing $`C`$, $`W_1`$ and
$`B`$ and differing only in $`W_2`$, with identical lawful sets (a linear code of $`16`$ sign
words closed under pointwise product) and identical pre-final data (the frames
$`C,\,D_gC,\,D_gW_1D_gC`$ on the lawful set and their Gram kernels), whose numbers of distinct
realizable targets are $`2`$ and $`16`$.*

Verified exhaustively in `verify_decisive_pair.py`: exactness residuals below
$`10^{-14}`$; the first invariant to differ is the last frame, $`\|K_3^A-K_3^B\|=2.00`$. The
same script checks Proposition 3.1 on two-query architectures with non-empty lawful sets
under several final layers.

**Remark 3.3 (what this motivates, and nothing more).** Neither fact is a bound. The first
says that at depth two capacity is a property of the pre-final map restricted to the lawful
set, an object whose Boolean structure we have no tool to count; the second says that beyond
depth two even that reduction fails, so any bound computed from pre-final data alone is at
best a maximum over completions. Both point away from combinatorial invariants of the lawful
set and toward a quantity that sees the operators $`T_{x,y}`$ themselves: their metric
entropy, which Section 4 controls through a Gaussian mean width. A $`2`$-versus-$`16`$ gap
at $`Q=8`$ establishes nothing asymptotically.

---

## 4. The width theorem

> **Terminology.** "Width" means the Gaussian mean width of a set of matrices, in the sense
> of convex geometry. [LMW24, Lem. 4.11] use "width" for an unrelated leverage-normalized
> overlap statistic.

Let $`G`$ be an $`N\times N`$ matrix with independent standard complex Gaussian entries,
$`\mathbb E|G_{ij}|^2=1`$.

**Theorem 4.1 (two-query Gaussian mean width, unimodular selectors).** *For every two-query
architecture $`(B,W,C)`$ of Definition 2.1,*

```math
w_2 \;:=\; \mathbb E_G\,\sup_{x,y\in\mathbb T^Q}\bigl|\mathrm{Tr}(G^{*}T_{x,y})\bigr|
\;\le\; \|B\|_F\,\|C\|_F\;\cdot\;2\sqrt{\log(2Q)}
\;\le\; 2\min(N,Q)\sqrt{\log(2Q)} .
```

*The supremum sits inside the expectation; no exchange of supremum and expectation occurs.*

*Proof.* Six steps; the first three are exact identities.

1. **Exact expansion.** Since $`D_x=\sum_p x_pe_pe_p^{\top}`$ and likewise for $`D_y`$,

   ```math
   \mathrm{Tr}(G^{*}T_{x,y})
   \;=\;\sum_{p,q}x_p\,y_q\,W_{pq}\,\bigl(c_q\,G^{*}\,b_p\bigr)
   \;=\; x^{\top}M\,y,
   \qquad M_{pq}:=W_{pq}\,(c_qG^{*}b_p).
   ```

   **The matrix $`M`$ does not depend on $`x`$ or $`y`$.** If $`\lambda_p=0`$ then
   $`b_p=0`$ and row $`p`$ of $`M`$ vanishes; if $`\mu_q=0`$ then column $`q`$ vanishes.
   Restricting to the supports is therefore exact.

2. **Leverage normalization.** Put $`K:=\Lambda^{-1/2}M\,M_\mu^{-1/2}`$ on the supports and
   $`u=\Lambda^{1/2}x`$, $`v=M_\mu^{1/2}y`$, so that $`x^{\top}My=u^{\top}Kv`$.

3. **Pointwise absorption.** For unimodular selectors the norms
   $`\|u\|^2=\sum_p\lambda_p|x_p|^2=\|B\|_F^2`$ and $`\|v\|^2=\|C\|_F^2`$ are **constant over
   the torus**. Hence, pointwise in $`G`$,

   ```math
   \sup_{x,y\in\mathbb T^Q}\bigl|x^{\top}My\bigr|
   \;\le\;\|B\|_F\,\|C\|_F\;\|K\|_{\mathrm{op}} .
   ```

   The supremum over the entire selector torus has been absorbed into one operator norm:
   no net, no chaining, no union bound. This is the heart of the proof.

   > The inequality of this step is not new. It is the diagonal-scaling relaxation of
   > quadratic maximization over the cube and the torus that runs from Delorme–Poljak and
   > Poljak–Rendl through Nesterov, Alon–Naor and Charikar–Wirth; the degree-$`t`$ structure
   > is the polynomial method [BBCMW98], of which [AK07, Thm 6.7] is the $`t=1`$ case. What is
   > ours is the composition: polynomial method, then this relaxation, then matrix
   > concentration, then Sudakov, run against a synthesis architecture, and the variance
   > computation of the next step that makes the composition close.

4. **$`K`$ is a matrix Gaussian series.** Writing $`G=\sum_{ij}G_{ij}e_ie_j^{\top}`$,

   ```math
   K=\sum_{i,j}\overline{G_{ij}}\,A_{ij},
   \qquad
   (A_{ij})_{pq}=\frac{W_{pq}\,(b_p)_i\,(c_q)_j}{\sqrt{\lambda_p\mu_q}} .
   ```

   Its two variance operators have closed forms. With
   $`\Gamma_B:=\Lambda^{-1/2}B^{*}B\Lambda^{-1/2}`$ and
   $`\Gamma_C:=M_\mu^{-1/2}CC^{*}M_\mu^{-1/2}`$, both positive semidefinite with unit diagonal
   on their supports,

   ```math
   \sum_{ij}A_{ij}A_{ij}^{*}=\overline{\Gamma_B}\circ\bigl(W P_\mu W^{*}\bigr),
   \qquad
   \sum_{ij}A_{ij}^{*}A_{ij}=\overline{\Gamma_C}\circ\bigl(W^{*}P_\lambda W\bigr),
   ```

   where $`\circ`$ is the entrywise product. Both closed forms are checked against Monte
   Carlo in `verify_complex_phase.py` (c).

5. **Both variance operators are dominated by the identity.** If $`\Gamma\succeq0`$ has
   $`\Gamma_{pp}\le1`$ then $`\|\Gamma\circ X\|_{\mathrm{op}}\le\|X\|_{\mathrm{op}}`$ for every
   $`X`$: writing $`\Gamma_{pp'}=u_p^{*}u_{p'}`$ with $`\|u_p\|\le1`$, the matrix
   $`\Gamma\circ X`$ is the compression of $`I\otimes X`$ by the contraction
   $`e_p\mapsto u_p\otimes e_p`$; $`\overline{\Gamma}`$ is positive semidefinite with the same
   diagonal, so the lemma applies to it as well. Since $`WP_\mu W^{*}`$ and $`W^{*}P_\lambda W`$ are positive
   semidefinite of norm at most one, both variance operators are positive semidefinite of
   norm at most one, and they are supported on $`P_\lambda`$ and $`P_\mu`$ respectively.
   Only $`\|W\|\le1`$ is used; the middle layer is otherwise arbitrary.

6. **Matrix Khintchine.** For a real Gaussian series $`Z=\sum_k\gamma_kA_k`$ with
   $`d_1\times d_2`$ coefficients, Tropp's bound gives
   $`\mathbb E\|Z\|\le\sqrt{2v\log(d_1+d_2)}`$ with $`v=\max\lbrace\|\sum A_kA_k^{*}\|,\|\sum A_k^{*}A_k\|\rbrace`$
   [Tro15, Thm 4.1.1]. Splitting each complex Gaussian into real and imaginary parts of
   variance $`1/2`$ gives two real series with variance parameter $`v/2\le1/2`$, hence
   $`\mathbb E\|K\|\le2\sqrt{\log(2Q)}`$ by the triangle inequality. (Treating the two
   parts as one real series of $`2N^2`$ terms gives the sharper $`\sqrt2`$; the constant 2
   is what the verifier tests and all that is needed.) Assemble with step 3. $`\square`$

**Corollary 4.2 (real sign words; the certified form).** *For the repeated real sign word,*

```math
\mathbb E_G\sup_{g\in\lbrace\pm1\rbrace^Q}\mathrm{Re}\,\mathrm{Tr}(G^{*}T_g)
\;\le\;2\min(N,Q)\sqrt{\log(2Q)} .
```

*Proof.* Take $`x=y=g`$ and use $`|\mathrm{Re}\,z|\le|z|`$. $`\square`$

The form of this corollary that was independently audited (see `STATUS.md`) carries the
constant $`8`$ and is proved by a symmetric weighting: with $`d_p=\lambda_p+\mu_p`$ and
$`D=\mathrm{diag}(d)`$ one writes $`\mathrm{Re}(g^{\top}Mg)=v^{*}Kv`$ with
$`v=D^{1/2}g`$, $`\|v\|^2=\mathrm{Tr}D\le2\min(N,Q)`$, and bounds $`\mathbb E\|K\|\le4\sqrt{\log2Q}`$.
The asymmetric route above was also verified in that audit, as the $`\sqrt{LR}`$
strengthening; both routes are valid and incomparable instance by instance. The
deterministic spine of the symmetric route is formalized in `lean/`.

**Remark 4.3 (the weighting is sufficient, not unique).** Under $`d\mapsto td`$ the
variance operators scale by $`t^{-2}`$, so any $`d'\succeq\lambda+\mu`$ works subject to the
trace budget. Balance is necessary: over random contractive architectures, $`d=\lambda`$
alone or $`d=\mu`$ alone fails the domination on a substantial fraction of instances, while
the balanced or the asymmetric normalization succeeds on all of them.

**Remark 4.4 (zero leverage is lossless).** $`\lambda_p=0`$ kills row $`p`$ of $`M`$ and
$`\mu_q=0`$ kills column $`q`$, so the supports may differ and $`K`$ is rectangular in
general. Nothing downstream changes; the concentration step stays at dimension at most
$`2Q`$. The same holds when the middle layer is itself rectangular, $`Q_1\times Q_2`$, with a
selector of length $`Q_1`$ on the left and $`Q_2`$ on the right: nothing in steps 1 to 6 uses
$`Q_1=Q_2`$, the endpoint norms are as before, and the Khintchine dimension becomes
$`Q_1+Q_2`$.

**Remark 4.5 (provenance of the logarithm).** $`\sqrt{\log(2Q)}`$ is the standard Khintchine
logarithm, $`2Q`$ being Tropp's $`d_1+d_2`$. It is the same factor appearing as
$`\sqrt{2\ln2M}`$ in [LMW24, Lem. 4.11]. Inherited, not discovered.

**Remark 4.6 (the selector need not be real).** An earlier draft proved only the real case
and disclaimed any extension to diagonal-unitary oracles, because its symmetric route
rewrote $`\mathrm{Re}(g^{\top}Mg)`$ as a Hermitian form, which uses $`g\in\mathbb R^Q`$.
The asymmetric route never does this: step 3 uses only that $`|x_p|=|y_q|=1`$. The real
restriction is gone, and with it the restriction to a single table.

---

## 5. From width to capacity, and the constant

**Corollary 5.1 (capacity).** *For fixed $`\rho>0`$ and $`Q\ge N`$,*

```math
\log\mathrm{Pack}_{\rho\sqrt N}\lbrace T_{x,y}\rbrace
\;\le\;\frac{2C_S^2\,w_2^2}{\rho^2N}
\;\le\;\frac{K}{\rho^2}\,N\log(2Q),
\qquad K:=8C_S^2,
```

*where $`C_S`$ is the constant in Sudakov's minoration for the process
$`X_T=\mathrm{Re}\,\mathrm{Tr}(G^{*}T)`$, whose intrinsic distance is $`\|S-T\|_F/\sqrt2`$.*

*Proof.* Sudakov's minoration [Ver18, Thm 7.4.1]: a family that is $`\epsilon`$-separated in
the process metric and lies in the index set of a Gaussian process with
$`\mathbb E\sup X_T=w`$ has $`\epsilon\sqrt{\log\mathrm{Pack}}\le C_Sw`$. At Frobenius
separation $`\rho\sqrt N`$ the process separation is $`\rho\sqrt{N/2}`$, and
$`\mathbb E\sup\mathrm{Re}\,\mathrm{Tr}(G^{*}T)\le w_2`$; insert Theorem 4.1. The hypothesis
$`Q\ge N`$ always holds, since the target register is part of the register the oracle acts
on. $`\square`$

> **Erratum, stated where a reader will quote it.** Corollary 5.1 is **false in raw
> Frobenius**, by a factor $`N`$. It is true only in the normalized metric of Definition 2.3.
> This is a corollary, not a link in the proof of Theorem 4.1, but it is the statement
> downstream consumers quote. Re-check the separation metric before instantiating it.

**Corollary 5.2 (qubit lower bound).** *For every $`0<\rho<\sqrt2`$,
$`\log\mathrm{Pack}_{\rho\sqrt N}U(N)\ge c\,(1-\rho^2/2)^2N^2`$ for an absolute $`c>0`$. Hence no
fixed clean two-query architecture realizes every unitary, exactly or with operator-norm
error below $`\rho/4`$, unless $`\log Q=\Omega_\rho(N)`$. In qubit count:*

```math
\boxed{\;w\;=\;\Omega(2^n)\;}
```

*for a universal clean two-query architecture on $`n`$ target qubits.*

*Proof.* For fixed $`U`$ and Haar-random $`V`$, $`\|U-V\|_F^2=2N-2\,\mathrm{Re}\,\mathrm{Tr}(U^{*}V)`$,
so $`V`$ lies in the Frobenius ball of radius $`\rho\sqrt N`$ about $`U`$ exactly when
$`\mathrm{Re}\,\mathrm{Tr}(U^{*}V)\ge N(1-\rho^2/2)`$. The function
$`V\mapsto\mathrm{Re}\,\mathrm{Tr}(U^{*}V)`$ has Haar mean zero and is $`\sqrt N`$-Lipschitz in the
Frobenius metric, and concentration of Lipschitz functions on $`U(N)`$ [Mec19, Ch. 5]
gives $`\Pr(f\ge\mathbb Ef+t)\le\exp(-c'Nt^2/L^2)`$; at $`t=N(1-\rho^2/2)`$, $`L=\sqrt N`$, the
ball has Haar measure at most $`\exp(-c'(1-\rho^2/2)^2N^2)`$. A maximal
$`\rho\sqrt N`$-separated set is a $`\rho\sqrt N`$-covering, so its size is at least the
reciprocal of that measure. If an architecture realizes each member $`U_i`$ of such a set to
operator-norm error below $`\rho/4`$, the realized operators are within Frobenius distance
$`\rho\sqrt N/4`$ of the $`U_i`$ and hence $`\rho\sqrt N/2`$-separated; Corollary 5.1 at
scale $`\rho/2`$ gives $`c(1-\rho^2/2)^2N^2\le4K\rho^{-2}N\log(2Q)`$, i.e.
$`\log Q=\Omega_\rho(N)`$. $`\square`$

For comparison, the trivial clean two-query lookup, which stores every entry of $`U`$ in the
table, uses $`\tilde\Theta(N^2)`$ qubits, so at this point the threshold sits in
$`[\Omega(N),\tilde O(N^2)]`$; Theorem 6.4 closes the gap.

**Proposition 5.3 (the constant matters).** *The permutation family of Section 6 forces
$`K\ge\tfrac14`$ asymptotically, the maximum being approached at separation $`\rho^2=1`$.*
Computed from Gilbert–Varshamov packings of permutation codes
(`verify_crho_requirement.py`), with $`d=N`$:

| $`d`$ | 16 | 256 | 4096 | 65536 | $`\to\infty`$ |
|---|---|---|---|---|---|
| required $`K`$ | 0.126 | 0.159 | 0.185 | 0.200 | **0.25** |

Against this, Sudakov's minoration holds with $`C_S\le2.09`$: comparing the process on the
packing set with independent $`\mathcal N(0,\epsilon^2/2)`$ variables through the
Sudakov–Fernique inequality [Ver18, Thm 7.2.11], and using that the expected maximum of
$`m\ge2`$ independent standard normals is at least $`0.677\sqrt{\log m}`$ (the minimum over
$`m`$ is at $`m=2`$, where the ratio is $`1/\sqrt{\pi\log2}`$), gives
$`\epsilon\sqrt{\log\mathrm{Pack}}\le2.09\,\mathbb E\sup X_T`$. Hence $`K=8C_S^2\le35`$ for
the constant-2 route of Theorem 4.1, and $`K\le560`$ for the audited constant-8 form. Either
way the requirement $`K\ge1/4`$ is met with two orders of magnitude to spare. We record
the requirement because it is real: at $`K<1/4`$ the theorem would be contradicted by a
published construction.

---

## 6. Tightness, address length, and the matching upper bound

### 6.1 Tight at workspace of order N squared

**Proposition 6.1.** *The clean two-query permutation-synthesis algorithm of [DLM26] lies
inside Definition 2.1 and attains the ceiling of Corollary 5.1 up to the constant factor
$`4K`$.*

[DLM26] synthesize any permutation unitary $`P_\pi`$ on $`d=N`$ dimensions in two clean
queries, querying $`\pi`$ then $`\pi^{-1}`$. Index addresses by $`(b,x,y)`$ with
$`b\in\lbrace0,1\rbrace`$ and $`x,y\in[d]`$, so $`Q=2d^2`$:

```math
g_{0,x,y}=(-1)^{y\cdot\pi(x)},\qquad g_{1,x,y}=(-1)^{x\cdot\pi^{-1}(y)},
```
```math
C:|x\rangle\mapsto d^{-1/2}\sum_y|0,x,y\rangle,\qquad
W=X_b\otimes H_X\otimes H_Y,\qquad
B:|b,x,y\rangle\mapsto\delta_{b,1}\,d^{-1/2}|y\rangle .
```

$`W`$ is fixed and independent of $`\pi`$; a single sign word serves both queries because the
state occupies block $`0`$ at the first query and block $`1`$ at the second. Tracing through,
$`T_g=P_\pi`$ exactly. Verified in `verify_permutation_family.py` over all $`24`$
permutations at $`d=4`$ and $`200`$ random permutations at $`d=8`$ with residuals below
$`10^{-14}`$. Since $`\|P_\pi-P_\sigma\|_F^2=2\,\mathrm{Ham}(\pi,\sigma)`$, permutation
codes at Hamming distance $`\alpha d`$ give $`\rho^2=2\alpha`$, and Gilbert–Varshamov gives
$`\log\mathrm{Pack}=\Theta(d\log d)=\Theta(N\log2Q)`$.

**Corollary 6.2 (tight in both factors at $`Q=\Theta(N^2)`$).** *Running Sudakov backwards
on the permutation family gives $`w_2\ge c\rho\sqrt d\sqrt{d\log d}=\Theta(N\sqrt{\log2Q})`$ at
$`Q=2N^2`$. So at that scale Theorem 4.1 is attained up to a constant, and neither the factor
$`N`$ nor the factor $`\sqrt{\log2Q}`$ can be improved by more than a constant there.* Padding
an architecture with unused workspace raises $`\log Q`$ without changing $`w_2`$, so this is a
statement about the scale $`Q=\Theta(N^2)`$, not about every polynomial. At exponential $`Q`$
the logarithm is necessary as well: by Theorem 6.4, at $`\log Q=\Theta(N\log(1/\varepsilon))`$
the class $`\varepsilon`$-covers $`U(N)`$ in operator norm, and since
$`|\mathrm{Tr}(G^{*}(T-U))|\le\|G\|_{S_1}\|T-U\|_{\mathrm{op}}`$ this gives
$`w_2\ge(1-\varepsilon)\,\mathbb E\|G\|_{S_1}=\Theta(N^{3/2})`$ there. No bound of the form
$`N\cdot\mathrm{polylog}(N)`$ can therefore hold uniformly in $`Q`$.

### 6.2 The oracle's input length alone

Theorem 4.1 bounds the total register dimension. The next theorem removes the workspace
from the accounting entirely. In this section and the next, $`M`$ denotes a number of
addresses; the coefficient matrix $`M`$ of Section 4 does not appear.

**Theorem 6.3 (address-length lower bound).** *Let the oracle act on the $`Q`$-dimensional
register through $`M`$ addresses of arbitrary rank,*

```math
D_x=\sum_{j=1}^{M}x_j\,P_j,\qquad x\in\mathbb T^M,
```

*with $`P_1,\dots,P_M`$ mutually orthogonal projections summing to $`I_Q`$ (for a Boolean
table on $`\ell`$ input bits, $`M=2^\ell`$ and $`P_j`$ projects onto address $`j`$ tensored
with the whole workspace). Then for every two-query architecture*

```math
\mathbb E_G\sup_{x,y\in\mathbb T^M}\bigl|\mathrm{Tr}(G^{*}B D_xWD_yC)\bigr|
\;\le\;2\,\|B\|_F\|C\|_F\sqrt{\log(2MN)}
\;\le\;2N\sqrt{\log(2MN)} ,
```

*and consequently universal clean two-query synthesis requires $`\log M=\Omega(N)`$, i.e.
an oracle input length of $`\Omega(2^n)`$ bits, regardless of the number of ancillas.*

*Proof.* Let $`L:=\sum_j\mathrm{ran}(P_jB^{*})`$ and $`R:=\sum_j\mathrm{ran}(P_jC)`$, subspaces of
$`\mathbb C^Q`$ of dimension at most $`MN`$ each. Because the $`P_j`$ are mutually
orthogonal, $`L`$ is the orthogonal direct sum of the $`\mathrm{ran}(P_jB^{*})`$ and is
invariant under every $`P_j`$; likewise $`R`$. Let $`V_L,V_R`$ be isometries onto $`L,R`$ and
$`\Pi_L=V_LV_L^{*}`$, $`\Pi_R=V_RV_R^{*}`$. Since $`(BP_j)^{*}=P_jB^{*}`$ has range in $`L`$,
$`BP_j=BP_j\Pi_L`$, and summing over $`j`$ gives $`B=B\Pi_L`$; similarly $`C=\Pi_RC`$. Hence

```math
BD_xWD_yC
\;=\;(BV_L)\Bigl(\sum_jx_j\,V_L^{*}P_jV_L\Bigr)\bigl(V_L^{*}WV_R\bigr)\Bigl(\sum_ky_k\,V_R^{*}P_kV_R\Bigr)(V_R^{*}C).
```

In an orthonormal basis of $`L`$ adapted to its decomposition, $`\sum_jx_jV_L^{*}P_jV_L`$ is
diagonal with each $`x_j`$ repeated $`\mathrm{rank}(P_j|_L)`$ times; it is a unimodular
diagonal, and likewise on $`R`$. The compressed middle $`V_L^{*}WV_R`$ is a contraction, and
$`\|BV_L\|_F^2=\mathrm{Tr}(B^{*}B\Pi_L)=\|B\|_F^2`$, $`\|V_R^{*}C\|_F=\|C\|_F`$. Theorem 4.1
applies to the compressed architecture: its statement has a square middle layer, but its
proof never uses squareness (Remark 4.4), and with a $`\dim L\times\dim R`$ middle the
Khintchine step gives $`\sqrt{\log(\dim L+\dim R)}\le\sqrt{\log(2MN)}`$. The supremum over
all unimodular diagonals on $`L`$ and $`R`$ dominates the supremum over the repeated
patterns. Corollary 5.1 then reads $`\log\mathrm{Pack}\le K\rho^{-2}N\log(2MN)`$, and the
$`\Omega(N^2)`$ demand of a Haar-scale code forces $`\log M=\Omega(N)`$. $`\square`$

For a bit-flip oracle $`O_f|x,b\rangle=|x,b\oplus f(x)\rangle`$ the eigen-decomposition of the
response bit gives $`M+1`$ addresses; nothing changes. The two selectors may also use
different address decompositions, as for two tables of lengths $`M_1`$ and $`M_2`$: the proof
is unchanged with $`\dim L\le M_1N`$ and $`\dim R\le M_2N`$. The factorization, the dimension
count, the preserved norms and the contractivity are checked on eighty random architectures
with address projectors of arbitrary rank, and with two different decompositions, in
`verify_address_compression.py`.

### 6.3 The matching upper bound

**Theorem 6.4 (permutation transport).** *For every $`\varepsilon\in(0,\tfrac1{10})`$ there
are a target-independent encoder and decoder such that every $`U\in U(N)`$ is implemented
by a clean two-query circuit with*

| resource | cost |
|---|---|
| Boolean tables | one |
| Boolean queries | two |
| oracle input length | $`O(N\log(1/\varepsilon))`$ bits |
| workspace | $`O(N\log(1/\varepsilon))`$ qubits |
| ordinary gates | $`\mathrm{poly}(N,1/\varepsilon)`$ |
| half-diamond error, all reference-entangled inputs | at most $`\varepsilon`$ |
| target-dependent advice, postselection | none |

*The construction restores its workspace to within the stated error, so it is clean; the
realized map on the target register is within operator-norm distance $`3\varepsilon/4`$ of
$`U`$; and the circuit is of the repeated-sign-word form $`T_g`$ of Definition 2.1, one
$`\pm1`$ table serving both queries.*

*Construction, condensed.* The full statement with all estimates is
`notes/permutation-transport.md`; the ingredients are checked at small $`N`$ in
`verify_koopman_transport.py`.

Encode a state as a function on a finite grid rather than as a list of columns. Quantize the
standard normal distribution into $`K`$ equiprobable intervals with conditional means
$`q_a`$, put $`s^2=\tfrac1K\sum_aq_a^2`$ and $`e^2=1-s^2`$; a direct estimate gives
$`e^2\le32K^{-1/3}`$ (natural logarithms throughout), so $`K=\mathrm{poly}(1/\varepsilon)`$ makes
$`e`$ as small as needed, independently of $`N`$. For addresses $`a\in[K]^{2N}`$, $`M=K^{2N}`$, define
$`w_a\in\mathbb C^N`$ by $`(w_a)_k=(q_{a_k}+iq_{a_{N+k}})/(s\sqrt2)`$ and the fixed encoding
$`(E\psi)_a=w_a^{*}\psi/\sqrt M`$. Independence, zero means and the normalization give
$`E^{*}E=I_N`$ exactly.

For a target $`U`$ with realification $`O_U`$, draw $`X\sim\mathcal N(0,I_{2N})`$ and quantize
coordinatewise: $`A=Q(X)`$, $`B=Q(O_UX)`$. Both are uniform on $`[K]^{2N}`$, so
$`D_{ba}=M\Pr(A=a,B=b)`$ is doubly stochastic and, by Birkhoff–von Neumann, a convex
combination $`\sum_{\ell=1}^{m}p_\ell P_{\pi_\ell}`$ of at most $`m\le M^2`$ permutation matrices. The quantization
residual $`r(X)=X-q_{Q(X)}`$ has covariance $`e^2I_{2N}`$, and
$`q_{Q(O_UX)}-O_Uq_{Q(X)}=O_Ur(X)-r(O_UX)`$, whence, using
$`(a-b)(a-b)^{\top}\preceq2aa^{\top}+2bb^{\top}`$,

```math
\sum_\ell p_\ell\,\Delta_{\pi_\ell}^{*}\Delta_{\pi_\ell}\;\preceq\;\frac{4e^2}{s^2}\,I_N,
\qquad \Delta_\pi:=P_\pi E-EU .
```

This is a uniform operator bound, valid for every input. To make the mixture a single
permutation, put $`\delta:=\varepsilon/2`$, take $`K`$ with $`e^2\le\delta^2/100`$, and add a
uniform coherent label register over $`R\ge100M^2/\delta^2`$ copies, a power of two fixed by
$`M`$ and $`\delta`$ alone and hence independent of the target, with integer multiplicities
$`k_\ell\approx Rp_\ell`$, $`\sum_\ell|p_\ell-k_\ell/R|\le2m/R`$; set $`\Pi_U(j,a)=(j,\pi_{\ell(j)}(a))`$
with the enlarged, still target-independent, encoding $`\tilde E\psi=|+_R\rangle\otimes E\psi`$.
Then

```math
\bigl\|\Pi_U\tilde E-\tilde EU\bigr\|_{\mathrm{op}}^2
=\Bigl\|\sum_\ell\tfrac{k_\ell}{R}\Delta_{\pi_\ell}^{*}\Delta_{\pi_\ell}\Bigr\|_{\mathrm{op}}
\;\le\;\frac{4e^2}{s^2}+\frac{8m}{R}\;<\;\delta^2 ,
```

so the label register returns to $`|+_R\rangle`$ up to $`\delta`$, uncorrelated with the
input. The permutation $`\Pi_U`$ and its inverse are stored in one Boolean table on
$`1+2t`$ input bits, $`t=\log_2(MR)=O(N\log(1/\varepsilon))`$, and are applied by two queries
in the standard way (load $`\Pi_U(z)`$ into a blank register, swap, unload with the inverse
table). The fixed encoder is an explicit polynomial-size circuit and its inverse decodes.
Total isometry error at most $`3\varepsilon/4`$; two isometries at operator distance
$`\alpha`$ induce channels at half-diamond distance at most $`\alpha`$. $`\square`$

The lower bound of Theorem 6.3 shows that for this mechanism, and indeed for any monomial
action between fixed encodings, the exponential address length is unavoidable; the note
records that converse and its one loophole, a target-dependent final ancilla state, which
does not affect the clean model treated here.

**Corollary 6.5 (the thresholds).** *For clean two-query synthesis of all $`n`$-qubit
unitaries at operator-norm error $`\varepsilon`$, $`\varepsilon`$ a small constant, both the qubit
count $`w^{\star}`$ and the oracle input length $`\ell^{\star}`$ satisfy*

```math
\Omega(2^n)\;\le\;w^{\star},\;\ell^{\star}\;\le\;O\bigl(2^n\log(1/\varepsilon)\bigr).
```

---

## 7. The garbage model: what the two-query problem costs there

In the garbage model the workspace may end in a target-dependent state
$`|j_g\rangle`$, independent of the input, or may simply be discarded. [LMW24] and [DLM26]
state their one-query bounds there. This section prices, exactly, what a two-query bound in
that model would require of the clean machinery.

**Setting.** A $`t`$-query circuit
$`V(g)=A_t(D_g\otimes I_K)A_{t-1}\cdots(D_g\otimes I_K)A_0`$ on $`\mathbb C^Q\otimes\mathbb C^K`$,
where $`Q`$ is the number of addresses and $`K`$ the workspace dimension, restricted to the
clean input, $`S_g:=V(g)\iota`$ with $`\iota\psi=\psi\otimes|0\rangle`$, is an isometry
$`\mathbb C^N\to\mathbb C^Q\otimes\mathbb C^K\cong\mathbb C^N\otimes\mathbb C^D`$, with
$`D=QK/N`$ the junk dimension, whose entries are multilinear polynomials of degree $`t`$ in
the signs. In the
**keep form** $`S_g\psi=(U_g\psi)\otimes|j_g\rangle`$. The **discard form** only requires the
induced channel to be within diamond distance $`\varepsilon`$ of $`U_g`$; since a unitary
channel has Kraus rank one, continuity of the Stinespring representation recovers the keep
form up to $`O(\sqrt\varepsilon)`$ in operator norm: [KSW08] gives the required unitary on a
common enlarged environment, which suffices here because embedding the junk register
isometrically in a larger one changes none of the conjugated objects below, and [vE23] gives
it on the junk register itself with constant $`\sqrt2`$. So everything below transfers to the
discard form at constant error.

**Lemma 7.1 (junk-blindness).** *Call a functional $`L`$ on isometries junk-blind if
$`L(S)=L((I\otimes W)S)`$ for every unitary $`W`$ on the junk register, equivalently if it
depends only on the induced channel. If $`L`$ is junk-blind and linear then $`L\equiv0`$.*

*Proof.* Take $`W=e^{i\theta}I`$: junk-blindness gives $`L(S)=L(e^{i\theta}S)=e^{i\theta}L(S)`$
for all $`\theta`$; for a real-linear $`L`$ take $`\theta=\pi`$. $`\square`$

The width-and-Sudakov machinery of Sections 4 and 5 prices a Gaussian process that is
*linear* in the realized operator. By the lemma, no junk-blind linear functional of a
garbage packet is nonzero; the invariants that factor through the channel are quadratic, the
basic one being $`S_g^{*}(X\otimes I)S_g=U_g^{*}XU_g`$, of degree $`2t`$ in the signs. This is
the exact sense in which cleanliness halves the degree, and it says what any junk-blind
width argument for garbage-$`t`$ must be: a clean bound at $`2t`$ insertions.

**Theorem 7.2 (conjugation transfer).** *Consider the $`2t`$-insertion class
$`\lbrace Y_{2t}D_gY_{2t-1}\cdots D_gY_0:\|Y_i\|\le1\rbrace`$ with $`D_g=\sum_{j\le Q}g_jP_j`$ for
mutually orthogonal projectors $`P_j`$ of arbitrary rank summing to the identity of a
register of arbitrary dimension, as in Theorem 6.3. Suppose it satisfies a capacity bound
$`\log\mathrm{Pack}_{\rho\sqrt N}\le C_\rho N\,\phi(Q)`$ for some function $`\phi`$ of the
address count. Then no universal garbage-$`t`$ synthesizer exists with $`\phi(Q)=o(N)`$, for
every workspace and junk dimension, with no change to the oracle, at constant diamond
error. At $`t=1`$, Theorem 6.3 gives $`\phi(Q)=\log(2QN)`$.*

*Proof.* Fix a balanced reflection $`R`$ ($`R=R^{*}`$, $`R^2=I`$, $`\mathrm{Tr}R=0`$) and
set $`T_g:=S_g^{*}(R\otimes I)S_g`$. Because $`D_g\otimes I`$ is Hermitian, expanding
$`S_g^{*}`$ by reversal writes $`T_g`$ exactly as a $`2t`$-insertion product of the same sign
word on $`\mathbb C^Q\otimes\mathbb C^K`$, with address projectors $`P_j\otimes I_K`$ and
contractions in every slot (the middle slot is $`A_t^{*}(R\otimes I)A_t`$; the product is
palindromic, $`T_g=Y_g^{*}\,A_t^{*}(R\otimes I)A_t\,Y_g`$ with
$`Y_g=(D_g\otimes I)A_{t-1}\cdots(D_g\otimes I)A_0\iota`$). In the
keep form $`T_g=\langle j_g|j_g\rangle\,U_g^{*}RU_g=U_g^{*}RU_g`$: the junk cancels by the
normalization of a unit vector, with no reference vector, no alignment and no dependence
on the junk dimension. The map $`U\mapsto U^{*}RU`$ lands in the balanced reflections, a
copy of the Grassmannian $`\mathrm{Gr}(N/2,N)`$ of real dimension $`N^2/2`$ on which every
point has Frobenius norm $`\sqrt N`$; volume counting there gives, for $`\rho`$ below an absolute constant, a
$`\rho\sqrt N`$-separated set of $`e^{c_\rho N^2}`$ points, each of the form $`U_i^{*}RU_i`$. A
universal synthesizer must in particular synthesize $`U_1,\dots,U_m`$, so the single family
$`\lbrace T_g\rbrace`$ contains an $`e^{c_\rho N^2}`$-point separated code, contradicting the
assumed capacity bound when $`\phi(Q)=o(N)`$. The discard form costs an additional $`O(\sqrt\varepsilon)`$ in
normalized Frobenius distance, absorbed by halving $`\rho`$. $`\square`$

The reversal identity, the junk cancellation and the orbit properties are checked in
`verify_conjugation_transfer.py`.

**Corollary 7.3 (one query, garbage model).** *At $`t=1`$ the hypothesis of Theorem 7.2 is
Theorem 6.3, applied to two insertions of a repeated real sign word. Hence no one-query
circuit, with garbage allowed and with any number of ancillas, synthesizes every unitary at
constant diamond error unless its oracle input length is $`\Omega(2^n)`$.* This recovers, by a
different route and at constant error only, the address-length form of the one-query bound
of [LMW24]; we claim no improvement on their statement, which also rules out vanishing
advantage.

**Remark 7.4 (what is not claimed).** At $`t=2`$ the hypothesis of Theorem 7.2 is a clean
**four**-insertion width bound, which is open. By the palindromic form it suffices to bound
the sandwich subclass $`Y_g^{*}XY_g`$ with $`Y_g`$ a two-insertion packet whose left endpoint
is the identity of the full register: exactly the case in which the endpoint budget that
makes step 3 of Theorem 4.1 work is unavailable, since that identity has squared Frobenius
norm equal to the register dimension. So the unrestricted garbage-model two-query problem is open, and this
paper does not claim it. A restricted two-query garbage bound, valid when the junk register
has fewer than $`2n-O(1)`$ qubits, exists in the author's working records but has not passed
an independent audit and is not included.

---

## 8. Scope, and depth three and beyond

**What Theorem 4.1 covers.** Exact and approximate clean two-query synthesis; arbitrary
contractive $`B,W,C`$; real sign oracles, diagonal-unitary oracles and two different tables;
ancillas explicitly inside $`\log(2Q)`$, and, via Theorem 6.3, the oracle input length alone.

**What it does not cover.**

1. **Three or more queries.** Nothing here bears on $`t\ge3`$. The natural conjecture, that the
   width of the $`k`$-query clean class is $`O(N\cdot k\cdot\mathrm{polylog}\,Q)`$, is pinned
   in exponent by [Ros26]: every unitary is implementable at $`k=\tilde O(\sqrt N)`$ with
   $`\log Q=\mathrm{poly}(n)`$, so the width must reach the $`\Theta(N^{3/2})`$ width of the
   unitary group there, which forces at least linear growth in $`k`$. A proof of the conjecture would
   resolve the Aaronson–Kuperberg question negatively. The two-query proof does not extend
   as it stands: with a third insertion the coefficient matrix depends on the selector.
2. **Parallel versus adaptive.** [LMW24] rule out polynomially many *parallel* queries. The
   difficulty at $`t\ge3`$ is therefore entirely in interleaving.
3. **One query.** [AK07, Thm 6.7] is stronger at $`t=1`$ in the exact model: $`4^N`$ with
   address width $`2N`$, independent of $`Q`$. The address-width concept is theirs.
4. **A remark on framing.** Quantum query algorithms are exactly the completely bounded
   multilinear forms, by the Christensen–Sinclair factorization [ABP18]. In that language
   a two-query architecture is a completely bounded bilinear map from the selector torus to
   $`N\times N`$ matrices, Theorem 4.1 bounds the Gaussian mean width of its range uniformly
   over such maps, and the $`k`$-query conjecture concerns the ranges of $`k`$-linear maps. We
   record this as a reformulation, not a result.

**Consistency checks.** The bound does not forbid the trivial clean two-query lookup
($`w=\Theta(d^2)`$), is attained up to $`\log(1/\varepsilon)`$ by Theorem 6.4, and does not
contradict [Ros26], whose algorithm uses $`t=O(\sqrt d)`$ queries. Any extension of this
method yielding $`t=\omega(\sqrt N)`$ in the clean model would contradict [Ros26] and be in
error.

---

## 9. Reproducibility

`verify/run_all.py` runs eight deterministic, numpy-only scripts, each exiting nonzero on any
failed check: Theorem 4.1 in both forms (exact expansion, pointwise absorption over the whole
cube and the torus, the closed-form variance operators against Monte Carlo and their
domination), Proposition 2.2, Propositions 3.1 and 3.2, Proposition 5.3, Proposition 6.1, the
factorization of Theorem 6.3, the ingredients of Theorem 6.4 at small $`N`$, and the
identities of Section 7. `lean/` holds a kernel-checked formalization of the deterministic
spine of Corollary 4.2, with the analytic inputs as explicit hypotheses; see `lean/README.md`
for exactly what is and is not certified.

---

## References

- **[AK07]** S. Aaronson, G. Kuperberg. *Quantum versus classical proofs and advice.*
  Theory of Computing 3(7):129–157, 2007. arXiv:quant-ph/0604056.
- **[LMW24]** A. Lombardi, F. Ma, J. Wright. *A one-query lower bound for unitary synthesis
  and breaking quantum cryptography.* STOC 2024, 979–990. arXiv:2310.08870.
- **[DLM26]** F. Dong, A. Lombardi, F. Ma. *Explicit separations for one-query unitary
  synthesis.* arXiv:2607.26478, 2026.
- **[BY26]** Z. Brakerski, H. Yuen. *On scalable pseudorandom unitaries and the unitary
  synthesis problem.* arXiv:2605.09957, 2026.
- **[Ros26]** G. Rosenthal. *Query and depth upper bounds for quantum unitaries via Grover
  search.* Quantum, 2026, article 2144. arXiv:2111.07992.
- **[ABP18]** S. Arunachalam, J. Briët, C. Palazuelos. *Quantum query algorithms are
  completely bounded forms.* ITCS 2018; SIAM J. Comput. 48(3):903–925, 2019. arXiv:1711.07285.
- **[KSW08]** D. Kretschmann, D. Schlingemann, R. F. Werner. *The information-disturbance
  tradeoff and the continuity of Stinespring's representation.* IEEE Trans. Inf. Theory
  54(4):1708–1717, 2008.
- **[vE23]** F. vom Ende. *Progress on the Kretschmann–Schlingemann–Werner conjecture.*
  arXiv:2308.15389, 2023.
- **[Tro15]** J. A. Tropp. *An introduction to matrix concentration inequalities.*
  Foundations and Trends in Machine Learning 8(1–2):1–230, 2015.
- **[Tal14]** M. Talagrand. *Upper and lower bounds for stochastic processes.* Ergebnisse der
  Mathematik 60, Springer, 2014.
- **[Ver18]** R. Vershynin. *High-Dimensional Probability.* Cambridge University Press, 2018.
  (Sudakov–Fernique, Thm 7.2.11; Sudakov minoration, Thm 7.4.1.)
- **[Mec19]** E. S. Meckes. *The Random Matrix Theory of the Classical Compact Groups.*
  Cambridge University Press, 2019. (Concentration of Lipschitz functions on $`U(N)`$, Ch. 5.)
- **[BBCMW98]** R. Beals, H. Buhrman, R. Cleve, M. Mosca, R. de Wolf. *Quantum lower bounds
  by polynomials.* FOCS 1998; J. ACM 48(4), 2001.
- **[DP93]** C. Delorme, S. Poljak. *Laplacian eigenvalues and the maximum cut problem.*
  Math. Programming 62, 1993.
- **[PR95]** S. Poljak, F. Rendl. *Nonpolyhedral relaxations of graph-bisection problems.*
  SIAM J. Optimization 5, 1995.
- **[Nes98]** Y. Nesterov. *Semidefinite relaxation and nonconvex quadratic optimization.*
  Optimization Methods and Software 9, 1998.
- **[AN06]** N. Alon, A. Naor. *Approximating the cut-norm via Grothendieck's inequality.*
  SIAM J. Computing 35(4), 2006.
- **[CW04]** M. Charikar, A. Wirth. *Maximizing quadratic programs: extending Grothendieck's
  inequality.* FOCS 2004.
- **[Bir46]** G. Birkhoff. *Tres observaciones sobre el algebra lineal.* Univ. Nac. Tucumán
  Rev. Ser. A 5, 1946.

---

## Acknowledgements and disclosure

Parts of this work were developed with AI assistance for exploration, drafting and
verification. Every proof was re-derived by the author and is machine-checked by the
included scripts; the audit status of each claim is recorded in `STATUS.md`.

## License

Text CC BY 4.0; code MIT. See `LICENSE` and `LICENSE-CODE`.
