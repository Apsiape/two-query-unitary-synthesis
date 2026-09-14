# Garbage-Model Bounds From Clean-Model Ones

*Companion to Section 7 of the paper. A sufficient conjugation transfer from
a clean insertion-class capacity bound to a garbage-model lower bound. No
minimal-degree or optimal-transfer theorem is asserted.*

---

## 1. The three models

Fix $`N=2^n`$, an oracle register of dimension $`Q`$, and a phase oracle acting as the
diagonal sign matrix $`D_g`$, $`g\in\lbrace\pm1\rbrace^Q`$. A $`t`$-query architecture is a
fixed sequence of unitaries interleaved with $`t`$ oracle calls,

```math
V(g)=A_t\,(D_g\otimes I_k)\,A_{t-1}\cdots(D_g\otimes I_k)\,A_0
```

on $`\mathbb C^Q\otimes\mathbb C^k`$, where $`k`$ is the workspace dimension, ancillas included. Only $`g`$
depends on the target. Write $`\iota\psi=\psi\otimes|0\rangle`$ for the clean input
embedding and

```math
S_g:=V(g)\,\iota:\ \mathbb C^N\to\mathbb C^N\otimes\mathbb C^D,\qquad D=\text{junk dimension},
```

for the isometry the architecture realizes on the clean input. The three models differ only
in what is demanded of $`S_g`$.

- **Clean.** Some $`g`$ gives $`\|S_g-V_U\|_{\mathrm{op}}\le\varepsilon`$ with
  $`V_U\psi=(U\psi)\otimes|0\rangle`$. The ancilla returns to the distinguished state.
- **Garbage, keep form.** $`S_g\psi=(U_g\psi)\otimes|j_g\rangle`$ with $`|j_g\rangle`$ a unit
  vector depending on the oracle but not on the input.
- **Garbage, discard form** (the model of [LMW24] and [DLM26]). Only the induced channel
  $`\mathrm{Tr}_{\mathrm{junk}}(S_g\,\cdot\,S_g^{*})`$ must be within diamond distance
  $`\varepsilon`$ of the unitary channel of $`U`$; nothing is posited about the discarded
  register.

**The bridge.** The keep form is derived from the discard form. A unitary channel has Kraus
rank one, so continuity of the Stinespring representation yields a unitary $`W`$ on the junk
register with

```math
\bigl\|S_g\psi-(U_g\psi)\otimes|j_g\rangle\bigr\|\le c\sqrt\varepsilon
\qquad\text{for every unit }\psi,\ \ |j_g\rangle:=W^{*}|e_0\rangle ,
\tag{1.1}
```

with an absolute constant $`c`$: [KSW08] gives $`W`$ on a common enlarged environment, which
suffices because embedding the junk register isometrically in a larger one changes none of
the conjugated objects below, and [vE23] (arXiv:2308.15389) gives $`W`$ on the junk register
itself with $`c=\sqrt2`$ for the Kraus-rank-one case. Everything below is proved for the keep
form and inherited by the discard form at the cost of $`O(\sqrt\varepsilon)`$ in operator
norm.

**The insertion classes.** Clean-model capacity theorems are proved for the abstract
$`s`$-insertion class

```math
\mathcal P_s=\bigl\lbrace\,g\mapsto Y_sD_gY_{s-1}\cdots D_gY_0:\ \|Y_i\|_{\mathrm{op}}\le1\,\bigr\rbrace ,
```

arbitrary contractions in every slot, rectangular allowed. Each member is an architecture
with the contractions and address projectors fixed; its output family varies over $`g`$
alone. Capacity bounds are uniform over architectures, not bounds on the union of their
outputs. We bound each output family's packing entropy at
the normalized Frobenius scale $`\rho\sqrt N`$. The sign word acts through address projectors
of arbitrary rank, $`D_g=\sum_{j\le Q}g_jP_j`$ on a register of arbitrary dimension, as in
Theorem 6.3 of the paper, so that a bound is stated in terms of the address count $`Q`$ and
not of the register dimension. The proved $`s=2`$ case is Theorem 6.3 of the paper:

```math
\log\mathrm{Pack}_{\rho\sqrt N}\lbrace T_g\rbrace\le C_\rho\,N\log(2QN).
\tag{1.2}
```

A **clean-$`s`$ capacity theorem** means a statement of the shape
$`\log\mathrm{Pack}_{\rho\sqrt N}\le C_\rho N\phi_N(Q)`$ for each architecture in $`\mathcal P_s`$,
where $`\phi_N`$ may depend on target dimension and address count but not workspace,
with constants uniform over the fixed contractions,
address projectors and workspace dimension.
Since $`D_g`$ is linear in $`g`$ and $`g_p^2=1`$, every entry of a member of $`\mathcal P_s`$
is a multilinear polynomial of degree at most $`s`$ in the signs: a $`t`$-query architecture
produces a degree-$`t`$ object, and clean-$`s`$ tools are degree-$`s`$ tools.

---

## 2. A quadratic observable removes the junk

In the clean model the invariant one studies is obtained by contracting the ancilla against
a fixed bra: $`(I_N\otimes\langle0|)S_g=U_g`$, linear in $`S_g`$, degree $`t`$. In the garbage
model that contraction is unavailable, because $`\langle0|j_g\rangle`$ is an unknown and
possibly vanishing scalar. What survives is the quadratic object

```math
S_g^{*}(X\otimes I_D)S_g=\langle j_g|j_g\rangle\,U_g^{*}XU_g=U_g^{*}XU_g ,
\tag{2.1}
```

quadratic in $`S_g`$, degree $`2t`$ in the signs.

**Lemma 2.1 (junk-blindness; Lemma 7.1 of the paper).** *Call a functional $`L`$ on
isometries junk-blind if $`L(S)=L((I\otimes W)S)`$ for every unitary $`W`$ on the junk
register, equivalently if it depends only on the induced channel. If $`L`$ is junk-blind and
linear, then $`L\equiv0`$.*

*Proof.* Take $`W=e^{i\theta}I`$: junk-blindness gives $`L(S)=L(e^{i\theta}S)=e^{i\theta}L(S)`$
for every $`\theta`$. For a real-linear $`L`$ take $`\theta=\pi`$. $`\square`$

The lemma rules out a nonzero junk-blind linear functional of the packet itself.
Conjugation provides a useful quadratic invariant with an explicit $`2t`$-insertion
representation. This proves sufficiency of the transfer below, not necessity of
oracle degree $`2t`$: substitution of an oracle circuit into an invariant can cause
degree cancellations. For example, $`D_gD_g=I`$ is a two-query packet with constant
conjugated observables. No lower bound on the degree of every possible invariant,
or on every width method, is asserted.

---

## 3. The conjugation transfer

**Theorem 3.1 (Theorem 7.2 of the paper).** *A clean-$`2t`$ capacity theorem with rate
$`\phi_N`$ implies that no universal garbage-$`t`$ synthesizer exists with $`\phi_N(Q)=o(N)`$. The
transfer depends on no workspace or junk dimension, changes no oracle, and tolerates a
constant diamond error. At $`t=2`$ the target is clean-4.*

*Proof.* Four steps.

**(a) The circuit is already a packet.** Restricting $`V(g)`$ to the clean input gives,
verbatim, $`S_g=A_t(D_g\otimes I)A_{t-1}\cdots(D_g\otimes I)A_0\iota\in\mathcal P_t`$ with
rectangular endpoint slots.

**(b) Conjugate.** Fix a contraction $`X`$ on $`\mathbb C^N`$ and set
$`T_g:=S_g^{*}(X\otimes I_D)S_g`$. Because $`D_g\otimes I`$ is Hermitian, expanding
$`S_g^{*}`$ by reversal gives

```math
T_g=\underbrace{\iota^{*}A_0^{*}}_{Y_{2t}}(D_g\otimes I)A_1^{*}\cdots(D_g\otimes I)\underbrace{A_t^{*}(X\otimes I_D)A_t}_{Y_t\ \text{(middle slot)}}(D_g\otimes I)\cdots A_1(D_g\otimes I)\underbrace{A_0\iota}_{Y_0},
\tag{3.1}
```

exactly a member of $`\mathcal P_{2t}`$ on the register $`\mathbb C^Q\otimes\mathbb C^k`$ with
address projectors $`P_j\otimes I_k`$: the same sign word inserted $`2t`$ times, every slot a
contraction. Nothing about the architecture is used beyond unitarity of the $`A_i`$. The
product is palindromic, $`T_g=Y_g^{*}\,A_t^{*}(X\otimes I)A_t\,Y_g`$ with
$`Y_g=(D_g\otimes I)A_{t-1}\cdots(D_g\otimes I)A_0\iota`$.
`verify_conjugation_transfer.py` checks (3.1) as a four-insertion product at $`t=2`$.

**(c) Junk cancels identically.** In the keep form, by (2.1),
$`T_g=\langle j_g|j_g\rangle U_g^{*}XU_g=U_g^{*}XU_g`$: no reference vector, no alignment, no
dependence on $`D`$. In the discard form, (1.1) gives
$`\|T_g-U_g^{*}XU_g\|_F\le2\sqrt N\|X\|\,c\sqrt\varepsilon`$, a constant budget at the
normalized Frobenius scale. More directly, channel duality gives
$`\|T_g-U_g^{*}XU_g\|_{\mathrm{op}}\le\varepsilon\|X\|`$ at unhalved
diamond error $`\varepsilon`$, so this transfer does not require (1.1).

**(d) The conjugation orbit has full packing entropy.** Take $`X=R`$ a balanced reflection
($`R=R^{*}`$, $`R^2=I`$, $`\mathrm{Tr}R=0`$). Then $`U\mapsto U^{*}RU`$ maps onto the balanced
reflections, a copy of the Grassmannian $`\mathrm{Gr}(N/2,N)`$ of real dimension $`N^2/2`$,
every image point having Frobenius norm $`\sqrt N`$. The explicit concentration-and-covering
argument in paper Section 7 gives, for
$`\rho`$ below an absolute constant, a $`\rho\sqrt N`$-separated set of $`e^{c_\rho N^2}`$
points, each of the form $`U_i^{*}RU_i`$. A universal garbage-$`t`$ synthesizer must
synthesize $`U_1,\dots,U_m`$, so the output family $`\lbrace T_g\rbrace`$ of a member of $`\mathcal P_{2t}`$
contains an $`e^{c_\rho N^2}`$-point separated code, while a clean-$`2t`$ capacity theorem
caps the same quantity by $`C_\rho N\phi_N(Q)=o(N^2)`$ when $`\phi_N(Q)=o(N)`$. $`\square`$

**Remarks.**

- The separated code is chosen *on the orbit* and then lifted. One may not push an
  arbitrary Haar-scale unitary code forward through $`U\mapsto U^{*}XU`$, which collapses
  the stabilizer of $`X`$.
- No circuit is appended and no second oracle is queried; the transfer re-reads the
  architecture's own data. It is therefore valid in fixed-oracle and random-oracle models
  alike.
- **Calibration at $`t=1`$.** There $`2t=2`$ and the hypothesis is Theorem 6.3 of the paper,
  so the transfer delivers an unrestricted garbage-model one-query bound at no overhead
  (Corollary 7.3 of the paper). This is consistent with the fact that the one-query bound
  of [LMW24] is stated in the discard model with clean-model tools: at $`t=1`$
  degree-doubled still means degree 2.

---

## 4. What a clean-4 capacity theorem means, and why it is open

The object required is a packing bound for $`\mathcal P_4`$ at the normalized Frobenius
scale: a quartic concentration statement. Nothing about unitary targets enters; the members
$`T_g`$ of (3.1) are Hermitian contractions.

The obstruction is concrete. In the degree-2 proof the sign supremum is absorbed because
the leverage-weighted selector vectors have selector-independent norms: every slot touches
an endpoint, so the whole coefficient mass is priced on the endpoint budget
$`\|B\|_F^2+\|C\|_F^2`$. By the palindromic form it suffices to bound the sandwich subclass
$`Y_g^{*}XY_g`$ with $`Y_g`$ a two-insertion packet whose left endpoint is the identity of
the full register, and there the endpoint budget is the register dimension itself: the two
inner sign layers are not priced by anything of order $`N`$. What is needed is a domination
certificate that reprices the inner layers on the target budget alone. This is open, and
the paper claims nothing about it.

---

## 5. Scope of the transfer

A clean four-insertion capacity bound suffices for a two-query garbage-model
lower bound. The linear-invariant lemma does not establish that four is necessary,
even within a broadly described family of width methods. We therefore do not claim
optimality of this route or an equivalence between the two problems.

Other garbage-removal compilers require their own oracle, accuracy and workspace
accounting. No unproved comparison with such a compiler is needed here.

---

## 6. The restricted route, and why it is not included

Breaking junk-blindness by a fixed unit reference $`|r\rangle`$ gives the exact degree-$`t`$
packet $`(I_N\otimes\langle r|)S_g=\langle r|j_g\rangle U_g`$, with the junk entering only as the
scalar $`\alpha_g=\langle r|j_g\rangle`$. A pigeonhole over a net of the junk sphere then
yields a two-query garbage-model bound valid when $`\min(Q,4D)=o(N^2)`$, i.e. when the
junk register has fewer than about $`2n`$ qubits. That statement rests on the paper's
Theorem 4.1 plus the net argument. It exists in the author's working records but has not
passed an independent audit, and it is **not** claimed in this repository. The clean-model
results of the paper do not depend on it.
