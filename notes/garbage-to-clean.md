# Garbage-Model Bounds From Clean-Model Ones

*Companion to Section 7 of `paper/paper.md`. Why an unrestricted garbage-model
$`t`$-query bound is a clean $`2t`$-insertion problem, why nothing below degree $`2t`$ can
work for the family of methods used here, and exactly what is and is not claimed at
$`t=2`$. Identities are checked in `verify/verify_conjugation_transfer.py`; see `STATUS.md`
for audit status.*

---

## 1. The three models

Fix $`N=2^n`$, an oracle register of dimension $`Q`$, and a phase oracle acting as the
diagonal sign matrix $`D_g`$, $`g\in\lbrace\pm1\rbrace^Q`$. A $`t`$-query architecture is a
fixed sequence of unitaries interleaved with $`t`$ oracle calls,

```math
V(g)=A_t\,(D_g\otimes I_K)\,A_{t-1}\cdots(D_g\otimes I_K)\,A_0
```

on $`\mathbb C^Q\otimes\mathbb C^K`$, where $`K`$ collects workspace and ancillas. Only $`g`$
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

arbitrary contractions in every slot, rectangular allowed, and bound its packing entropy at
the normalized Frobenius scale $`\rho\sqrt N`$. The sign word acts through address projectors
of arbitrary rank, $`D_g=\sum_{j\le Q}g_jP_j`$ on a register of arbitrary dimension, as in
Theorem 6.3 of the paper, so that a bound is stated in terms of the address count $`Q`$ and
not of the register dimension. The proved $`s=2`$ case is Theorem 6.3 of the paper:

```math
\log\mathrm{Pack}_{\rho\sqrt N}\lbrace T_g\rbrace\le C_\rho\,N\log(2QN).
\tag{1.2}
```

A **clean-$`s`$ capacity theorem** means a statement of the shape
$`\log\mathrm{Pack}_{\rho\sqrt N}\le C_\rho N\phi(Q)`$ for $`\mathcal P_s`$, with $`\phi`$ a
function of the address count.
Since $`D_g`$ is linear in $`g`$ and $`g_p^2=1`$, every entry of a member of $`\mathcal P_s`$
is a multilinear polynomial of degree at most $`s`$ in the signs: a $`t`$-query architecture
produces a degree-$`t`$ object, and clean-$`s`$ tools are degree-$`s`$ tools.

---

## 2. Cleanliness halves the degree

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

**Consequence.** The width-and-Sudakov machinery of the paper prices a Gaussian process
*linear* in the realized operator. By the lemma no junk-blind linear functional of a garbage
packet is nonzero; the invariants that factor through the channel are quadratic, the basic
one being (2.1), of degree $`2t`$. Hence no
junk-blind width argument can attack garbage-$`t`$ below degree $`2t`$. In particular there
is no route from the degree-2 theorem straight to an unrestricted garbage-2 bound by any
argument of this type. The only way to stay at degree $`t`$ is to break junk-blindness by
fixing a reference vector in the junk register, which is what Section 5 does, at a price.

---

## 3. The conjugation transfer

**Theorem 3.1 (Theorem 7.2 of the paper).** *A clean-$`2t`$ capacity theorem with rate
$`\phi`$ implies that no universal garbage-$`t`$ synthesizer exists with $`\phi(Q)=o(N)`$. The
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

exactly a member of $`\mathcal P_{2t}`$ on the register $`\mathbb C^Q\otimes\mathbb C^K`$ with
address projectors $`P_j\otimes I_K`$: the same sign word inserted $`2t`$ times, every slot a
contraction. Nothing about the architecture is used beyond unitarity of the $`A_i`$. The
product is palindromic, $`T_g=Y_g^{*}\,A_t^{*}(X\otimes I)A_t\,Y_g`$ with
$`Y_g=(D_g\otimes I)A_{t-1}\cdots(D_g\otimes I)A_0\iota`$.
`verify_conjugation_transfer.py` checks (3.1) as a four-insertion product at $`t=2`$.

**(c) Junk cancels identically.** In the keep form, by (2.1),
$`T_g=\langle j_g|j_g\rangle U_g^{*}XU_g=U_g^{*}XU_g`$: no reference vector, no alignment, no
dependence on $`D`$. In the discard form, (1.1) gives
$`\|T_g-U_g^{*}XU_g\|_F\le2\sqrt N\|X\|\,c\sqrt\varepsilon`$, a constant budget at the
normalized Frobenius scale, absorbed by halving $`\rho`$.

**(d) The conjugation orbit has full packing entropy.** Take $`X=R`$ a balanced reflection
($`R=R^{*}`$, $`R^2=I`$, $`\mathrm{Tr}R=0`$). Then $`U\mapsto U^{*}RU`$ maps onto the balanced
reflections, a copy of the Grassmannian $`\mathrm{Gr}(N/2,N)`$ of real dimension $`N^2/2`$,
every image point having Frobenius norm $`\sqrt N`$. Volume counting there gives, for
$`\rho`$ below an absolute constant, a $`\rho\sqrt N`$-separated set of $`e^{c_\rho N^2}`$
points, each of the form $`U_i^{*}RU_i`$. A universal garbage-$`t`$ synthesizer must
synthesize $`U_1,\dots,U_m`$, so the single family $`\lbrace T_g\rbrace\subset\mathcal P_{2t}`$
contains an $`e^{c_\rho N^2}`$-point separated code, while a clean-$`2t`$ capacity theorem
caps the same quantity by $`C_\rho N\phi(Q)=o(N^2)`$ when $`\phi(Q)=o(N)`$. $`\square`$

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

## 5. Exactness of the clean-4 location

"Strength" is operational: what must be proved to obtain the impossibility.

- **Clean-4 suffices** (Theorem 3.1 at $`t=2`$): a capacity theorem for $`\mathcal P_4`$
  yields the unrestricted garbage-2 impossibility, for every junk dimension, with no oracle
  augmentation, at constant diamond error.
- **Nothing below degree 4 can suffice, for this family of methods** (Lemma 2.1): every
  junk-blind invariant of a garbage-2 family is at least quadratic in the packet, so no
  junk-blind linear-process argument operates below degree 4.

**The honest boundary.** The second statement concerns a method class and is proved as
such; it is not an unconditional converse. Unconditionally, only the trivial containment
holds: a clean algorithm is a garbage algorithm with $`|j_g\rangle=|0\rangle`$, so a garbage-2
impossibility implies the clean-2 one. Unrestricted garbage-2 therefore sits between clean-2
and clean-4 in strength, and the junk-blindness lemma closes the interval at the clean-4
end for every method used here. It is not claimed that a garbage-2 bound would imply a
clean-4 bound.

**Comparison with erasure.** A different transfer, which appends a circuit querying a
modified oracle, turns garbage-$`t`$ into clean-$`(t+4)`$. Best transfer depth is therefore
$`\min(2t,t+4)`$: conjugation wins for $`t\le3`$, the two tie at $`t=4`$, and erasure wins
beyond. Only conjugation is used in this repository.

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
