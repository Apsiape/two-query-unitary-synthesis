# Permutation Transport: the Matching Upper Bound in Full

*Companion to Theorem 6.4 of `paper/paper.md`. This note gives the complete construction
with every estimate, the converse for fixed encodings, and the exact scope of what the
converse leaves open. The small-case ingredients are checked in
`verify/verify_koopman_transport.py`; see `STATUS.md` for the audit status.*

---

## 0. Statement

Put $`N=2^n`$. For every $`0<\varepsilon<1/10`$ there are a target-independent encoder and
decoder such that every $`U\in U(N)`$ is implemented with the resources below.

| resource | this construction |
|---|---|
| Boolean tables | one |
| Boolean queries | two |
| complete oracle input length | $`O(N\log(1/\varepsilon))`$ bits |
| workspace | $`O(N\log(1/\varepsilon))`$ qubits |
| ordinary gates | $`\mathrm{poly}(N,1/\varepsilon)`$; at fixed error, $`N\,\mathrm{poly}(n)`$ |
| half-diamond error, all reference-entangled inputs | at most $`\varepsilon`$ |
| target-dependent advice, postselection | none |

The construction is clean up to the stated error: every auxiliary register is returned to
a fixed, target-independent state.

Moreover (Section 6), any universal implementation by a permutation with arbitrary phases
between two *fixed* encodings, with a fixed final ancilla state, needs at least
$`(N-1)\log_2(1/\varepsilon)`$ encoded qubits. So the optimal width of this mechanism is
$`\Theta(N\log(1/\varepsilon))`$.

---

## 1. The idea: encode states as functions, not as columns

Let $`z\in\mathbb C^N`$ be a standard complex Gaussian vector, $`\mathbb E\,zz^{*}=I_N`$.
Encode a state as the linear function

```math
\psi\longmapsto f_\psi(z)=z^{*}\psi .
```

This is an isometry, since $`\mathbb E|f_\psi(z)|^2=\|\psi\|^2`$. A classical change of
coordinates $`z\mapsto Uz`$ acts on functions by composition, $`(T_Uf)(z)=f(U^{*}z)`$, and

```math
(T_Uf_\psi)(z)=(U^{*}z)^{*}\psi=z^{*}U\psi=f_{U\psi}(z),
\qquad\text{i.e.}\qquad T_UE=EU .
```

The unknown amplitudes are never inspected; the machine relabels the points at which the
encoded linear function is evaluated. This is the standard observation that a change of
variables induces a linear action on functions (the Koopman representation), and the
representation itself is not claimed as new. What is new is making it finite without
assuming a quantum implementation of the rotation, and doing so with a permutation that
restores its own auxiliary register.

---

## 2. Quantizing the Gaussian

Choose a power of two $`J`$. Partition the real standard normal distribution into $`J`$
intervals of equal probability. Let $`\theta(X)\in[J]`$ name the interval containing $`X`$ and
put (the verifier's code writes `K` for $`J`$, `Q` for $`\theta`$, `A, B` for the labels
$`\alpha,\beta`$ below and `D` for the coupling $`\Xi`$)

```math
q_a=\mathbb E[X\mid \theta(X)=a],\qquad
\frac1J\sum_aq_a=0,\qquad
s^2:=\frac1J\sum_aq_a^2,\qquad
e^2:=\mathbb E\bigl(X-q_{\theta(X)}\bigr)^2=1-s^2 .
```

**Lemma 2.1 (quantizer error).** *For $`J\ge4`$, $`e^2\le32J^{-1/3}`$ (natural logarithms
throughout).*

*Proof.* Write $`g(u)=\Phi^{-1}(u)`$ for the normal quantile function and $`g_T`$ for its
clipping to $`[-T,T]`$. The clipped function is Lipschitz with constant
$`L_T=\sqrt{2\pi}\,e^{T^2/2}`$. Conditional means minimize squared error among functions
constant on the $`J`$ equal-probability intervals, so
$`e\le\|g-g_T\|_2+L_T/J`$. For $`T\ge1`$,
$`\|g-g_T\|_2^2\le\mathbb E[X^2\mathbf 1_{|X|>T}]\le2(T+1)\phi(T)`$. Taking
$`T=\sqrt{\log J}`$ gives

```math
e^2\le\frac{4(\sqrt{\log J}+1)}{\sqrt{2\pi J}}+\frac{4\pi}{J}\le32J^{-1/3}.
\qquad\square
```

So $`J=\mathrm{poly}(1/\varepsilon)`$, independent of $`N`$, makes $`e`$ as small as
required; concretely a power of two with $`J\ge(3200/\delta^2)^3`$ gives
$`e^2\le\delta^2/100`$, and $`\log J=O(\log(1/\delta))`$.
`verify_koopman_transport.py` check (1) evaluates $`s^2,e^2`$ for several $`J`$ against
the bound.

---

## 3. The fixed encoding

Use $`2N`$ real coordinates. Let $`\Omega=[J]^{2N}`$, $`M=|\Omega|=J^{2N}`$, and for
$`a\in\Omega`$ define

```math
w_a=\frac{1}{s\sqrt2}\begin{pmatrix}q_{a_1}+i\,q_{a_{N+1}}\\ \vdots\\ q_{a_N}+i\,q_{a_{2N}}\end{pmatrix}\in\mathbb C^N,
\qquad
(E\psi)_a=\frac{w_a^{*}\psi}{\sqrt M}.
```

Independence of the coordinates, zero means and the normalization by $`s`$ give
$`\frac1M\sum_aw_aw_a^{*}=I_N`$, hence

```math
E^{*}E=I_N\quad\text{exactly.}
```

Check (2) of the verifier confirms this at $`N=1,2`$, $`J=4`$ to machine precision.

**Circuit for the encoder.** Define two orthonormal states on $`\log_2J`$ qubits,
$`|g\rangle=J^{-1/2}\sum_a|a\rangle`$ and $`|h\rangle=(s\sqrt J)^{-1}\sum_aq_a|a\rangle`$. Then

```math
E|j\rangle=\frac1{\sqrt2}\Bigl(|h\rangle_j\bigotimes_{k\ne j}|g\rangle_k\;-\;i\,|h\rangle_{N+j}\bigotimes_{k\ne N+j}|g\rangle_k\Bigr)
```

on $`2N`$ coordinate registers. A circuit: (i) convert the binary input label coherently
into a single excitation among $`2N`$ registers, the two locations $`j`$ and $`N+j`$
carrying relative phase $`-i`$; (ii) erase the binary label reversibly by reconstructing it
from the excitation's position; (iii) on every coordinate register apply the same fixed
$`J`$-dimensional unitary taking $`|0\rangle\mapsto|g\rangle`$, $`|1\rangle\mapsto|h\rangle`$;
(iv) prepare the label register of Section 5 with Hadamards. Step (i) takes
$`O(N(n+1)^2)`$ reversible gates; step (iii), after finite-precision compilation of one
fixed $`J`$-dimensional unitary, $`O(NJ^2\,\mathrm{poly}(\log J,\log(N/\varepsilon)))`$
gates. The quantiles and the resulting gates are numerical constants depending only on
$`J`$. Compile the encoder to operator error at most $`\varepsilon/8`$ and decode with its
circuit inverse. No target-dependent gate angle is supplied anywhere.

---

## 4. Transport as a doubly stochastic coupling

Fix the target $`U`$ and its realification

```math
O_U=\begin{pmatrix}\mathrm{Re}\,U&-\mathrm{Im}\,U\\ \mathrm{Im}\,U&\mathrm{Re}\,U\end{pmatrix}\in O(2N).
```

Draw $`X\sim\mathcal N(0,I_{2N})`$ and form the grid labels $`\alpha=\theta(X)`$, $`\beta=\theta(O_UX)`$ with
coordinatewise quantization. Both are uniform on $`\Omega`$, so

```math
\Xi_{ba}=M\,\Pr(\alpha=a,\,\beta=b)
```

is doubly stochastic and admits a Birkhoff decomposition
$`\Xi=\sum_{\ell=1}^{m}p_\ell P_{\pi_\ell}`$, $`m\le M^2`$, into permutation matrices. (Hall's
condition holds on the support; pick a perfect matching, subtract its smallest entry, and
repeat. Row and column sums stay equal and one positive entry disappears each round, so at
most $`M^2`$ rounds occur and the subtracted weights sum to one.)

**Lemma 4.1 (uniform coupling error).** *With $`\Delta_\pi:=P_\pi E-EU`$,*

```math
\sum_\ell p_\ell\,\Delta_{\pi_\ell}^{*}\Delta_{\pi_\ell}\;\preceq\;\frac{4e^2}{s^2}\,I_N .
\tag{1}
```

*Proof.* Let $`r(X)=X-q_{\theta(X)}`$ be the quantization residual. Its covariance is
$`e^2I_{2N}`$, and the same holds for $`r(O_UX)`$ since $`O_UX`$ is again standard normal.
Also $`q_{\theta(O_UX)}-O_Uq_{\theta(X)}=O_Ur(X)-r(O_UX)`$. Using
$`(a-b)(a-b)^{\top}\preceq2aa^{\top}+2bb^{\top}`$,

```math
\mathbb E\bigl[(q_{\theta(O_UX)}-O_Uq_{\theta(X)})(q_{\theta(O_UX)}-O_Uq_{\theta(X)})^{\top}\bigr]\preceq4e^2I_{2N}.
```

For real $`x,y`$, the map $`(x,y)\mapsto(x+iy)/(s\sqrt2)`$ has squared
real-linear norm $`1/(2s^2)`$. The covariance calculation instead uses the
complex-linear extension $`L=[I_N\;\;iI_N]/(s\sqrt2)`$, which satisfies
$`LL^*=I_N/s^2`$. Applying $`\Sigma\mapsto L\Sigma L^*`$ to the real covariance
bound therefore gives
$`\mathbb E\,(w_\beta-Uw_\alpha)(w_\beta-Uw_\alpha)^{*}\preceq(4e^2/s^2)I_N`$. At output coordinate
$`b=\pi(a)`$ the entry of $`\Delta_\pi\psi`$ is $`(w_a^{*}\psi-w_b^{*}U\psi)/\sqrt M`$;
averaging its squared magnitude gives the quadratic form of
$`U^{*}\mathbb E[(w_\beta-Uw_\alpha)(w_\beta-Uw_\alpha)^{*}]U`$.
Unitary conjugation preserves the preceding scalar identity bound. $`\square`$

Equation (1) is an operator bound, valid for every input at once, not an average-case
statement. No measurement of the input has entered. Check (3) of the verifier estimates the
left side by Monte Carlo at $`N=2`$, $`J=4`$ and confirms both the operator bound and the
identity between the coupling reading and the operator reading; check (5) forms the left
side from a real coupling and a real Birkhoff decomposition at $`N=1`$, $`J=4`$.

---

## 5. One coherent permutation from the mixture

Set $`\delta=\varepsilon/2`$ and take $`J`$ with $`e^2\le\delta^2/100`$. Choose a power of
two $`R\ge100M^2/\delta^2`$; since $`m\le M^2`$, $`R`$ is fixed by $`M`$ and $`\delta`$ alone
and does not depend on the target, which is what keeps the enlarged encoding below
target-independent. Approximate the weights by integer multiplicities $`k_\ell/R`$ with
$`\sum_\ell k_\ell=R`$ and $`\sum_\ell|p_\ell-k_\ell/R|\le2m/R`$ (round each weight down and
distribute the remaining slots). Make a list of $`R`$ permutations containing $`k_\ell`$
copies of $`\pi_\ell`$ and define one permutation of $`[R]\times\Omega`$:

```math
\Pi_U(j,a)=(j,\pi_{\ell(j)}(a)).
```

The enlarged, still target-independent, encoding is
$`\tilde E\psi=|+_R\rangle\otimes E\psi`$ with $`|+_R\rangle=R^{-1/2}\sum_j|j\rangle`$.

**Lemma 5.1 (coherent-label identity).**

```math
\bigl\|\Pi_U\tilde E-\tilde EU\bigr\|_{\mathrm{op}}^2
=\Bigl\|\sum_\ell\frac{k_\ell}{R}\,\Delta_{\pi_\ell}^{*}\Delta_{\pi_\ell}\Bigr\|_{\mathrm{op}} .
```

*Proof.* $`(\Pi_U\tilde E-\tilde EU)\psi=R^{-1/2}\sum_j|j\rangle\otimes\Delta_{\pi_{\ell(j)}}\psi`$,
and the $`|j\rangle`$ are orthonormal. $`\square`$

Since $`\|\Delta_\pi\|\le2`$, replacing $`p_\ell`$ by $`k_\ell/R`$ in (1) changes the left
side by at most $`4\sum_\ell|p_\ell-k_\ell/R|\le8m/R`$ in operator norm, so

```math
\bigl\|\Pi_U\tilde E-\tilde EU\bigr\|_{\mathrm{op}}^2\le\frac{4e^2}{s^2}+\frac{8m}{R}<\delta^2 .
\tag{2}
```

The label $`j`$ is never measured. Equation (2) says it returns, up to $`\delta`$, to the
same uniform state $`|+_R\rangle`$, uncorrelated with the input, and the inverse encoder
restores the coordinate registers. All target dependence sits in the classical permutation
$`\Pi_U`$. Check (4) of the verifier confirms Lemma 5.1 on a synthetic Birkhoff mixture
with $`M=256`$, $`R=8`$, and check (5) confirms it, together with the inequality above, on a
real coupling, a real Birkhoff decomposition and the real $`\Pi_U`$ at $`N=1`$, $`J=4`$,
$`R=2^{14}`$.

---

## 6. Two Boolean queries to one table

Let $`t=\log_2(MR)=O(N\log(1/\varepsilon))`$ be the encoded register width. Define one
Boolean table on a direction bit $`d`$ and two $`t`$-bit strings $`z,y`$:

```math
f_U(d,z,y)=\begin{cases}\langle\Pi_U(z),y\rangle_{\mathbb F_2},&d=0,\\ \langle\Pi_U^{-1}(z),y\rangle_{\mathbb F_2},&d=1,\end{cases}
\qquad\text{input length }\ell=1+2t .
```

With the answer qubit in $`|-\rangle`$ one query supplies the phase
$`(-1)^{\langle\Pi_U(z),y\rangle}`$; Hadamards before and after on the $`y`$ register
turn it into $`|z\rangle|b\rangle\mapsto|z\rangle|b\oplus\Pi_U(z)\rangle`$, coherently for
arbitrary superpositions. Starting from $`b=0`$: query the forward table to obtain
$`|z,\Pi_U(z)\rangle`$; swap the registers; query the inverse table to obtain
$`|\Pi_U(z),0\rangle`$. That is $`\Pi_U`$ with two Boolean queries, $`O(t)`$ extra gates
and $`O(t)`$ scratch qubits, both directions being tagged regions of the same table. Because
the direction bit places the state in the $`d=0`$ half of the table at the first query and in
the $`d=1`$ half at the second, one $`\pm1`$ sign word over $`(d,z,y)`$ serves both queries, and
the whole circuit is of the repeated-sign-word form $`T_g=BD_gWD_gC`$ of the paper's
Definition 2.1 with $`B,W,C`$ independent of $`U`$.

**Total error.** (2), the encoder error and the decoder error give isometry error at most
$`\delta+2(\varepsilon/8)=3\varepsilon/4`$; compressing to the target register is a
contraction, so the clean-output compression is within operator-norm distance
$`3\varepsilon/4`$ of $`U`$ as well. Two isometries at operator distance $`\gamma`$ induce
channels at half-diamond distance at most $`\gamma`$: tensoring with a reference preserves
the operator bound, the trace distance of the resulting pure states is at most their vector
distance, tracing out the environment does not increase it, and convexity covers mixed
inputs. The table is fixed before the
input arrives; the online circuit receives no advice and no target-dependent gate.

---

## 7. Converse for fixed encodings

**Theorem 7.1.** *Let $`E,F:\mathbb C^N\to\mathbb C^D`$ be arbitrary fixed isometries and let
the target-dependent action be any monomial unitary $`P|i\rangle=e^{i\theta_i}|\pi(i)\rangle`$.
If for every $`U\in U(N)`$ some $`P`$ brings $`PE`$ within half-diamond distance
$`\varepsilon`$ of $`FU`$, then $`\log_2D\ge(N-1)\log_2(1/\varepsilon)`$.*

*Proof.* Let $`e_i=E^{*}|i\rangle`$, $`f_j=F^{*}|j\rangle`$ and
$`a(U)=\max_{i,j:\,e_i,f_j\ne0}|\langle Ue_i/\|e_i\|,\,f_j/\|f_j\|\rangle|`$. For every
monomial $`P`$ simultaneously,

```math
|\mathrm{Tr}(U^{*}F^{*}PE)|\le a(U)\sum_i\|e_i\|\,\|f_{\pi(i)}\|\le a(U)\sqrt{\textstyle\sum_i\|e_i\|^2}\sqrt{\textstyle\sum_i\|f_{\pi(i)}\|^2}=N\,a(U),
```

so the normalized Choi vectors of $`PE`$ and $`FU`$ have squared overlap at most $`a(U)^2`$.
Channels within half-diamond distance $`\varepsilon`$ have pure Choi states within trace
distance $`\varepsilon`$, whence $`a(U)^2\ge1-\varepsilon^2`$. For fixed unit $`v,w`$ and
Haar-random $`U`$, $`\Pr(|\langle Uv,w\rangle|^2\ge1-\varepsilon^2)=\varepsilon^{2(N-1)}`$
(the first squared coordinate of a Haar vector is $`\mathrm{Beta}(1,N-1)`$, with upper tail
$`x^{N-1}`$ at $`1-x`$). A union bound over the at most $`D^2`$ pairs of row directions gives

```math
\Pr_{U\sim\mathrm{Haar}}\Bigl[\exists P:\ \tfrac12\|\mathrm{Ad}_{PE}-\mathrm{Ad}_{FU}\|_\diamond\le\varepsilon\Bigr]\le D^2\varepsilon^{2(N-1)} ,
```

and universality forces the left side to be one. $`\square`$

The bound makes no computational assumption on the encodings and already permits
continuous target-dependent phases. At $`\varepsilon=0.01`$ it is about $`6.64(N-1)`$
encoded qubits. The mechanism is transparent: a permutation can only rearrange the finite
collection of row directions the encoding supplies, and reproducing every rotation while
returning to a fixed output encoding needs exponentially many directions at the
prescribed error-dependent angular resolution.

---

## 8. What the converse leaves open

1. **A target-dependent final ancilla state.** Theorem 7.1 assumes a fixed final encoding
   including a fixed final ancilla state. It does not exclude a polynomial-width construction
   of the form $`WP_{f_U}E\approx(\psi\mapsto U\psi\otimes|\xi_U\rangle)`$ with $`W,E`$ fixed and
   $`|\xi_U\rangle`$ depending on the table but not on the input. There the output row
   directions depend on $`\xi_U`$ and the fixed-pair Haar computation no longer applies.
   Neither a construction nor a converse is known for that variant. This variant is a
   garbage-model question (paper, Section 7) and does not affect the clean model.
2. **Interleaved target-dependent mixing.** After earlier oracle-dependent operations the
   intermediate encoding may depend on $`U`$, and the union bound cannot simply be reapplied.
   In the clean two-query model this is exactly the case handled by the width theorem of the
   paper, Theorems 4.1 and 6.3, which need no fixed-encoding assumption.

The clean-model conclusion is therefore closed from both sides: Theorem 6.3 of the paper
gives $`\Omega(2^n)`$ address bits for any clean two-query circuit whatsoever, and Sections
1 to 6 above attain $`O(2^n\log(1/\varepsilon))`$.
