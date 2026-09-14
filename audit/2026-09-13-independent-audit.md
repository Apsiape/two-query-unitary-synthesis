# Independent mathematical release audit — 2026-09-13

Scope: the current working tree, principally `paper/paper.tex`, its Markdown companion,
both `notes/` companions, `lean/TwoQuery/Q2Width.lean`, and the relevant verifiers.
This was a separate AI agent's fresh mathematical review within the same release
workflow, not external human peer review. Existing audit labels and prior review
verdicts were not premises. No manuscript or implementation file was edited by this
auditor; corrections described below were made by the main release agent and inspected
here. There was no `AGENTS.md` in the audited repository or its parent directory.

## Verdict

The headline clean two-query results survive this audit: the complex-selector width
bound, workspace-independent address compression, and the finite permutation-transport
upper bound have coherent analytic proofs with the claimed asymptotic resources.
I found no counterexample or missing load-bearing mathematical dependency in those
arguments. The unrestricted two-query garbage-model question remains open, as stated.

Two corrections were required before treating the public package as reviewed: an
ambiguous/wrong literal quantification of the transfer theorem's capacity hypothesis,
and a non-rigorous derangement approximation described as a rigorous finite bound.
Both were reported immediately and repaired during this audit. Neither requires changing
the headline thresholds. Additional wording issues and validation limits are recorded
below. This is an analytic audit, not a full formalization or an exhaustive priority search.

## Findings and disposition

### F1 — Fix the architecture before packing its outputs

Original locations: `paper/paper.tex:722`, Theorem 7.2; `paper/paper.md:669`;
`notes/garbage-to-clean.md`, sections 1 and 3.

The original displayed insertion class allowed the contractions `Y_i` to vary and then
asserted a capacity bound for that class. On the literal union-of-outputs reading, its
claimed application at `t=1` is false: take one address, `P_1=I`, all slots but `Y_0`
equal to the identity, and let `Y_0` range over `U(N)`. Two repeated signs cancel, so
the union contains all of `U(N)` with only one address. Its normalized packing entropy
is order `N^2`, not `O(N log N)`.

The required hypothesis is instead: for every fixed compatible tuple of contractions
and fixed address projectors, the family indexed only by `g` satisfies the bound, with
constants uniform over those architectural choices and workspace dimension. The proof
itself fixes those choices and is correct under that hypothesis. The conditional
implication was not refuted; its literal hypothesis and its advertised instantiation
were inconsistent.

Disposition: the revised tex/md statement explicitly fixes the tuple, packs over `g`
alone, and requires uniform constants. The note now distinguishes the class of maps
from each map's output family. These changes were inspected and resolve the substantive
issue. Also use `phi_N(Q)`, or explicitly permit dependence on `N`: the original wording
called `phi` a function only of address count but instantiated it as `log(2QN)`.
In the note, `{T_g} subset P_{2t}` should mean the output family of a member of `P_{2t}`,
since the latter is a collection of maps. These last two changes are notation repairs.

Release classification: required theorem-statement correction; no new analytic proof
needed. Resolved: the fixed-architecture quantifier, `phi_N` dependence, and the
note's map/output-family typing were all inspected after correction.

### F2 — The finite permutation calibration used a false derangement formula

Original locations: `verify/verify_crho_requirement.py:33` and `:52`;
paper Proposition 5.3's finite calibration table.

The old implementation used `log(!n)=log(n!)-1` for every positive `n`, while calling
the resulting Gilbert–Varshamov expression a rigorous lower bound. This is neither
exact nor consistently conservative: `!1=0`, `!2=1`, whereas the approximation gives
`1/e` and `2/e`. For `d=4`, the Hamming ball of radius two has exactly
`1 + binom(4,2)=7` permutations. The old code returned `1+16/e`, approximately 6.886,
and therefore overstated the corresponding greedy-packing lower bound.

Disposition: inspected the replacement recurrence
`!n=(n-1)(!(n-1)+!(n-2))`, initialized by `!0=1,!1=0` and evaluated with a stable
log-sum. It is correct. The code now explicitly distinguishes the exact recurrence
from floating-point evaluation without certified enclosures, and checks the small
derangements and the seven-element ball. The main release agent reports unchanged
table entries at the displayed precision; that numerical rerun belongs to the main
release verification, not to the independent test list below.

The asymptotic claim `kappa >= 1/4` does not depend on this numerical approximation.
For radius `r` proportional to `d`, the elementary upper bound
`sum_{j<=r} binom(d,j)!j <= (r+1)d!/(d-r)!` gives a GV logarithm at least
`log((d-r)!)-log(r+1)=(1-alpha)d log d-O(d)`. With `rho^2=2alpha` and
`Q=2d^2`, the required constant tends to at least `alpha(1-alpha)`;
choosing `alpha=1/2` gives `1/4`.

Release classification: required correction to the verifier's finite-bound claim;
resolved by inspected implementation and scope wording. No headline theorem defect.

### F3 — Label heuristic and physical objects accurately

Nonblocking presentation corrections:

- `verify/verify_complex_phase.py:16,163` called alternating maximization “the torus
  supremum.” It supplies candidate optima, with no global optimality certificate.
  Call them candidate values found by alternating maximization. The analytic
  pointwise inequality proves the true supremum bound independently of that routine.
- `paper/paper.tex:615`, `paper/paper.md:562`, and
  `notes/permutation-transport.md:261` call the operator-norm-close matrix the
  “realized map on the target register.” Precisely, this matrix is the clean-output
  compression `T_g`; the actual reduced physical evolution is a quantum channel.
  The separate full-isometry and diamond guarantees already handle physical leakage.
- In the transport construction, choosing the *least* power of two above the
  stated lower bound for `R` makes the claimed upper resource estimate immediate.
  Existence already permits this choice; it is not a gap in the upper-bound theorem.
- The transport note's phrase “directions exponentially close to every rotated row”
  (`notes/permutation-transport.md:304` before edits) is misleading at fixed error.
  The required closeness is set by epsilon; it is the number of directions that is
  exponential in `N`. The preceding formal converse is correct.

At the final read, the candidate-optimum wording and the tex/md clean-compression
wording were corrected. The transport note's two prose phrasings and the optional
least-dyadic clarification were left as nonblocking recommendations.

## Independent re-derivations of the load-bearing arguments

### Theorem 4.1 — complex selectors and constants

For each Gaussian realization, expanding the trace gives
`M_pq=W_pq sum_ij conjugate(G_ij) B_ip C_qj`. There is no selector in `M`.
After restricting to nonzero endpoint leverage supports, set
`K_pq=M_pq/sqrt(lambda_p mu_q)`. The bilinear form is `u^T K v`;
complex Cauchy–Schwarz applies as `(conjugate u)^* K v`, so it does not require
real selectors. On the full torus the two vector norms are exactly `||B||_F`
and `||C||_F`. Taking a supremum therefore costs no net or union bound.

For the Gaussian coefficient matrices
`(A_ij)_pq=W_pq B_ip C_qj/sqrt(lambda_p mu_q)`, direct coefficient summation gives

```
sum_ij A_ij A_ij* = conjugate(Gamma_B) o (W P_mu W*)
sum_ij A_ij* A_ij = conjugate(Gamma_C) o (W* P_lambda W).
```

Both conjugations are necessary and correct. A positive semidefinite correlation
matrix induces a contractive Schur multiplier: its Gram-factor dilation compresses
`I tensor X`. The two variance matrices are thus positive semidefinite and bounded
by their support identities. This remains valid with distinct supports, rectangular
middle matrices, strictly contractive endpoints, and endpoint rank below `N`.
Zero endpoints give the zero family and require no inverse leverage.

The complex Gaussian normalization is `E|G_ij|^2=1`. Splitting each coefficient into
real and imaginary Gaussian series of variance `1/2`, and using Tropp for each,
gives the stated factor `2 sqrt(log(2Q))`. Combining both real series at once gives
`sqrt(2 log(2Q))`, so the advertised constant 2 is conservative. Tropp's theorem
allows complex rectangular coefficient matrices with real Gaussian scalars; no
unmentioned realification dimension factor is necessary. I checked the actual
statement in [Tropp, Theorem 4.1.1](https://tropp.caltech.edu/books/Tro14-Introduction-Matrix-preprint.pdf).

### Corollaries 5.1 and 5.2 — metric, packing, and quantifiers

For `X_T=Re Tr(G*T)`, the increment variance is `||S-T||_F^2/2`. Consequently
Sudakov at Frobenius separation `rho sqrt(N)` gives
`log Pack <= 2 C_S^2 w_2^2/(rho^2 N)`, hence `kappa=8 C_S^2` with the quoted
width constant. This is a statement about each fixed architecture's outputs;
it is not a count over architectures, advice-dependent gates, or arbitrary metrics.

For Haar unitary packing, `f(V)=Re Tr(U*V)` has mean zero, Lipschitz constant
`sqrt(N)`, and the relevant ball is `f >= N(1-rho^2/2)`. The cited
[Meckes Theorem 5.17](https://case.edu/artsci/math/esmeckes/Haar_book.pdf)
does include `U(N)` and gives exponent `(N-2)t^2/(24L^2)`. For `N>=3`,
`N-2>=N/3` proves the manuscript's bound with an absolute constant.
Maximal separated sets cover; Haar invariance supplies a uniform ball measure.
Operator error below `rho/4` consumes at most half the pairwise separation.
Thus `log Q=Omega_rho(N)` follows without a change of metric or a hidden
supremum/expectation interchange. The constants are fixed before `N` grows.

This checks the stated operator/full-isometry model. A channel-distance claim
should use phase-invariant packing or the separate channel argument, not silently
identify matrices differing by a scalar phase. The `N=1` finite example is properly
identified as a count of matrices, not of channels.

### Theorem 6.3 — address compression independent of workspace

For fixed address projectors, the spaces
`L=direct-sum_j ran(P_j B*)` and `R=direct-sum_j ran(P_j C)` have dimensions
at most `MN`. Their definitions are target/selector independent. They are invariant
under every projector, contain `ran(B*)` and `ran(C)`, and satisfy
`B=B Pi_L`, `C=Pi_R C`. Inserting these projections on the two sides of `W`
therefore preserves the product exactly. The compressed middle is a contraction;
the endpoint Frobenius norms are preserved, not merely bounded.

Choosing address-adapted bases repeats each phase along its compressed block.
Relaxing those repeated phases to independent coordinates only enlarges the
supremum, so the rectangular width argument applies with dimension sum at most
`2MN`. The consequence `log M=Omega(N)` is asymptotic: the extra `log N` in
`log(2MN)` is lower order. It does not depend on the original register dimension.
Distinct address decompositions work with dimensions `M_1 N` and `M_2 N`.
For a bit-flip query, the response `|+>` sector contributes one fixed-phase
projector and the `|->` sector contributes `M` address projectors, giving the
claimed `M+1`. No workspace factor re-enters the address count.

### Theorem 6.4 and Appendix A — transport, uniformity, physical circuit

The quantizer estimate follows from best `L2` cellwise approximation, clipping,
and the Gaussian tail second moment. Substitution `T=sqrt(log J)` yields the
displayed bound; its elementary constant estimate is conservative. A least dyadic
`J >= (3200/delta^2)^3` suffices, independent of `N`. Normalization gives
`s^2=1-e^2>0`; coordinate independence and zero means make `E*E=I` exactly.

The binary-label-to-single-excitation construction is a reversible implementation
of that isometry. The signs `1/sqrt(2), -i/sqrt(2)` agree with the conjugated row
`w_a*`. A fixed unitary can map the two initial coordinate basis states to the
orthonormal pair `g,h`; standard finite-dimensional synthesis and compilation use
`O(N J^2 poly(log J,log(N/epsilon)))` gates. Preparation, binary-label erasure,
and their inverses need only the stated polynomial resources. None depends on `U`.

For the coupling, both quantized marginals are uniform because the realification
`O_U` is orthogonal. The finite coupling is doubly stochastic and has a Birkhoff
decomposition with at most `M^2` terms. Its residual vectors need not be independent:
the positive-semidefinite inequality for `(a-b)(a-b)^T` already gives covariance
at most `4 e^2 I`. Complexification has squared norm `1/s^2` and gives the
claimed complex covariance bound. At output address `b=pi(a)`, the discrepancy is
`w_a* psi - w_b* U psi`; this checks the orientation `P_pi E - E U`, including
the adjoints. Taking its second moment gives the operator bound simultaneously
for every input. No passage from average input error to operator error is used.

Round the mixture using a common dyadic `R` chosen from `M,delta` only. The
total variation in weights is at most `2m/R`; each squared error operator has
norm at most four. Orthogonal coherent labels turn the rounded average of squared
error operators into exactly the squared operator norm of one encoded permutation.
With `e^2<=delta^2/100`, `R>=100M^2/delta^2`, and `delta<0.05`, the squared
error is less than `0.121 delta^2`, comfortably below `delta^2`.
The label is uniform before encoding and approximately restored after transport;
it is not chosen according to the target's mixture weights and is not discarded.

For `t=log2(MR)`, the one Boolean table takes the complete address `(d,z,y)` of
length `1+2t`. Character orthogonality proves that Hadamards around its phase
query implement XOR of `Pi_U(z)` for any initial answer string. The physical
sequence is forward XOR, swap, inverse XOR; the second register becomes
`z XOR Pi_U^{-1}(Pi_U(z))=0`. The tag changes by a fixed gate, and a Boolean
bit-query answer in `|->` can be prepared and restored by fixed gates. Thus
there are exactly two Boolean queries to one table, not an assumed direct
permutation oracle. The second string, response qubit, tag, encoder scratch,
and label all fit `O(N log(1/epsilon))` qubits; the complete address has the
same order. No efficient offline construction or small truth-table size is claimed.

Encoder and decoder compilation add at most `epsilon/8` each to the transport
error `delta=epsilon/2`. The resulting full isometry is within `3epsilon/4`
of `psi -> U psi tensor |0...0>`. Tensoring with a reference preserves operator
norm, and pure-state trace distance is at most vector distance. Contractivity
and convexity give the stated half-diamond bound for arbitrary entangled inputs.
Approximate cleanup is proved in the stronger full-output metric.

The fixed-encoding monomial converse also checks: its normalized Choi overlap is
at most the largest normalized row-direction overlap `a(U)`. Half-diamond error
epsilon requires `a(U)^2>=1-epsilon^2`; the Haar coordinate tail and a union
bound over at most `D^2` fixed pairs give `D^2 epsilon^{2(N-1)}>=1`.
The rows cease to be fixed if final junk depends on the target; that restriction
is explicit and essential.

### Garbage transfer and higher-depth scope

The conjugated observable is `S_g*(R tensor I)S_g`. Expanding the adjoint reverses
the circuit and gives `2t` appearances of the same Hermitian sign oracle, with
fixed contractive slots. Normalized target-dependent junk cancels exactly.
For approximate discard synthesis, channel duality directly gives
`||Phi_g*(R)-U* R U||op <= epsilon` at unhalved diamond error epsilon; no
dimension-dependent Stinespring estimate is needed for the transfer.

For even `N`, balanced reflections form an orbit of real dimension `N^2/2`.
The function `Tr(R V* R V)` has mean zero and Lipschitz constant `2sqrt(N)`.
The same concentration/covering proof yields order-`exp(c N^2)` packing directly
on that orbit. This avoids incorrectly pushing an arbitrary unitary packing
through a map with a large stabilizer. With the repaired fixed-architecture
hypothesis, the conditional contradiction follows. At `t=1` compression supplies
the hypothesis; at `t=2` it would require a four-insertion estimate not proved here.
The phase argument for vanishing junk-blind linear functionals does not establish
a minimum oracle degree for all invariants; the current manuscript correctly
disclaims that conclusion.

## Primary-source and formalization checks

Targeted checks support the current attribution/model comparisons:

- [LMW, sections 1 and 3 and Appendix B](https://arxiv.org/html/2310.08870v1)
  supply the one-query address-size lower bound, the parallel-query reduction,
  and the discussion of limitations of simple counting. The present clean
  threshold should not be advertised as resolving their discard-model two-query problem.
- [DLM, Question 1.4](https://arxiv.org/html/2607.26478v1) asks whether `t`
  alternations need `t` sequential queries and proves the `t=2` case. The
  manuscript's indexing explanation is correct and does not claim their explicit family.
- [Brakerski–Yuen, introduction and section 5](https://arxiv.org/html/2605.09957v1)
  support the all-query-count oracle input lower bound and the cited synthesis model.
- [Rosenthal, Theorems 2.2 and 3.2 and their composition](https://arxiv.org/html/2111.07992v5)
  support `O(sqrt(N))` classical queries with polynomial-in-`n` workspace and
  approximate cleanup. The `2^{O(n)}` ancillas in the separate circuit-depth
  result must not be substituted for this query result. The
  [Quantum publication record](https://quantum-journal.org/papers/q-2026-06-30-2144/)
  confirms the 2026 article citation.
- [Banerjee's report](https://www.scottaaronson.com/showcase5/synthesis.pdf)
  explicitly assumes exact cleanup for every oracle and states its sharpened
  `2^N` count. [Huang's record](https://arxiv.org/abs/2508.13215) confirms the
  cited 1.5-query paper. Neither was independently re-proved here or used as
  a dependency. This audit does not certify exhaustive novelty or validate all
  proofs of every cited paper.

Reading the Lean declarations confirms that `width_of_khintchine` and
`asymmetric_width_assembly` are scalar implications from explicit absorption and
concentration hypotheses. They do not formalize those analytic inputs, the Schur
variance computation, address compression, packing, or transport. The seven scoped
declarations are consistent with the disclosure in `lean/README.md`.
I did not rerun Lean's kernel in this audit; recorded axiom output is a
historical artifact, not a fresh check by this auditor. The main release build
and kernel-check status must remain separately reported.

## Bounded diagnostics actually run in this audit

All processes used `../ak-working/tools/run_ak_capped.py` (relative to the audited
repository; an external local safety wrapper, not a package dependency), 512 MiB and a
60-second timeout, one diagnostic at a time. Every process exited successfully;
none remained running.

1. `verify/verify_exact_controls.py`: exact rational 256-word control, query
   cancellation control, and non-scalar `N=2` transport/Boolean lookup controls passed.
2. `verify/verify_address_compression.py`: 80 projector decompositions and 40
   distinct-decomposition cases passed; maximum residuals were approximately
   `7.3e-16` and `6.8e-16`.
3. A separate inline NumPy diagnostic directly formed all Gaussian coefficient
   matrices and summed both variance operators, rather than relying on Monte Carlo.
   Cases `(N,L,R)=(1,1,1),(2,3,7),(5,2,3),(3,8,5)`, twelve trials each,
   included distinct zero-leverage supports and rectangular middle layers.
   Both conjugated Schur identities, variance domination, and sampled complex
   pointwise absorption passed; maximum coefficient-identity residual was
   `6.66e-16`. This is a finite algebra check, not an asymptotic certificate.

No full numerical suite was rerun by this auditor. In particular, successful samples
or alternating maximization cannot certify a continuous supremum or the transport
asymptotics; the analytic arguments above carry those conclusions.

## Reviewed snapshot provenance

The initial independently read working-tree versions had these SHA-256 hashes:

```
paper/paper.tex                  E39591FDCFB5F31A95A65DD9D06C914797F4E144F1D8649E64F509F3C293C9BB
paper/paper.md                   967E95892A65BA41D526C6DF25759148447D4081DA45BE725CC12BF03DE69CBE
notes/permutation-transport.md   323A74EFD1A06A0CDB9580AE738C8B1C54218587C6793369C015A30C2BA5644B
notes/garbage-to-clean.md         34BC78A1A8568ED1CA539A8A4D656C539133EAB18734561716EF656BDF2FFAF3
lean/TwoQuery/Q2Width.lean        26714C25FF96331C44AD83D15A93FFA2154BB23C9DD5F93836B44310319B5A5E
```

After the inspected substantive corrections, no mathematical release blocker remains
within this audit's scope. The main release agent separately reported a successful
nine-script suite and Lean build; these are not represented as tests run by this
auditor. The Lean axiom-query result was still pending when that report was received.

The corrected versions inspected before final handoff included:

```
paper/paper.tex                  A2F8A1B417410A4014E5D2BE4FA354D5BE9A9E0A6B7750529FA5C569C4DE0263
paper/paper.md                   8A6DD25E0CC37E1E3D77AEC4EF6466E9D77D753C5EE1338A485C246A49853CCB
notes/garbage-to-clean.md         8DB80B79E10A4E84E1963ABEB91F359A2E71652F3746D70324B5C53FE9E9D9BB
verify/verify_crho_requirement.py 597B40A63675627964FD1A96A981B276F264D64DD773C1BC7DA138DEA3128207
```

The release agent's scoped corrections were inspected as described above.
Later edits require their own change-driven check; this report does not certify
an unseen final PDF, archive, commit, repository publication, or Zenodo record.
