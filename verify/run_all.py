"""Run every verifier in this repository.  Exit 0 iff all pass.

    python -u verify/run_all.py

Every script is deterministic (fixed seeds), numpy-only, and finishes in seconds to a
couple of minutes on a laptop.  Each prints PASS/FAIL per check and exits nonzero on any
failure.  The suite is the machine-check layer of the paper; see STATUS.md for what each
script establishes and what it does not.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SUITE = [
    ("verify_proof_steps.py",
     "Corollary 4.2, audited real-selector form: exact expansion, pointwise absorption over "
     "the whole sign cube, Schur domination; Proposition 2.2 (corrected unitarity "
     "criterion); zero-leverage lemma (Remark 4.4)"),
    ("verify_complex_phase.py",
     "Theorem 4.1, unimodular-selector form: expansion, torus absorption, closed-form "
     "variance operators vs Monte Carlo and their domination, corrected T*T identity"),
    ("verify_permutation_family.py",
     "Proposition 6.1 -- the Dong-Lombardi-Ma permutation family lies in the class and "
     "fills a constant fraction of the ceiling"),
    ("verify_crho_requirement.py",
     "Proposition 5.3 -- the constant the permutation family forces"),
    ("verify_decisive_pair.py",
     "Proposition 3.1 -- the decisive pair (2 vs 16 targets, identical pre-final data)"),
    ("verify_address_compression.py",
     "Theorem 6.3 -- the address-respecting compression behind the input-length bound"),
    ("verify_koopman_transport.py",
     "Theorem 6.4 -- the matching upper bound: quantizer, isometry, coupling error, "
     "coherent-label identity, at small N"),
    ("verify_conjugation_transfer.py",
     "Section 7 -- junk-blindness lemma and the conjugation transfer identities"),
]

failed = []
for script, desc in SUITE:
    print("=" * 74)
    print(f"RUN  {script}")
    print(f"     {desc}")
    print("=" * 74)
    r = subprocess.run([sys.executable, "-u", os.path.join(HERE, script)], cwd=HERE)
    status = "PASS" if r.returncode == 0 else f"FAIL (exit {r.returncode})"
    print(f"\n>>> {script}: {status}\n")
    if r.returncode != 0:
        failed.append(script)

print("=" * 74)
if failed:
    print("SUITE FAILED: " + ", ".join(failed))
    sys.exit(1)
print(f"SUITE PASSED -- {len(SUITE)}/{len(SUITE)} verifiers green")
sys.exit(0)
