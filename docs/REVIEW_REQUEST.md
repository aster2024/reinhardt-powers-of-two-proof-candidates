# Request for expert review

These proof candidates would settle previously open cases of a century-old extremal geometry problem if correct. They have not yet received independent expert mathematical verification. The most useful review is therefore adversarial, line-specific, and should treat every new lemma as untrusted until checked.

## Highest-priority mathematical questions

### Common to n=16, 32, and 64

1. **Difference-body reconstruction.** Does the closure condition on the selected antipodal edge vectors suffice exactly as stated to reconstruct a strictly convex small polygon with the claimed difference body?
2. **Perturbation feasibility.** When vertices of the difference body are perturbed subject to closure, are cyclic order, genuine-vertex inequalities, and feasibility for the original polygon preserved in every case used by the proof?
3. **Saturation.** Is the argument excluding the final interior half-vertex valid under the stated MFCQ/KKT hypotheses, including cyclic endpoint cases?
4. **Near-regular localization.** Are the pointwise angle-gap conclusions genuinely implied by the tiny global perimeter deficit with the constants used?
5. **Fixed-code uniqueness.** Is the strong-convexity/twisted-Dirichlet argument sufficient to prove that the high-perimeter KKT point is unique, rather than merely locally isolated?
6. **Symmetry quotient.** Do the dihedral code orbits correspond exactly to congruence classes under the allowed relabelings and reflections?

### n=16

- Replace or independently certify the geometric feasibility of Bingane's explicit \(C_{16}\) lower-bound construction from coordinates, not only its perimeter formula.

### n=32

- Reimplement the narrow interval inequalities with Arb/FLINT, MPFI, or another directed-rounding system. The actual gap is about \(1.33546\times10^{-13}\) against a \(1.35\times10^{-13}\) screen.

### n=64

- Check the proof that the ternary reflection-pair encoding is a bijective cover of all \(2^{64}\) half-codes and contains no hidden axial-symmetry assumption.
- Bundle and validate the generator for the 64 fixed-point trigonometric weights used by the C++ scan.
- Confirm the corrected residual padding: the core post-screen now uses `256`, which dominates the generic bound `128*sqrt(2)` and passes all five exclusions.
- Scrutinize the normal-cone localization and twisted discrete spectral bound, which are the most novel analytic steps.

## How to report an issue

Please include:

- case and file path;
- lemma/equation or code line;
- a precise objection, counterexample, or failed command;
- platform, compiler/interpreter version, and complete output for computational failures.

A confirmed flaw should be recorded even if it appears repairable. Silence or successful execution is not a substitute for reviewing the mathematical reduction.
