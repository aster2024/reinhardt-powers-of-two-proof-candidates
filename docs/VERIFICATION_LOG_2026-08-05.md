# Verification log — 2026-08-05

The following checks were executed in the packaging environment:

- `make verify-fast`: passed;
- n=16 exact certificate: passed, 32768 normalized codes, 16 survivors;
- n=32 exact certificate: passed, 2^31 normalized codes covered, 96 survivors;
- n=64 analytic certificate: passed;
- n=64 C++ exhaustive scan: passed, 896 survivors;
- n=64 original SHA-256 manifest: passed;
- n=64 post-screen certificate: passed, six orbits and all exclusion/uniqueness assertions passed.

The standalone n=64 scan used approximately 1.47 GB peak resident memory and completed in about 22 seconds on the packaging host. Runtime is machine-dependent and is not a proof assertion.
