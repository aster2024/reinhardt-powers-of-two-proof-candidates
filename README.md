# Computer-assisted proof candidates for Reinhardt polygons at n=16, 32, and 64

This repository contains candidate computer-assisted proofs and executable certificates for the previously open power-of-two cases

\[
n=16,\qquad n=32,\qquad n=64
\]

of Reinhardt's maximum-perimeter problem for convex small polygons of diameter at most one.

> **Status (7 August 2026):** these are proof candidates, not peer-reviewed theorems. The packaged core certificates and separately organized computational cross-checks pass, while independent expert review of the analytic bridges is still needed.

中文说明：[README.zh-CN.md](README.zh-CN.md)

Archived releases are available through the Zenodo concept DOI
[10.5281/zenodo.21796494](https://doi.org/10.5281/zenodo.21796494).
The current archived version is
[10.5281/zenodo.21796495](https://doi.org/10.5281/zenodo.21796495).
Changes in this working tree postdate that immutable snapshot and require a new
Zenodo version before release.

## Results included

| Case | Candidate maximum perimeter | Exhaustive search | Survivors |
|---|---:|---:|---:|
| \(n=16\) | `3.136547716486607386085967...` | all \(2^{15}=32768\) normalized codes | 16, one dihedral orbit |
| \(n=32\) | `3.1403311569546193658254013805774586723...` | all \(2^{31}\) normalized codes | 96, three orbits |
| \(n=64\) | `3.1412772509327728680619914155024682980...` | all \(2^{64}\) half-codes via exact ternary meet-in-the-middle | 896, six orbits |

The candidate proofs claim uniqueness up to congruence, reflection, and cyclic relabeling.

## Run the certificates

Requirements:

- Python 3.10 or later;
- a C++17 compiler for the full \(n=64\) scan;
- about 2 GB of available RAM for the full \(n=64\) scan.

Lightweight verification:

```bash
make verify-fast
```

Full verification, including regeneration of the 896 \(n=64\) survivors:

```bash
make verify-all
```

The core certificate decisions use exact rational/integer interval arithmetic. For \(n=64\), the analytic core certifies all 64 C++ fixed-point weights and both endpoint Jensen bounds, while the post-screen uses a conservative 256-unit residual padding. Files under `audits/` are separate implementation and high-precision cross-checks.

## Repository layout

```text
cases/      proof candidates, verifiers, and recorded outputs
audits/     separate scans and high-precision checks
docs/       status, review targets, AI disclosure, and references
scripts/    reproducibility helpers
```

## Review status

The computational results have been reproduced, but the most important remaining task is expert review of the analytic bridges connecting the finite certificates to unconditional global optimality and uniqueness. See [docs/REVIEW_REQUEST.md](docs/REVIEW_REQUEST.md).

## Provenance

Generative-AI tools provided substantial assistance during proof exploration, implementation of searches and verifiers, computational cross-checking, and manuscript preparation. The human authors directed the research, assessed the resulting materials, and take responsibility for the released claims and artifacts. AI-assisted checks are not treated as independent peer review. See [docs/AUTHORSHIP_AND_AI_DISCLOSURE.md](docs/AUTHORSHIP_AND_AI_DISCLOSURE.md).

## Authors and maintainer

- **Jizhou Guo** — [ORCID](https://orcid.org/0009-0001-0699-9164), [Google Scholar](https://scholar.google.com/citations?user=fcBDdsYAAAAJ), [DBLP](https://dblp.org/pid/378/4049.html), [homepage](https://aster2024.github.io/), [X / Twitter](https://twitter.com/TheOsmanthus)
- **Yitao Luo** — [luoyt.kd@mail.ustc.edu.cn](mailto:luoyt.kd@mail.ustc.edu.cn)

Jizhou Guo maintains the repository.

## License

MIT License. See [LICENSE](LICENSE).

## References

- K. Reinhardt, *Extremale Polygone gegebenen Durchmessers*, 1922.
- C. Bingane, *Maximal perimeter and maximal width of a convex small polygon*, arXiv:2106.11831.
- B. Mulansky and A. Potschka, *A zonogon approach for computing small convex polygons of maximum perimeter*, *Mathematical Programming* (2025), DOI: 10.1007/s10107-025-02244-x, with correction DOI: 10.1007/s10107-025-02257-6.

Machine-readable references are in [docs/references.bib](docs/references.bib).
