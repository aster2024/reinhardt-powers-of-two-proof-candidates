# Audit artifacts

The files in this directory are separately organized computational cross-checks, not the original formal certificates and not independent human peer review.

- `audit_report_n16_n32_zh.md`: consolidated audit report for the n=16 and n=32 candidates.
- `n16/independent_scan.py`: a direct NumPy rescan of the normalized sign codes.
- `n32/independent_scan.py`: a differently partitioned meet-in-the-middle scan.
- `n32/high_precision_analysis.py`: high-precision closure, spectrum, and KKT reproduction.
- `n64/independent_scan.cpp`: an even/odd split rather than the original consecutive split.
- `n64/supplemental_checks.py`: fixed-point weight, padding, and constant checks.
- `n64/kkt_high_precision.py`: 100-digit full KKT reproduction.

Install optional dependencies with:

```bash
make audit-deps
```

High-precision numerical agreement is supporting evidence. The exact proof decisions remain in `cases/`.
