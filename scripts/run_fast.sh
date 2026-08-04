#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
(cd "$ROOT/cases/n16" && python3 reinhardt_n16_verified_audit.py)
(cd "$ROOT/cases/n32" && python3 reinhardt_n32_candidate_verifier.py)
(cd "$ROOT/cases/n64" && python3 n64_analytic_verifier.py)
(cd "$ROOT/cases/n64" && python3 n64_post_verifier.py)
printf '\nFAST VERIFICATIONS COMPLETED\n'
