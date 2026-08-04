#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# Run the memory-heavy scan first and buffer its small log so low-memory
# environments do not retain previous verifier output while the scan peaks.
printf '\n== n=64 exhaustive C++ scan ==\n'
(
  cd "$ROOT/cases/n64"
  g++ -O3 -std=c++17 n64_code_fixed.cpp -o n64_code_fixed
  MALLOC_ARENA_MAX=1 MALLOC_TRIM_THRESHOLD_=0 ./n64_code_fixed \
    >"$TMP/n64_scan.stdout" 2>"$TMP/n64_scan.stderr"
)
cat "$TMP/n64_scan.stderr"
cat "$TMP/n64_scan.stdout"

printf '\n== n=64 original manifest check ==\n'
(cd "$ROOT/cases/n64" && sha256sum -c SHA256SUMS)

printf '\n== n=64 analytic certificate ==\n'
(cd "$ROOT/cases/n64" && python3 n64_analytic_verifier.py)

printf '\n== n=64 post-screen certificate ==\n'
(cd "$ROOT/cases/n64" && python3 n64_post_verifier.py)

printf '\n== n=16 exact certificate ==\n'
(cd "$ROOT/cases/n16" && python3 reinhardt_n16_verified_audit.py)

printf '\n== n=32 exact certificate ==\n'
(cd "$ROOT/cases/n32" && python3 reinhardt_n32_candidate_verifier.py)

printf '\nALL CORE VERIFICATIONS COMPLETED\n'
