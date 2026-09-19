#!/usr/bin/env bash
set -euo pipefail

PROJECT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT"
RSCRIPT="${RSCRIPT:-/home/research/miniconda3/envs/paper2-r442/bin/Rscript}"
PYTHON="${PYTHON:-python3}"
RUN_DIR="results/paper_run_$(date +%Y%m%d_%H%M%S)"

"$PYTHON" scripts/run_paper.py \
  --project-root "$PROJECT" \
  --output-dir "$RUN_DIR" \
  --rscript "$RSCRIPT" \
  --dry-run

"$PYTHON" scripts/run_paper.py \
  --project-root "$PROJECT" \
  --output-dir "$RUN_DIR" \
  --rscript "$RSCRIPT" \
  --bootstrap 2000 \
  --selection-bootstrap 1000 \
  --seed 20260915

printf '\nResults: %s/%s/paper_results.zip\n' "$PROJECT" "$RUN_DIR"
