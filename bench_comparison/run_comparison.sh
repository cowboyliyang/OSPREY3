#!/bin/bash
# ===================================================================
# Full comparison pipeline: K* vs MARK* vs GNN S9
#
# Usage:
#   bash run_comparison.sh d004           # single design
#   bash run_comparison.sh d001 d004 d019 # multiple designs
#   bash run_comparison.sh small          # preset: 5 small designs
#   bash run_comparison.sh medium         # preset: 10 medium designs
# ===================================================================
set -euo pipefail
BENCHDIR="$(cd "$(dirname "$0")" && pwd)"

case "${1:-small}" in
    small)
        # Small designs (2-3 mutable, few flex) — fast to complete
        DESIGNS=(d001 d004 d010 d035 d036)
        ;;
    medium)
        # Medium designs (4-6 mutable+flex)
        DESIGNS=(d005 d007 d008 d009 d012 d016 d020 d034 d047 d052)
        ;;
    large)
        # Large designs (7+ residues, 10^7+ conformations)
        DESIGNS=(d019 d023 d026 d041 d048 d049 d050 d051)
        ;;
    *)
        DESIGNS=("$@")
        ;;
esac

echo "============================================="
echo "  K* vs MARK* vs GNN S9 Comparison"
echo "  Designs: ${DESIGNS[*]}"
echo "============================================="

# Step 1: Train GNN for all designs
echo ""
echo "=== Step 1: Submit GNN training ==="
for did in "${DESIGNS[@]}"; do
    bash "$BENCHDIR/train_gnn.sh" "$did"
done

echo ""
echo "=== Step 2: Submit K* and MARK* baselines (can run in parallel with training) ==="
for did in "${DESIGNS[@]}"; do
    bash "$BENCHDIR/bench_kstar.sh" "$did"
    bash "$BENCHDIR/bench_markstar.sh" "$did"
done

echo ""
echo "=== Step 3: GNN S9 must wait for training to finish ==="
echo "After training completes, run:"
echo "  for did in ${DESIGNS[*]}; do bash $BENCHDIR/bench_gnn_s9.sh \$did; done"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Results: /usr/xtmp/lz280/bench_comparison/results/"
