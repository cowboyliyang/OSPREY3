#!/bin/bash

JOB_ID=10340523

echo "=========================================="
echo "Monitoring DP-MARKStar Large Test"
echo "Job ID: $JOB_ID"
echo "=========================================="
echo ""

# Check job status
echo "=== Job Status ==="
squeue -j $JOB_ID || echo "Job completed or not found"
echo ""

# Show latest log
LOG_FILE="dp_comparison_${JOB_ID}.log"
if [ -f "$LOG_FILE" ]; then
    echo "=== Latest Log Output (last 30 lines) ==="
    tail -30 "$LOG_FILE"
else
    echo "Log file not created yet: $LOG_FILE"
fi

echo ""
echo "=========================================="
echo "To keep watching, run: watch -n 10 ./watch_large_test.sh"
echo "=========================================="
