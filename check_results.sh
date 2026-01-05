#!/bin/bash

echo "Checking DP-MARKStar test results..."
echo ""

# Find the latest log file
LOG_FILE=$(ls -t dp_comparison_*.log 2>/dev/null | head -1)

if [ -z "$LOG_FILE" ]; then
    echo "No log file found yet."
    echo ""
    echo "Job status:"
    squeue -u $USER | grep dp_mark
else
    echo "Log file: $LOG_FILE"
    echo ""
    echo "=== Latest output ==="
    tail -50 "$LOG_FILE"
    echo ""
    echo "=== Job status ==="
    squeue -u $USER | grep dp_mark || echo "Job completed or not found"
fi

# Check if there's an error file
ERR_FILE=$(ls -t dp_comparison_*.err 2>/dev/null | head -1)
if [ -n "$ERR_FILE" ] && [ -s "$ERR_FILE" ]; then
    echo ""
    echo "=== Errors ==="
    tail -20 "$ERR_FILE"
fi
