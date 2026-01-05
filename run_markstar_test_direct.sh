#!/bin/bash
cd /home/users/lz280/IdeaProjects/OSPREY3

echo "=========================================="
echo "Testing MARK* vs K* Partition Function"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Run the test with console output
./gradlew test --tests "edu.duke.cs.osprey.kstar.TestMARKStarVsKStarPartitionFunction.testSmallSystem3Flex" --console=plain 2>&1

echo ""
echo "=========================================="
echo "End time: $(date)"
echo "Exit code: $?"
echo "=========================================="
