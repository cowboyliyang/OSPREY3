#!/bin/bash
# Compare the completed single-node and sequence-sharded multi-node PACK*
# outputs for the same design/budget. Run only after both formal jobs finish.
#
#SBATCH --job-name=packstar_compare
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodelist=jerry1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

SINGLE=${SINGLE:-/usr/xtmp/lz280/packstar_production/formal_singlenode_gpu/1a0r_rb0.5_12279812/rank_0/1a0r_packstar.csv}
MULTI0=${MULTI0:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_0/1a0r_packstar.csv}
MULTI1=${MULTI1:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_1/1a0r_packstar.csv}
REPORT=${REPORT:-/usr/xtmp/lz280/packstar_production/compare_1a0r_rb0.5_single_vs_multi.tsv}

for f in "$SINGLE" "$MULTI0" "$MULTI1"; do
    [ -f "$f" ] || { echo "missing CSV: $f" >&2; exit 2; }
done

awk -F, -v single="$SINGLE" -v multi0="$MULTI0" -v multi1="$MULTI1" '
function finite(x) { return x ~ /^[-+]?[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ }
function abs(x) { return x < 0 ? -x : x }
function add(path,    k) {
    k = $2
    if (path == single) {
        sScore[k]=$5; sLB[k]=$6; sUB[k]=$7; sProtEps[k]=$11; sCompEps[k]=$33
        sSeen[k]=1
    } else if (!(k in mSeen)) {
        mScore[k]=$5; mLB[k]=$6; mUB[k]=$7; mProtEps[k]=$11; mCompEps[k]=$33
        mSeen[k]=1
    }
}
FNR == 1 { next }
{
    path = FILENAME
    add(path)
}
END {
    shared=0; singleOnly=0; multiOnly=0; scoreCompared=0; scoreDiffs=0
    maxScoreDiff=0; overlap=0; finiteIntervals=0; protEpsDiffs=0; compEpsDiffs=0
    for (k in sSeen) {
        if (!(k in mSeen)) { singleOnly++; continue }
        shared++
        sd=""; if (finite(sScore[k]) && finite(mScore[k])) {
            sd=abs(sScore[k]-mScore[k]); scoreCompared++
            if (sd > maxScoreDiff) maxScoreDiff=sd
            if (sd > 0.1) scoreDiffs++
        }
        if (finite(sLB[k]) && finite(sUB[k]) && finite(mLB[k]) && finite(mUB[k])) {
            finiteIntervals++
            if (sLB[k] <= mUB[k] && mLB[k] <= sUB[k]) overlap++
        }
        if (finite(sProtEps[k]) && finite(mProtEps[k]) && abs(sProtEps[k]-mProtEps[k]) > 0.02) protEpsDiffs++
        if (finite(sCompEps[k]) && finite(mCompEps[k]) && abs(sCompEps[k]-mCompEps[k]) > 0.02) compEpsDiffs++
        printf "%s\t%s\t%s\t%s\t%s\t%s\n", k, sScore[k], mScore[k], sd, sLB[k] "/" sUB[k], mLB[k] "/" mUB[k] > "/dev/stderr"
    }
    for (k in mSeen) if (!(k in sSeen)) multiOnly++
    print "metric\tvalue"
    print "singleSequenceRows\t" length(sSeen)
    print "multiUniqueSequenceRows\t" length(mSeen)
    print "sharedSequences\t" shared
    print "singleOnlySequences\t" singleOnly
    print "multiOnlySequences\t" multiOnly
    print "finiteScoreComparisons\t" scoreCompared
    print "scoreDiffOver0.1Log10\t" scoreDiffs
    print "maxAbsoluteScoreDiffLog10\t" maxScoreDiff
    print "finiteIntervalComparisons\t" finiteIntervals
    print "overlappingLBUBIntervals\t" overlap
    print "proteinEpsDiffOver0.02\t" protEpsDiffs
    print "complexEpsDiffOver0.02\t" compEpsDiffs
}
' "$SINGLE" "$MULTI0" "$MULTI1" > "$REPORT"

echo "comparison report=$REPORT"
cat "$REPORT"
