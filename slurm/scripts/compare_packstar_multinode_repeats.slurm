#!/bin/bash
# Compare two identical sequence-sharded multi-node PACK* runs.  Rank and
# wall-clock fields are ignored; all estimator/interval fields are compared
# by sequence key.
#
#SBATCH --job-name=packstar_mn_repeat
#SBATCH --account=grisman
#SBATCH --partition=grisman
#SBATCH --nodelist=grisman-37
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/usr/xtmp/lz280/slurm_logs/%x_%j.out
#SBATCH --error=/usr/xtmp/lz280/slurm_logs/%x_%j.err

set -euo pipefail

A0=${A0:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_0/1a0r_packstar.csv}
A1=${A1:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_1/1a0r_packstar.csv}
B0=${B0:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279817/rank_0/1a0r_packstar.csv}
B1=${B1:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279817/rank_1/1a0r_packstar.csv}
REPORT=${REPORT:-/usr/xtmp/lz280/packstar_production/packstar_multinode_repeat_1a0r_rb0.5.tsv}
DETAIL=${DETAIL:-/usr/xtmp/lz280/packstar_production/packstar_multinode_repeat_1a0r_rb0.5_detail.tsv}

for f in "$A0" "$A1" "$B0" "$B1"; do
    [ -f "$f" ] || { echo "missing CSV: $f" >&2; exit 2; }
done

awk -F, -v a0="$A0" -v a1="$A1" -v b0="$B0" -v b1="$B1" \
    -v report="$REPORT" -v detail="$DETAIL" '
    function finite(x) { return x ~ /^[-+]?[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ }
    function abs(x) { return x < 0 ? -x : x }
    FNR == 1 { next }
    {
        k=$2
        if (FILENAME == a0 || FILENAME == a1) {
            if (!(k in a)) a[k]=$0
        } else if (!(k in b)) {
            b[k]=$0
        }
    }
    END {
        shared=0; exact=0; scoreN=0; scoreDiffN=0; maxScoreDiff=0
        lbN=0; maxLBDiff=0; ubN=0; maxUBDiff=0
        intervalN=0; overlapN=0; jSum=0; jMin=1e300
        statusChange=0; epsN=0; epsDiffN=0; maxEpsDiff=0
        print "sequence\texactFields4to39\tscoreDiff\tlbDiff\tubDiff\tintervalJaccard\tACompStatus\tBCompStatus\tACompEps\tBCompEps" > detail
        for (k in a) {
            if (!(k in b)) continue
            shared++
            split(a[k], av, ",")
            split(b[k], bv, ",")
            same=1
            for (i=4; i<=39; i++) if (av[i] != bv[i]) same=0
            if (same) exact++
            sd=ld=ud="NA"
            if (finite(av[5]) && finite(bv[5])) { sd=abs(av[5]-bv[5]); scoreN++; if (sd>0) scoreDiffN++; if (sd>maxScoreDiff) maxScoreDiff=sd }
            if (finite(av[6]) && finite(bv[6])) { ld=abs(av[6]-bv[6]); lbN++; if (ld>maxLBDiff) maxLBDiff=ld }
            if (finite(av[7]) && finite(bv[7])) { ud=abs(av[7]-bv[7]); ubN++; if (ud>maxUBDiff) maxUBDiff=ud }
            jac="NA"
            if (finite(av[6]) && finite(av[7]) && finite(bv[6]) && finite(bv[7])) {
                lo=(av[6]>bv[6] ? av[6] : bv[6])
                hi=(av[7]<bv[7] ? av[7] : bv[7])
                inter=(hi>lo ? hi-lo : 0)
                uni=(av[7]>bv[7] ? av[7] : bv[7])-(av[6]<bv[6] ? av[6] : bv[6])
                jac=(uni>0 ? inter/uni : 1)
                intervalN++; if (inter>0) overlapN++; jSum+=jac; if (jac<jMin) jMin=jac
            }
            if (av[32] != bv[32]) statusChange++
            if (finite(av[33]) && finite(bv[33])) { epsN++; ed=abs(av[33]-bv[33]); if (ed>0.02) epsDiffN++; if (ed>maxEpsDiff) maxEpsDiff=ed } else ed="NA"
            print k "\t" same "\t" sd "\t" ld "\t" ud "\t" jac "\t" av[32] "\t" bv[32] "\t" av[33] "\t" bv[33] >> detail
        }
        print "metric\tvalue" > report
        print "runAUniqueSequences\t" length(a) >> report
        print "runBUniqueSequences\t" length(b) >> report
        print "sharedSequences\t" shared >> report
        print "exactEstimatorFieldsRows\t" exact >> report
        print "scoreComparisons\t" scoreN >> report
        print "scoreNonzeroDiffRows\t" scoreDiffN >> report
        print "maxAbsoluteScoreDiffLog10\t" maxScoreDiff >> report
        print "lbComparisons\t" lbN >> report
        print "maxAbsoluteLBDiffLog10\t" maxLBDiff >> report
        print "ubComparisons\t" ubN >> report
        print "maxAbsoluteUBDiffLog10\t" maxUBDiff >> report
        print "finiteIntervalPairs\t" intervalN >> report
        print "overlappingIntervalPairs\t" overlapN >> report
        print "meanIntervalJaccard\t" (intervalN ? jSum/intervalN : "NA") >> report
        print "minIntervalJaccard\t" (intervalN ? jMin : "NA") >> report
        print "complexStatusChanges\t" statusChange >> report
        print "complexEpsComparisons\t" epsN >> report
        print "complexEpsDiffOver0.02\t" epsDiffN >> report
        print "maxAbsoluteComplexEpsDiff\t" maxEpsDiff >> report
    }
' "$A0" "$A1" "$B0" "$B1"

echo "repeat report=$REPORT"
cat "$REPORT"
echo "repeat detail=$DETAIL"
