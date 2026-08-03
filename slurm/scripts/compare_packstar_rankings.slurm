#!/bin/bash
# Compare global PACK* ranking order between one single-node CSV and the
# union of two sequence-sharded multi-node CSVs.  The Java writer sorts by
# lowerBoundLog10, so this report deliberately ignores the rank column in
# each shard and reconstructs the global order from lb_log10.
#
#SBATCH --job-name=packstar_rankcmp
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

SINGLE=${SINGLE:-/usr/xtmp/lz280/packstar_production/formal_singlenode_gpu/1a0r_rb0.5_12279812/rank_0/1a0r_packstar.csv}
MULTI0=${MULTI0:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_0/1a0r_packstar.csv}
MULTI1=${MULTI1:-/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_1/1a0r_packstar.csv}
REPORT=${REPORT:-/usr/xtmp/lz280/packstar_production/packstar_ranking_1a0r_rb0.5_single_vs_multi.tsv}
DETAIL=${DETAIL:-/usr/xtmp/lz280/packstar_production/packstar_ranking_1a0r_rb0.5_single_vs_multi_detail.tsv}

for f in "$SINGLE" "$MULTI0" "$MULTI1"; do
    [ -f "$f" ] || { echo "missing CSV: $f" >&2; exit 2; }
done

tmp=$(mktemp -d /usr/xtmp/lz280/packstar_rankcmp.XXXXXX)
trap 'rm -rf "$tmp"' EXIT

awk -F, '
    FNR == 1 { next }
    { print $2 "\t" $6 "\t" $5 "\t" $32 "\t" $33 }
' "$SINGLE" > "$tmp/single.raw"

awk -F, '
    FNR == 1 { next }
    !seen[$2]++ { print $2 "\t" $6 "\t" $5 "\t" $32 "\t" $33 }
' "$MULTI0" "$MULTI1" > "$tmp/multi.raw"

tab=$(printf '\t')
sort -t "$tab" -k2,2gr -k1,1 "$tmp/single.raw" |
    awk -F'\t' '{ print $1 "\t" NR "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' > "$tmp/single.rank"
sort -t "$tab" -k2,2gr -k1,1 "$tmp/multi.raw" |
    awk -F'\t' '{ print $1 "\t" NR "\t" $2 "\t" $3 "\t" $4 "\t" $5 }' > "$tmp/multi.rank"

awk -F'\t' -v report="$REPORT" -v detail="$DETAIL" '
    NR == FNR {
        sr[$1]=$2; slb[$1]=$3; ss[$1]=$4; sst[$1]=$5; se[$1]=$6; sn++
        next
    }
    {
        mr[$1]=$2; mlb[$1]=$3; ms[$1]=$4; mst[$1]=$5; me[$1]=$6; mn++
    }
    function abs(x) { return x < 0 ? -x : x }
    END {
        shared=0; exactRank=0; sumD2=0; orderSame=0; orderPairs=0
        top1=0; top5S=0; top5M=0; top5Intersect=0; top10Intersect=0
        multiNotEstimated=0; singleNotEstimated=0
        print "sequence\tsingleRank\tmultiRank\trankDelta\tsingleLB\tmultiLB\tsingleScore\tmultiScore\tsingleCompStatus\tmultiCompStatus\tsingleCompEps\tmultiCompEps" > detail
        for (k in sr) {
            if (!(k in mr)) continue
            shared++
            d=sr[k]-mr[k]; sumD2 += d*d
            if (d == 0) exactRank++
            if (sr[k] <= 5) top5S++
            if (mr[k] <= 5) top5M++
            if (sr[k] <= 10 && mr[k] <= 10) top10Intersect++
            if (sr[k] <= 5 && mr[k] <= 5) top5Intersect++
            if (sr[k] == 1 && mr[k] == 1) top1=1
            if (sst[k] != "Estimated") singleNotEstimated++
            if (mst[k] != "Estimated") multiNotEstimated++
            print k "\t" sr[k] "\t" mr[k] "\t" d "\t" slb[k] "\t" mlb[k] "\t" ss[k] "\t" ms[k] "\t" sst[k] "\t" mst[k] "\t" se[k] "\t" me[k] >> detail
        }
        for (a in sr) if (a in mr) for (b in sr) if (b in mr && a < b) {
            orderPairs++
            if ((sr[a] < sr[b]) == (mr[a] < mr[b])) orderSame++
        }
        rho = (shared > 1 ? 1 - 6*sumD2/(shared*(shared*shared-1)) : 1)
        print "metric\tvalue" > report
        print "singleRows\t" sn >> report
        print "multiUniqueRows\t" mn >> report
        print "sharedSequences\t" shared >> report
        print "exactRankPositions\t" exactRank >> report
        print "exactRankFraction\t" (shared ? exactRank/shared : 0) >> report
        print "sameTop1\t" (top1 ? "true" : "false") >> report
        print "top5Intersection\t" top5Intersect >> report
        print "top10Intersection\t" top10Intersect >> report
        print "pairwiseOrderAgreement\t" (orderPairs ? orderSame/orderPairs : 1) >> report
        print "spearmanRankCorrelation\t" rho >> report
        print "singleNonEstimatedOverall\t" singleNotEstimated >> report
        print "multiNonEstimatedOverall\t" multiNotEstimated >> report
    }
' "$tmp/single.rank" "$tmp/multi.rank"

echo "ranking report=$REPORT"
cat "$REPORT"
echo "ranking detail=$DETAIL"
