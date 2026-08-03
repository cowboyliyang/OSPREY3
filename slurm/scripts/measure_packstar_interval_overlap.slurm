#!/bin/bash
# Measure LB/UB interval overlap for the completed single-node and first
# multi-node 1a0r/rb0.5 runs.
#
#SBATCH --job-name=packstar_overlap
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

SINGLE=/usr/xtmp/lz280/packstar_production/formal_singlenode_gpu/1a0r_rb0.5_12279812/rank_0/1a0r_packstar.csv
MULTI0=/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_0/1a0r_packstar.csv
MULTI1=/usr/xtmp/lz280/packstar_production/formal_multinode_gpu/1a0r_rb0.5_12279803/rank_1/1a0r_packstar.csv
REPORT=/usr/xtmp/lz280/packstar_production/interval_overlap_1a0r_rb0.5.tsv

awk -F, -v single="$SINGLE" -v multi0="$MULTI0" -v multi1="$MULTI1" '
function finite(x) { return x ~ /^[-+]?[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ }
function max(a,b) { return a > b ? a : b }
function min(a,b) { return a < b ? a : b }
function add(path, k) {
    k=$2
    if (path == single) { sLo[k]=$6; sHi[k]=$7; sSeen[k]=1 }
    else if (!(k in mSeen)) { mLo[k]=$6; mHi[k]=$7; mSeen[k]=1 }
}
FNR == 1 { next }
{ add(FILENAME) }
END {
    print "sequence\tsingleLB\tsingleUB\tmultiLB\tmultiUB\tintersection\tunion\tjaccard\tsingleCoverage\tmultiCoverage"
    n=0; sumI=0; sumJ=0; sumSC=0; sumMC=0; minI=1e300; maxI=-1; minJ=1e300; maxJ=-1
    for (k in sSeen) {
        if (!(k in mSeen)) continue
        if (!finite(sLo[k]) || !finite(sHi[k]) || !finite(mLo[k]) || !finite(mHi[k])) continue
        lo=max(sLo[k],mLo[k]); hi=min(sHi[k],mHi[k]); inter=max(0,hi-lo)
        uni=max(sHi[k],mHi[k])-min(sLo[k],mLo[k])
        sw=sHi[k]-sLo[k]; mw=mHi[k]-mLo[k]
        j=(uni>0 ? inter/uni : 1); sc=(sw>0 ? inter/sw : 1); mc=(mw>0 ? inter/mw : 1)
        printf "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", k,sLo[k],sHi[k],mLo[k],mHi[k],inter,uni,j,sc,mc
        n++; sumI+=inter; sumJ+=j; sumSC+=sc; sumMC+=mc
        if (inter<minI) minI=inter; if (inter>maxI) maxI=inter
        if (j<minJ) minJ=j; if (j>maxJ) maxJ=j
    }
    print "metric\tvalue"
    print "finiteIntervalPairs\t" n
    print "allFinitePairsOverlap\t" (n>0 && minI>=0 ? "true" : "false")
    print "intersectionMinLog10\t" (n>0 ? minI : "NA")
    print "intersectionMeanLog10\t" (n>0 ? sumI/n : "NA")
    print "intersectionMaxLog10\t" (n>0 ? maxI : "NA")
    print "jaccardMin\t" (n>0 ? minJ : "NA")
    print "jaccardMean\t" (n>0 ? sumJ/n : "NA")
    print "singleCoverageMean\t" (n>0 ? sumSC/n : "NA")
    print "multiCoverageMean\t" (n>0 ? sumMC/n : "NA")
}
' "$SINGLE" "$MULTI0" "$MULTI1" > "$REPORT"

echo "interval overlap report=$REPORT"
cat "$REPORT"
