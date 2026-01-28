/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.markstar.framework;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeTools;

public class MARKStarBound implements PartitionFunction.WithConfDB {

    // PartialFixCache (Phase 4): Threshold for skipping minimization based on tight bounds
    // If (quickUpperBound - confCorrection) < SKIP_THRESHOLD, skip full minimization
    // Typical RT at 300K ≈ 0.6 kcal/mol, so 1.0 kcal/mol is a reasonable threshold
    private static final double PARTIALFIX_SKIP_THRESHOLD = 1.0;  // kcal/mol

    protected double targetEpsilon = 1;
    public boolean debug = false;
    public boolean profileOutput = false;
    private Status status = null;
    private Values values = null;

    // the number of full conformations minimized
    private int numConfsEnergied = 0;
    // max confs minimized, -1 means infinite.
    private int maxNumConfs = -1;

    protected int maxMinimizations = 1;


    // the number of full conformations scored OR energied
    private int numConfsScored = 0;

    protected int numInternalNodesProcessed = 0;

    private boolean printMinimizedConfs;
    private MARKStarProgress progress;
    public String stateName = String.format("%4f",Math.random());
    private int numPartialMinimizations;
    private ArrayList<Integer> minList;
    protected double internalTimeAverage;
    protected double leafTimeAverage;
    private double cleanupTime;
    private boolean nonZeroLower;
    protected static TaskExecutor loopTasks;

    private ConfDB confDB;
    private ConfDB.Key confDBKey;

    public void setCorrections(UpdatingEnergyMatrix cachedCorrections) {
        correctionMatrix = cachedCorrections;
    }

    // Overwrite the computeUpperBound and computeLowerBound methods
    public static class Values extends PartitionFunction.Values {

        public Values ()
        {
            pstar = MathTools.BigPositiveInfinity;
        }
        @Override
        public BigDecimal calcUpperBound() {
            return pstar;
        }

        @Override
        public BigDecimal calcLowerBound() {
            return qstar;
        }

        @Override
        public double getEffectiveEpsilon() {
            return MathTools.bigDivide(pstar.subtract(qstar), pstar, decimalPrecision).doubleValue();
        }
    }

    public void setRCs(RCs rcs) {
        RCs = rcs;
    }

    public void setReportProgress(boolean showPfuncProgress) {
        this.printMinimizedConfs = true;
    }

    @Override
    public void setConfListener(ConfListener val) {

    }

    @Override
    public void setStabilityThreshold(BigDecimal threshold) {
        stabilityThreshold = threshold;
    }

    public void setMaxNumConfs(int maxNumConfs) {
        this.maxNumConfs = maxNumConfs;
    }

    @Override
	public void setConfDB(ConfDB confDB, ConfDB.Key key) {
    	this.confDB = confDB;
    	this.confDBKey = key;
	}

	private ConfDB.ConfTable confTable() {
    	if (confDB == null || confDBKey == null) {
    		return null;
		}
    	return confDB.get(confDBKey);
	}

    @Override
    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        status = Status.Estimating;
        values = new Values();
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        targetEpsilon = epsilon;
        status = Status.Estimating;
        values = new Values();
        this.stabilityThreshold = stabilityThreshold;
    }


    @Override
    public Status getStatus() {
        return status;
    }

    @Override
    public PartitionFunction.Values getValues() {
        return values;
    }

    @Override
    public int getParallelism() {
        return 0;
    }

    @Override
    public int getNumConfsEvaluated() {
        return numConfsEnergied;
    }

    public int getNumConfsScored() {
        return numConfsScored;
    }

    private int workDone() {
        return numInternalNodesProcessed + numConfsEnergied + numConfsScored + numPartialMinimizations ;
    }

    @Override
    public void compute(int maxNumConfs) {
        debugPrint("Num conformations: "+rootNode.getConfSearchNode().getNumConformations());
        double lastEps = 1;

        int previousConfCount = workDone();

        if(!nonZeroLower) {
            runUntilNonZero();
            updateBound();
            // DEBUG OUTPUT DISABLED
            // System.out.println("========== MARK* After runUntilNonZero ==========");
            // System.out.println("  epsilonBound:     " + epsilonBound);
            // System.out.println("  targetEpsilon:    " + targetEpsilon);
            // System.out.println("  lower bound:      " + rootNode.getLowerBound());
            // System.out.println("  upper bound:      " + rootNode.getUpperBound());
            // System.out.println("  Will enter while loop? " + (epsilonBound > targetEpsilon));
            // System.out.println("==================================================");
        }
        while (epsilonBound > targetEpsilon &&
                workDone()-previousConfCount < maxNumConfs
                && isStable(stabilityThreshold)) {
            debugPrint("Tightening from epsilon of "+epsilonBound);
            if(debug)
                debugHeap(queue);
            tightenBoundInPhases();
            debugPrint("Errorbound is now "+epsilonBound);
            if(lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println("Error. Bounds got looser.");
                //System.exit(-1);
            }
            lastEps = epsilonBound;
        }
        if(!isStable(stabilityThreshold))
            status = Status.Unstable;
        loopTasks.waitForFinish();
        minimizingEcalc.tasks.waitForFinish();
        BigDecimal averageReduction = BigDecimal.ZERO;
        int totalMinimizations = numConfsEnergied + numPartialMinimizations;
        if(totalMinimizations> 0)
            averageReduction = cumulativeZCorrection
                .divide(new BigDecimal(totalMinimizations), new MathContext(BigDecimal.ROUND_HALF_UP));
        debugPrint(String.format("Average Z reduction per minimization: %12.6e",averageReduction));
        values.pstar = rootNode.getUpperBound();
        values.qstar = rootNode.getLowerBound();
        // qprime should be the gap (upper - lower), not the upper bound itself
        // This matches the implementation in GradientDescentPfunc.java:486-489
        values.qprime = rootNode.getUpperBound().subtract(rootNode.getLowerBound());
        if(epsilonBound < targetEpsilon) {
            status = Status.Estimated;
            if(values.qstar.compareTo(BigDecimal.ZERO) == 0) {
                status = Status.Unstable;
            }
            //rootNode.printTree(stateName, minimizingEcalc.confSpace);
        }
    }

    protected void debugPrint(String s) {
        if(debug)
            System.out.println(s);
    }

    protected void profilePrint(String s) {
        if(profileOutput)
            System.out.println(s);
    }

    public void compute() {
        compute(Integer.MAX_VALUE);
    }

    @Override
    public Result makeResult() {
        // Calculate the upper bound z reductions from conf lower bounds, since we don't explicitly record these
        lowerReduction_ConfUpperBound = rootNode.getLowerBound().subtract(startLowerBound).subtract(lowerReduction_FullMin);
        // Calculate the lower bound z reductions from conf upper bounds, since we don't explicitly record these
        upperReduction_ConfLowerBound = startUpperBound.subtract(rootNode.getUpperBound()).subtract(upperReduction_FullMin).subtract(upperReduction_PartialMin);

        PartitionFunction.Result result = new PartitionFunction.Result(getStatus(), getValues(), getNumConfsEvaluated());
        /*
        result.setWorkInfo(numPartialMinimizations, numConfsScored,minList);
        result.setZInfo(lowerReduction_FullMin, lowerReduction_ConfUpperBound, upperReduction_FullMin, upperReduction_PartialMin, upperReduction_ConfLowerBound);
        result.setOrigBounds(startUpperBound, startLowerBound);
        result.setTimeInfo(stopwatch.getTimeNs());
        result.setMiscInfo(new BigDecimal(rootNode.getNumConfs()));
        */
        return result;
    }


    /**
     * TODO: 1. Make MARKStarBounds use and update a queue.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */
    // We keep track of the root node for computing our K* bounds
    protected MARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    protected final Queue<MARKStarNode> queue;
    protected double epsilonBound = Double.POSITIVE_INFINITY;
    private ConfIndex confIndex;
    public final AStarOrder order;
    // TODO: Implement new AStarPruner for MARK*?
    public final AStarPruner pruner;
    protected RCs RCs;
    protected Parallelism parallelism;
    private ObjectPool<ScoreContext> contexts;
    private MARKStarNode.ScorerFactory gscorerFactory;
    private MARKStarNode.ScorerFactory hscorerFactory;

    public boolean reduceMinimizations = true;
    private ConfAnalyzer confAnalyzer;
    EnergyMatrix minimizingEmat;
    EnergyMatrix rigidEmat;
    UpdatingEnergyMatrix correctionMatrix;
    ConfEnergyCalculator minimizingEcalc;
    private edu.duke.cs.osprey.ematrix.SubtreeDOFCache tripleDOFCache; // For storing triple DOF values
    private edu.duke.cs.osprey.ematrix.PartialFixCache partialFixCache; // PartialFixCache (Phase 4): L-set/M-set caching
    private Stopwatch stopwatch = new Stopwatch().start();
    // Variables for reporting pfunc reductions more accurately
    BigDecimal startUpperBound = null; //can't start with infinity
    BigDecimal startLowerBound = BigDecimal.ZERO;
    BigDecimal lowerReduction_FullMin = BigDecimal.ZERO; //Pfunc lower bound improvement from full minimization
    BigDecimal lowerReduction_ConfUpperBound = BigDecimal.ZERO; //Pfunc lower bound improvement from conf upper bounds
    BigDecimal upperReduction_FullMin = BigDecimal.ZERO; //Pfunc upper bound improvement from full minimization
    BigDecimal upperReduction_PartialMin = BigDecimal.ZERO; //Pfunc upper bound improvement from partial minimization corrections
    BigDecimal upperReduction_ConfLowerBound = BigDecimal.ZERO; //Pfunc upper bound improvement from conf lower bounds

    BigDecimal cumulativeZCorrection = BigDecimal.ZERO;//Pfunc upper bound improvement from partial minimization corrections
    BigDecimal ZReductionFromMin = BigDecimal.ZERO;//Pfunc lower bound improvement from full minimization
    BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private boolean computedCorrections = false;
    private long loopPartialTime = 0;
    private Set<String> correctedTuples = Collections.synchronizedSet(new HashSet<>());
    private BigDecimal stabilityThreshold;
    private double leafTimeSum = 0;
    private double internalTimeSum = 0;
    private int numLeavesScored = 0;
    private int numInternalScored = 0;

    public static MARKStarBound makeFromConfSpaceInfo(BBKStar.ConfSpaceInfo info, RCs rcs) {
        throw new UnsupportedOperationException("MARK* is not yet integrated into BBK*. Coming soon!");
        /*
        ConfEnergyCalculator minimizingConfEcalc = info.confEcalcMinimized;
        return new MARKStarBound(info.confSpace, info.ematRigid, info.ematMinimized, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
        */
    }

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                         ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        this.queue = new PriorityQueue<>();
        this.minimizingEcalc = minimizingConfEcalc;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        MPLPPairwiseHScorer scorer = new MPLPPairwiseHScorer(updater, minimizingEmat, 1, 0.0001);
        hscorerFactory = (emats) -> new MPLPPairwiseHScorer(updater, emats, 1, 0.0001);//TraditionalPairwiseHScorer(emats, rcs);//

        rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
                gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat),
                gscorerFactory.make(rigidEmat),
                new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);
        confIndex = new ConfIndex(rcs.getNumPos());
        this.minimizingEmat = minimizingEmat;
        this.rigidEmat = rigidEmat;
        this.RCs = rcs;
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat),hscorerFactory.make(minimizingEmat));
        this.pruner = null;

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.gscorer = gscorerFactory.make(minimizingEmat);
            context.hscorer = hscorerFactory.make(minimizingEmat);
            context.rigidscorer = gscorerFactory.make(rigidEmat);
            /** These scoreres should match the scorers in the MARKStarNode root - they perform the same calculations**/
            context.negatedhscorer = new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs); //this is used for upper bounds, so we want it rigid
            context.ecalc = minimizingConfEcalc;
            return context;
        });

        progress = new MARKStarProgress(RCs.getNumPos());
        //confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        setParallelism(parallelism);
        updateBound();
        // Recording pfunc starting bounds
        this.startLowerBound = rootNode.getLowerBound();
        this.startUpperBound = rootNode.getUpperBound();
        this.minList = new ArrayList<Integer>(Collections.nCopies(rcs.getNumPos(),0));
    }

    protected static class ScoreContext {
        public ConfIndex index;
        public AStarScorer gscorer;
        public AStarScorer hscorer;
        public AStarScorer negatedhscorer;
        public AStarScorer rigidscorer;
        public ConfEnergyCalculator ecalc;
    }



    public void setTripleDOFCache(edu.duke.cs.osprey.ematrix.SubtreeDOFCache cache) {
        this.tripleDOFCache = cache;
    }

    public void setPartialFixCache(edu.duke.cs.osprey.ematrix.PartialFixCache cache) {
        this.partialFixCache = cache;
    }

    public void setParallelism(Parallelism val) {

        if (val == null) {
            val = Parallelism.makeCpu(1);
        }

        parallelism = val;
        //loopTasks = minimizingEcalc.tasks;
        if(loopTasks == null)
            loopTasks = parallelism.makeTaskExecutor(1000);
        contexts.allocate(parallelism.getParallelism());
    }

    private void debugEpsilon(double curEpsilon) {
        if(debug && curEpsilon < epsilonBound) {
            System.err.println("Epsilon just got bigger.");
        }
    }

    protected boolean shouldMinimize(Node node) {
        return node.getLevel() == RCs.getNumPos() && !node.isMinimized();
    }

    protected void recordCorrection(double lowerBound, double correction) {
        BigDecimal upper = bc.calc(lowerBound);
        BigDecimal corrected = bc.calc(lowerBound + correction);
        cumulativeZCorrection = cumulativeZCorrection.add(upper.subtract(corrected));
        upperReduction_PartialMin = upperReduction_PartialMin.add(upper.subtract(corrected));
    }
    private void recordReduction(double lowerBound, double upperBound, double energy) {
        BigDecimal lowerBoundWeight = bc.calc(lowerBound);
        BigDecimal upperBoundWeight = bc.calc(upperBound);
        BigDecimal energyWeight = bc.calc(energy);
        ZReductionFromMin = ZReductionFromMin.add(lowerBoundWeight.subtract(upperBoundWeight));
        upperReduction_FullMin = upperReduction_FullMin.add(lowerBoundWeight.subtract(energyWeight));
        lowerReduction_FullMin = lowerReduction_FullMin.add(energyWeight.subtract(upperBoundWeight));

    }

    private void debugBreakOnConf(int[] conf) {
        int[] confOfInterest = new int[]{4,5,8,18};
        if(conf.length != confOfInterest.length)
            return;
        boolean match = true;
        for(int i = 0; i < confOfInterest.length; i++) {
            if(conf[i] != confOfInterest[i]) {
                match = false;
                break;
            }
        }
        if(match)
            System.out.println("Matched "+SimpleConfSpace.formatConfRCs(conf));
    }

    // We want to process internal nodes without worrying about the bound too much until we have
    // a nonzero lower bound. We have to have a nonzero lower bound, so we have to have at least
    // one node with a negative conf upper bound.
    private void runUntilNonZero() {
        // System.out.println("Running until leaf is found...");
        double bestConfUpper = Double.POSITIVE_INFINITY;

        List<MARKStarNode> newNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        boundLowestBoundConfUnderNode(rootNode,newNodes);
        queue.addAll(newNodes);


        newNodes.clear();
        // System.out.println("Found a leaf!");
        nonZeroLower = true;
    }

    protected void tightenBoundInPhases() {
        // System.out.println(String.format("Current overall error bound: %12.10f, spread of [%12.6e, %12.6e]",epsilonBound, rootNode.getLowerBound(), rootNode.getUpperBound()));
        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        int numNodes = 0;
        Stopwatch loopWatch = new Stopwatch();
        loopWatch.start();
        Stopwatch internalTime = new Stopwatch();
        Stopwatch leafTime = new Stopwatch();
        double leafTimeSum = 0;
        double internalTimeSum = 0;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ,leafZ};
        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        updateBound();
        debugPrint(String.format("After corrections, bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
        internalZ = ZSums[0];// MathTools.bigDivide(ZSums[0], new BigDecimal(Math.max(1,internalTimeAverage*internalNodes.size())), PartitionFunction.decimalPrecision);
        leafZ = ZSums[1]; //MathTools.bigDivide(ZSums[1], new BigDecimal(Math.max(1,leafTimeAverage)), PartitionFunction.decimalPrecision);
        // System.out.println(String.format("Z Comparison: %12.6e, %12.6e", internalZ, leafZ));
        if(MathTools.isLessThan(internalZ, leafZ)) {
            numNodes = leafNodes.size();
            // System.out.println("Processing "+numNodes+" leaf nodes...");
            leafTime.reset();
            leafTime.start();
            for(MARKStarNode leafNode: leafNodes) {
                processFullConfNode(newNodes, leafNode, leafNode.getConfSearchNode());
                leafNode.markUpdated();
                debugPrint("Processing Node: " + leafNode.getConfSearchNode().toString());
            }
            loopTasks.waitForFinish();
            leafTime.stop();
            leafTimeAverage = leafTime.getTimeS();
            // System.out.println("Processed "+numNodes+" leaves in "+leafTimeAverage+" seconds.");
            queue.addAll(internalNodes);
            if(maxMinimizations < parallelism.numThreads)
                maxMinimizations++;
        }
        else {
            numNodes = internalNodes.size();
            // System.out.println("Processing "+numNodes+" internal nodes...");
            internalTime.reset();
            internalTime.start();
            for (MARKStarNode internalNode : internalNodes) {
                if(!MathTools.isGreaterThan(internalNode.getLowerBound(),BigDecimal.ONE) &&
                    MathTools.isGreaterThan(
                            MathTools.bigDivide(internalNode.getUpperBound(),rootNode.getUpperBound(),
                                    PartitionFunction.decimalPrecision),
                            new BigDecimal(1-targetEpsilon))
                ) {
                    loopTasks.submit(() -> {
                        boundLowestBoundConfUnderNode(internalNode, newNodes);
                        return null;
                    }, (ignored) -> {
                    });
                }
                else {
                    processPartialConfNode(newNodes, internalNode, internalNode.getConfSearchNode());
                }
                internalNode.markUpdated();
            }
            loopTasks.waitForFinish();
            internalTime.stop();
            internalTimeSum=internalTime.getTimeS();
            internalTimeAverage = internalTimeSum/Math.max(1,internalNodes.size());
            debugPrint("Internal node time :"+internalTimeSum+", average "+internalTimeAverage);
            queue.addAll(leafNodes);
            numInternalNodesProcessed+=internalNodes.size();
        }
        if (epsilonBound <= targetEpsilon)
            return;
        loopCleanup(newNodes, loopWatch, numNodes);
    }

    protected void debugHeap(Queue<MARKStarNode> queue) {
        int maxNodes = 10;
        System.out.println("Node heap:");
        List<MARKStarNode> nodes = new ArrayList<>();
        while(!queue.isEmpty() && nodes.size() < 10)
        {
            MARKStarNode node = queue.poll();
            System.out.println(node.getConfSearchNode());
            nodes.add(node);
        }
        queue.addAll(nodes);
    }


    boolean isStable(BigDecimal stabilityThreshold) {
        return numConfsEnergied <= 0 || stabilityThreshold == null
                || MathTools.isGreaterThanOrEqual(rootNode.getUpperBound(), stabilityThreshold);
    }


    protected void populateQueues(Queue<MARKStarNode> queue, List<MARKStarNode> internalNodes, List<MARKStarNode> leafNodes, BigDecimal internalZ,
                                BigDecimal leafZ, BigDecimal[] ZSums) {
        List<MARKStarNode> leftoverLeaves = new ArrayList<>();
        int maxNodes = 1000;
        if(leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int)Math.floor(0.1*leafTimeAverage/internalTimeAverage));

        // LOG: Priority queue state before popping
        if (numConfsEnergied <= 5) {  // Only log for first few minimizations
            System.out.println("[QUEUE_STATE] size=" + queue.size() + " internal_nodes=" + internalNodes.size() + " leaf_nodes=" + leafNodes.size());
            List<MARKStarNode> peek = new ArrayList<>();
            int count = 0;
            while (!queue.isEmpty() && count < 10) {
                MARKStarNode n = queue.poll();
                peek.add(n);
                Node nd = n.getConfSearchNode();
                if (nd.getLevel() == RCs.getNumPos()) {  // leaf node
                    System.out.println("  [QUEUE_LEAF] rank=" + (count+1)
                        + " conf=" + SimpleConfSpace.formatConfRCs(nd.assignments)
                        + " errorBound=" + String.format("%.6e", n.getErrorBound())
                        + " confLower=" + String.format("%.6f", nd.getConfLowerBound())
                        + " gscore=" + String.format("%.6f", nd.gscore));
                }
                count++;
            }
            queue.addAll(peek);
        }

        while(!queue.isEmpty() && internalNodes.size() < maxNodes){
            MARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);
            double correctgscore = correctionMatrix.confE(node.assignments);
            double hscore = node.getConfLowerBound() - node.gscore;
            double confCorrection = Math.min(correctgscore, node.rigidScore) + hscore;

            // PartialFixCache: Try to tighten upper bound (symmetric to triple correction for lower bound)
            double quickUpperBound = tryQuickUpperBound(node);

            if(!node.isMinimized() && node.getConfLowerBound() < confCorrection
                    && node.getConfLowerBound() - confCorrection > 1e-5) {
                if(confCorrection < node.getConfLowerBound()) {
                    System.out.println("huh!?");
                }
                recordCorrection(node.getConfLowerBound(), correctgscore - node.gscore);

                node.gscore = correctgscore;
                if (confCorrection > node.rigidScore) {
                    System.out.println("Overcorrected"+SimpleConfSpace.formatConfRCs(node.assignments)+": " + confCorrection + " > " + node.rigidScore);
                    node.gscore = node.rigidScore;
                    confCorrection = node.rigidScore + hscore;
                }
                // Use tightened upper bound from PartialFixCache
                node.setBoundsFromConfLowerAndUpper(confCorrection, quickUpperBound);
                curNode.markUpdated();
                leftoverLeaves.add(curNode);
                continue;
            }

            BigDecimal diff = curNode.getUpperBound().subtract(curNode.getLowerBound());
            if (node.getLevel() < RCs.getNumPos()) {
                internalNodes.add(curNode);
                internalZ = internalZ.add(diff);
            }
            else if(shouldMinimize(node) && !correctedNode(leftoverLeaves, curNode, node)) {
                // LOG: Why this node will be minimized
                if (numConfsEnergied <= 5) {
                    System.out.println("[WILL_MINIMIZE] conf=" + SimpleConfSpace.formatConfRCs(node.assignments)
                        + " confLower=" + String.format("%.6f", node.getConfLowerBound())
                        + " confUpper=" + String.format("%.6f", node.getConfUpperBound())
                        + " gscore=" + String.format("%.6f", node.gscore)
                        + " rigidScore=" + String.format("%.6f", node.rigidScore)
                        + " correctgscore=" + String.format("%.6f", correctgscore)
                        + " confCorrection=" + String.format("%.6f", confCorrection)
                        + " errorBound=" + String.format("%.6e", curNode.getErrorBound()));
                }
                if(leafNodes.size() < maxMinimizations) {
                    leafNodes.add(curNode);
                    leafZ = leafZ.add(diff);
                }
                else
                    leftoverLeaves.add(curNode);
            }

        }
        ZSums[0] = internalZ;
        ZSums[1] = leafZ;
        queue.addAll(leftoverLeaves);
    }

    protected void loopCleanup(List<MARKStarNode> newNodes, Stopwatch loopWatch, int numNodes) {
        for(MARKStarNode node: newNodes) {
            if(node != null)
                queue.add(node);
        }
        loopWatch.stop();
        double loopTime = loopWatch.getTimeS();
        profilePrint("Processed "+numNodes+" this loop, spawning "+newNodes.size()+" in "+loopTime+", "+stopwatch.getTime()+" so far");
        loopWatch.reset();
        loopWatch.start();
        processPreminimization(minimizingEcalc);
        profilePrint("Preminimization time : "+loopWatch.getTime(2));
        double curEpsilon = epsilonBound;
        //rootNode.updateConfBounds(new ConfIndex(RCs.getNumPos()), RCs, gscorer, hscorer);
        updateBound();
        loopWatch.stop();
        cleanupTime = loopWatch.getTimeS();
        //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);
        // DEBUG OUTPUT DISABLED
        // System.out.println(String.format("Loop complete. Bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
    }

    protected boolean correctedNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        assert(curNode != null && node != null);
        double confCorrection = correctionMatrix.confE(node.assignments);
        if((node.getLevel() == RCs.getNumPos() && node.getConfLowerBound()< confCorrection)
                || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);

            // PartialFixCache (Phase 4): Tighten upper bound (symmetric to triple correction)
            double quickUpperBound = tryQuickUpperBound(node);

            node.setBoundsFromConfLowerAndUpper(node.getConfLowerBound() - oldg + confCorrection, quickUpperBound);
            curNode.markUpdated();
            newNodes.add(curNode);
            return true;
        }
        return false;
    }

    private MARKStarNode drillDown(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        try(ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            ConfIndex confIndex = context.index;
            node.index(confIndex);
            // which pos to expand next?
            int nextPos = order.getNextPos(confIndex, RCs);
            assert (!confIndex.isDefined(nextPos));
            assert (confIndex.isUndefined(nextPos));

            // score child nodes with tasks (possibly in parallel)
            List<MARKStarNode> children = new ArrayList<>();
            double bestChildLower = Double.POSITIVE_INFINITY;
            MARKStarNode bestChild = null;
            for (int nextRc : RCs.get(nextPos)) {

                if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                    continue;
                }

                // if this child was pruned dynamically, then don't score it
                if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                    continue;
                }
                Stopwatch partialTime = new Stopwatch().start();
                Node child = node.assign(nextPos, nextRc);
                double confLowerBound = Double.POSITIVE_INFINITY;

                // score the child node differentially against the parent node
                if (child.getLevel() < RCs.getNumPos()) {
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double diff = confCorrection;
                    double rigiddiff = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double hdiff = context.hscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    child.gscore = diff;
                    //Correct for incorrect gscore.
                    rigiddiff = rigiddiff - node.gscore + node.rigidScore;
                    child.rigidScore = rigiddiff;

                    confLowerBound = child.gscore + hdiff;
                    double confUpperbound = rigiddiff + maxhdiff;
                    child.computeNumConformations(RCs);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }

                    // PartialFixCache (Phase 4): For now, use original upper bound
                    // TODO: Implement quick upper bound tightening via PartialFixCache
                    double quickUpperBound = confUpperbound;

                    child.setBoundsFromConfLowerAndUpper(confLowerBound, quickUpperBound);
                    progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                }
                if (child.getLevel() == RCs.getNumPos()) {
                    double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    confRigid=confRigid-node.gscore+node.rigidScore;

                    child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double lowerbound = minimizingEmat.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    checkBounds(confCorrection, confRigid);

                    // PartialFixCache (Phase 4): For now, use rigid upper bound for leaf
                    // TODO: Implement quick upper bound tightening via PartialFixCache
                    double quickUpperBound = confRigid;

                    child.setBoundsFromConfLowerAndUpper(confCorrection, quickUpperBound);
                    child.gscore = child.getConfLowerBound();
                    confLowerBound = lowerbound;
                    child.rigidScore = confRigid;

                    // LOG: Leaf node created with pairwise energy (DISABLED)
                    // if (numConfsScored < 10) {
                    //     System.out.println("[LEAF_CREATED] conf=" + SimpleConfSpace.formatConfRCs(child.assignments)
                    //         + " pairwise=" + String.format("%.6f", lowerbound)
                    //         + " correction=" + String.format("%.6f", confCorrection)
                    //         + " rigid=" + String.format("%.6f", confRigid)
                    //         + " gscore=" + String.format("%.6f", child.gscore));
                    // }

                    numConfsScored++;
                    progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                }
                partialTime.stop();
                loopPartialTime += partialTime.getTimeS();


                if (Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                MARKStarNodeChild.markUpdated();
                if (confLowerBound < bestChildLower) {
                    bestChild = MARKStarNodeChild;
                }
                // collect the possible children
                if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                    children.add(MARKStarNodeChild);
                }
                newNodes.add(MARKStarNodeChild);

            }
            return bestChild;
        }
    }

    protected void boundLowestBoundConfUnderNode(MARKStarNode startNode, List<MARKStarNode> generatedNodes) {
        Comparator<MARKStarNode> confBoundComparator = Comparator.comparingDouble(o -> o.getConfSearchNode().getConfLowerBound());
        PriorityQueue<MARKStarNode> drillQueue = new PriorityQueue<>(confBoundComparator);
        drillQueue.add(startNode);

        List<MARKStarNode> newNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        while(!drillQueue.isEmpty()) {
            numNodes++;
            MARKStarNode curNode = drillQueue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);

            if (node.getLevel() < RCs.getNumPos()) {
                MARKStarNode nextNode = drillDown(newNodes, curNode, node);
                newNodes.remove(nextNode);
                drillQueue.add(nextNode);
            }else if(RCs.getNumPos() == 0){
                // weird case where there are no residues?
                curNode.recomputeEpsilon();
            }
            else {
                newNodes.add(curNode);
            }

            //debugHeap(drillQueue, true);
            if(leafLoop.getTimeS() > 1) {
                leafLoop.stop();
                leafLoop.reset();
                leafLoop.start();
                System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",numNodes, overallLoop.getTime(2),rootNode.getLowerBound(),rootNode.getUpperBound()));
            }
        }
        generatedNodes.addAll(newNodes);

    }

    protected void processPartialConfNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        // which pos to expand next?
        node.index(confIndex);
        int nextPos = order.getNextPos(confIndex, RCs);
        assert (!confIndex.isDefined(nextPos));
        assert (confIndex.isUndefined(nextPos));

        // score child nodes with tasks (possibly in parallel)
        List<MARKStarNode> children = new ArrayList<>();
        for (int nextRc : RCs.get(nextPos)) {

            if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                continue;
            }

            // if this child was pruned dynamically, then don't score it
            if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                continue;
            }

            loopTasks.submit(() -> {

                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                    Stopwatch partialTime = new Stopwatch().start();
                    ScoreContext context = checkout.get();
                    node.index(context.index);
                    Node child = node.assign(nextPos, nextRc);

                    // score the child node differentially against the parent node
                    if (child.getLevel() < RCs.getNumPos()) {
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double diff = confCorrection;
                        double rigiddiff = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double hdiff = context.hscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        child.gscore = diff;
                        //Correct for incorrect gscore.
                        rigiddiff=rigiddiff-node.gscore+node.rigidScore;
                        child.rigidScore = rigiddiff;

                        double confLowerBound = child.gscore + hdiff;
                        double confUpperbound = rigiddiff + maxhdiff;
                        child.computeNumConformations(RCs);
                        double lowerbound = minimizingEmat.confE(child.assignments);
                        if(diff < confCorrection) {
                            recordCorrection(confLowerBound, confCorrection - diff);
                            confLowerBound = confCorrection + hdiff;
                        }

                        // PartialFixCache (Phase 4): For now, use original upper bound (parallel version)
                        // TODO: Implement quick upper bound tightening via PartialFixCache
                        double quickUpperBound = confUpperbound;

                        child.setBoundsFromConfLowerAndUpper(confLowerBound, quickUpperBound);
                        progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                    }
                    if (child.getLevel() == RCs.getNumPos()) {
                        double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        confRigid=confRigid-node.gscore+node.rigidScore;

                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double lowerbound = minimizingEmat.confE(child.assignments);

                        if(lowerbound < confCorrection) {
                            recordCorrection(lowerbound, confCorrection - lowerbound);
                        }
                        checkBounds(confCorrection,confRigid);

                        // PartialFixCache (Phase 4): For now, use rigid upper bound for leaf (parallel version)
                        // TODO: Implement quick upper bound tightening via PartialFixCache
                        double quickUpperBound = confRigid;

                        child.setBoundsFromConfLowerAndUpper(confCorrection, quickUpperBound);
                        child.gscore = confCorrection;
                        child.rigidScore = confRigid;
                        numConfsScored++;
                        progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                    }
                    partialTime.stop();
                    loopPartialTime+=partialTime.getTimeS();


                    return child;
                }

            }, (Node child) -> {
                if(Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                    // collect the possible children
                    if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                        children.add(MARKStarNodeChild);
                    }
                    if (!child.isMinimized()) {
                        newNodes.add(MARKStarNodeChild);
                    }
                    else
                        MARKStarNodeChild.computeEpsilonErrorBounds();

                curNode.markUpdated();
            });
        }
    }


    protected void processFullConfNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        double confCorrection = correctionMatrix.confE(node.assignments);
        double pairwiseLowerBound = node.getConfLowerBound();
        double currentGscore = node.gscore;
        double currentUpper = node.getConfUpperBound();

        // PartialFixCache (Phase 4): Calculate quick upper bound BEFORE making skip decision
        double quickUpperBound = tryQuickUpperBound(node);

        // Decision 1: Triple correction - skip if correction gives tighter lower bound
        boolean shouldSkipDueToLowerBound =
            (pairwiseLowerBound < confCorrection || currentGscore < confCorrection);

        // Decision 2: PartialFixCache - skip if quick optimization gives tighter upper bound (symmetric to triple correction)
        boolean shouldSkipDueToUpperBound =
            (partialFixCache != null && quickUpperBound < currentUpper);

        // Decision 3: PartialFixCache - skip if bounds gap is already tight enough (cost-benefit analysis)
        // TEMPORARILY DISABLED: First test cache functionality alone
        double boundsGap = quickUpperBound - confCorrection;
        boolean shouldSkipDueToTightGap = false;
            // (partialFixCache != null && boundsGap >= 0 && boundsGap < PARTIALFIX_SKIP_THRESHOLD);

        // LOG: Decision point - every time we consider minimizing a conformation
        boolean hasCorrection = confCorrection > Double.NEGATIVE_INFINITY;
        boolean willSkipMinimization = shouldSkipDueToLowerBound || shouldSkipDueToUpperBound || shouldSkipDueToTightGap;

        String skipReason = "no_sufficient_correction";
        if (shouldSkipDueToLowerBound && shouldSkipDueToUpperBound) {
            skipReason = "both_bounds_tightened";
        } else if (shouldSkipDueToLowerBound && shouldSkipDueToTightGap) {
            skipReason = "lower_tightened_and_tight_gap";
        } else if (shouldSkipDueToUpperBound && shouldSkipDueToTightGap) {
            skipReason = "upper_tightened_and_tight_gap";
        } else if (shouldSkipDueToLowerBound) {
            skipReason = "triple_correction_lower";
        } else if (shouldSkipDueToUpperBound) {
            skipReason = "partialfix_upper";
        } else if (shouldSkipDueToTightGap) {
            skipReason = "tight_gap";
        }

        double upperImprovement = currentUpper - quickUpperBound;
        double lowerImprovement = confCorrection - pairwiseLowerBound;

        // System.out.println("[MINIMIZATION_DECISION] conf=" + SimpleConfSpace.formatConfRCs(node.assignments)
        //     + " pairwise=" + String.format("%.6f", pairwiseLowerBound)
        //     + " gscore=" + String.format("%.6f", currentGscore)
        //     + " correction=" + (hasCorrection ? String.format("%.6f", confCorrection) : "NONE")
        //     + " lowerImprove=" + String.format("%.6f", lowerImprovement)
        //     + " currentUpper=" + String.format("%.6f", currentUpper)
        //     + " quickUpper=" + String.format("%.6f", quickUpperBound)
        //     + " upperImprove=" + String.format("%.6f", upperImprovement)
        //     + " boundsGap=" + String.format("%.6f", boundsGap)
        //     + " decision=" + (willSkipMinimization ? "SKIP" : "MINIMIZE")
        //     + " reason=" + skipReason);

        if(shouldSkipDueToLowerBound || shouldSkipDueToUpperBound || shouldSkipDueToTightGap) {
            double oldg = node.gscore;
            node.gscore = confCorrection;

            if (shouldSkipDueToLowerBound) {
                recordCorrection(oldg, confCorrection - oldg);
            }

            // PartialFixCache (Phase 4): Use quick upper bound (already calculated above)
            node.setBoundsFromConfLowerAndUpper(confCorrection, quickUpperBound);
            curNode.markUpdated();
            newNodes.add(curNode);

            // LOG: Minimization skipped
            double finalGap = quickUpperBound - confCorrection;
            System.out.println("[SKIP_MINIMIZATION] conf=" + SimpleConfSpace.formatConfRCs(node.assignments)
                + " reason=" + skipReason
                + " lower: " + String.format("%.6f", pairwiseLowerBound) + " -> " + String.format("%.6f", confCorrection)
                + " (+" + String.format("%.6f", confCorrection - pairwiseLowerBound) + ")"
                + " upper: " + String.format("%.6f", currentUpper) + " -> " + String.format("%.6f", quickUpperBound)
                + " (-" + String.format("%.6f", currentUpper - quickUpperBound) + ")"
                + " finalGap=" + String.format("%.6f", finalGap));
            return;
        }
        loopTasks.submit(() -> {
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                node.index(context.index);

                ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, node.getConfLowerBound());
                Stopwatch minimizationTimer = new Stopwatch().start();
                ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
                minimizationTimer.stop();

                // LOG: Detailed minimization info
                // double pairwiseEstimate = node.getConfLowerBound();
                // double minimizedEnergy = analysis.epmol.energy;
                // double energyGap = minimizedEnergy - pairwiseEstimate;
                // System.out.println("[MINIMIZE] conf=" + SimpleConfSpace.formatConfRCs(node.assignments)
                //     + " pairwise=" + String.format("%.6f", pairwiseEstimate)
                //     + " minimized=" + String.format("%.6f", minimizedEnergy)
                //     + " gap=" + String.format("%.6f", energyGap)
                //     + " time=" + String.format("%.2f", minimizationTimer.getTimeMs()) + "ms");
                
                // record the conf energy in the ConfDB, if needed
                ConfDB.ConfTable confTable = confTable();
                if (confTable != null) {
                	long timestamp = TimeTools.getTimestampNs();
                	confTable.setLowerBound(conf.getAssignments(), conf.getScore(), timestamp);
                    confTable.setUpperBound(conf.getAssignments(), analysis.epmol.energy, timestamp);
                }
                
                Stopwatch correctionTimer = new Stopwatch().start();
                computeEnergyCorrection(analysis, conf, context.ecalc);

                double energy = analysis.epmol.energy;
                double newConfUpper = energy;
                double newConfLower = energy;
                // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
                double oldConfUpper = node.getConfUpperBound();
                double oldConfLower = node.getConfLowerBound();
                checkConfLowerBound(node, energy);
                if (newConfUpper > oldConfUpper) {
                    System.err.println("Upper bounds got worse after minimization:" + newConfUpper
                            + " > " + (oldConfUpper)+". Rejecting minimized energy.");
                    System.err.println("Node info: "+node);

                    newConfUpper = oldConfUpper;
                    newConfLower = oldConfUpper;
                }

                // RESTORED: Use minimized energy directly (original logic)
                // Do NOT re-apply correction after minimization
                curNode.setBoundsFromConfLowerAndUpper(newConfLower, newConfUpper);
                double oldgscore = node.gscore;
                node.gscore = newConfLower;

                // LOG: Final bounds after minimization
                // System.out.println("[MINIMIZED_FINAL] conf=" + SimpleConfSpace.formatConfRCs(node.assignments)
                //     + " minimized=" + String.format("%.6f", energy)
                //     + " finalLower=" + String.format("%.6f", newConfLower)
                //     + " finalUpper=" + String.format("%.6f", newConfUpper)
                //     + " oldGscore=" + String.format("%.6f", oldgscore));

                String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (node.getConfLowerBound()) + "," + (node.getConfUpperBound()) + "]";
                debugPrint(out);
                curNode.markUpdated();
                synchronized(this) {
                    numConfsEnergied++;
                    minList.set(conf.getAssignments().length-1,minList.get(conf.getAssignments().length-1)+1);
                    recordReduction(oldConfLower, oldConfUpper, energy);
                    printMinimizationOutput(node, newConfLower, oldgscore);

                    // PRIORITY QUEUE DUMP: After first minimization
                    if (numConfsEnergied == 1) {
                        dumpPriorityQueueAfterFirstMinimization(node, energy, newConfLower, newConfUpper);
                    }
                }


            }
            return null;
        },
                // Dummy function. We're not doing anything here.
                (Node child) -> {
                    progress.reportLeafNode(node.gscore, queue.size(), epsilonBound);
                    if(!node.isMinimized())
                        newNodes.add(curNode);

                });
    }

    private void printMinimizationOutput(Node node, double newConfLower, double oldgscore) {
        // if (printMinimizedConfs) {
        //     System.out.println("[" + SimpleConfSpace.formatConfRCs(node.assignments) + "]"
        //             + String.format("conf:%4d, score:%12.6f, lower:%12.6f, corrected:%12.6f energy:%12.6f"
        //                     +", bounds:[%12e, %12e], delta:%12.6f, time:%10s",
        //             numConfsEnergied, oldgscore, minimizingEmat.confE(node.assignments),
        //             correctionMatrix.confE(node.assignments), newConfLower,
        //             rootNode.getConfSearchNode().getSubtreeLowerBound(),rootNode.getConfSearchNode().getSubtreeUpperBound(),
        //             epsilonBound, stopwatch.getTime(2)));

        // }
    }

    /**
     * Dump detailed priority queue state after first minimization
     * This shows what nodes are in the queue and their priorities
     * DISABLED: Too verbose for production use
     */
    private void dumpPriorityQueueAfterFirstMinimization(Node minimizedNode, double energy, double newConfLower, double newConfUpper) {
        // DEBUG OUTPUT DISABLED
        /*
        System.out.println("\n" + "=".repeat(100));
        System.out.println("[PRIORITY_QUEUE_DUMP] After minimizing FIRST conformation");
        System.out.println("=".repeat(100));

        // Information about the minimized node
        System.out.println("\n[MINIMIZED_CONF]");
        System.out.println("  Conf: " + SimpleConfSpace.formatConfRCs(minimizedNode.assignments));
        System.out.println("  Minimized energy: " + String.format("%.6f", energy));
        System.out.println("  New lower bound: " + String.format("%.6f", newConfLower));
        System.out.println("  New upper bound: " + String.format("%.6f", newConfUpper));
        System.out.println("  gscore: " + String.format("%.6f", minimizedNode.gscore));

        // Dump the entire priority queue (up to 100 nodes)
        System.out.println("\n[QUEUE_CONTENTS] Total size: " + queue.size());
        System.out.println("Format: Rank | Level | ErrorBound | ConfLower | ConfUpper | gscore | rigidScore | Assignments");
        System.out.println("-".repeat(100));

        List<MARKStarNode> tempList = new ArrayList<>();
        int rank = 1;
        int maxToDump = Math.min(100, queue.size());

        while (!queue.isEmpty() && rank <= maxToDump) {
            MARKStarNode mNode = queue.poll();
            tempList.add(mNode);
            Node node = mNode.getConfSearchNode();

            boolean isLeaf = (node.getLevel() == RCs.getNumPos());
            String nodeType = isLeaf ? "LEAF" : "INTERNAL";

            System.out.println(String.format("  %3d | %-8s | %12.6e | %10.6f | %10.6f | %10.6f | %10.6f | %s",
                rank,
                nodeType,
                mNode.getErrorBound(),
                node.getConfLowerBound(),
                node.getConfUpperBound(),
                node.gscore,
                node.rigidScore,
                SimpleConfSpace.formatConfRCs(node.assignments)));

            rank++;
        }

        // Restore the queue
        queue.addAll(tempList);

        System.out.println("-".repeat(100));
        System.out.println("[QUEUE_DUMP_END] Dumped " + (rank-1) + " nodes out of " + queue.size() + " total");
        System.out.println("=".repeat(100) + "\n");
        */
    }

    private void checkConfLowerBound(Node node, double energy) {
        if(energy < node.getConfLowerBound()) {
            System.err.println("Bounds are incorrect:" + (node.getConfLowerBound()) + " > "
                    + energy);
            if (energy < 10)
                System.err.println("The bounds are probably wrong.");
            //System.exit(-1);
        }
    }


    private void checkBounds(double lower, double upper)
    {
        if (upper < lower && upper - lower > 1e-5 && upper < 10)
            debugPrint("Bounds incorrect.");
    }

    private void computeEnergyCorrection(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf,
                                                  ConfEnergyCalculator ecalc) {
        if(conf.getAssignments().length < 3)
            return;

        // Track corrections generated from this conformation
        int triplesGenerated = 0;

        //System.out.println("Analysis:"+analysis);
        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(minimizingEmat);
        Stopwatch correctionTime = new Stopwatch().start();
        //System.out.println("Energy Analysis: "+energyAnalysis);
        //System.out.println("Score Analysis: "+scoreAnalysis);
        EnergyMatrix diff = energyAnalysis.diff(scoreAnalysis);
        //System.out.println("Difference Analysis " + diff);
        List<TupE> sortedPairwiseTerms2 = new ArrayList<>();
        for (int pos = 0; pos < diff.getNumPos(); pos++)
        {
            for (int rc = 0; rc < diff.getNumConfAtPos(pos); rc++)
            {
                for (int pos2 = 0; pos2 < diff.getNumPos(); pos2++)
                {
                    for (int rc2 = 0; rc2 < diff.getNumConfAtPos(pos2); rc2++)
                    {
                        if(pos >= pos2)
                            continue;
                        double sum = 0;
                        sum+=diff.getOneBody(pos, rc);
                        sum+=diff.getPairwise(pos, rc, pos2, rc2);
                        sum+=diff.getOneBody(pos2,rc2);
                        TupE tupe = new TupE(new RCTuple(pos, rc, pos2, rc2), sum);
                        sortedPairwiseTerms2.add(tupe);
                    }
                }
            }
        }
        Collections.sort(sortedPairwiseTerms2);

        double threshhold = 0.1;
        double minDifference = 0.9;
        double triplethreshhold = 0.3;
        double maxDiff = sortedPairwiseTerms2.get(0).E;
        for(int i = 0; i < sortedPairwiseTerms2.size(); i++)
        {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if(pairDiff < minDifference &&  maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, tupe.E);
            int pos1 = tupe.tup.pos.get(0);
            int pos2 = tupe.tup.pos.get(1);
            int localMinimizations = 0;
            for(int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                double tupleBounds = rigidEmat.getInternalEnergy(tuple) - minimizingEmat.getInternalEnergy(tuple);
                if(tupleBounds < triplethreshhold)
                    continue;
                minList.set(tuple.size()-1,minList.get(tuple.size()-1)+1);
                computeDifference(tuple, minimizingEcalc);
                localMinimizations++;
                triplesGenerated++;
            }
            numPartialMinimizations+=localMinimizations;
            progress.reportPartialMinimization(localMinimizations, epsilonBound);
        }
        correctionTime.stop();
        ecalc.tasks.waitForFinish();

        // LOG: Summary of corrections generated
        // System.out.println("[CORRECTION_SUMMARY] conf=" + SimpleConfSpace.formatConfRCs(conf.getAssignments())
        //     + " triples_generated=" + triplesGenerated
        //     + " computation_time=" + String.format("%.2f", correctionTime.getTimeMs()) + "ms");
    }




    private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        computedCorrections = true;
        if(correctedTuples.contains(tuple.stringListing()))
            return;
        correctedTuples.add(tuple.stringListing());
        if(correctionMatrix.hasHigherOrderTermFor(tuple))
            return;
        minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
            double tripleEnergy = minimizedTuple.energy;

            double lowerbound = minimizingEmat.getInternalEnergy(tuple);
            double correction = tripleEnergy - lowerbound;

            // LOG: Triple correction details (DISABLED)
            // System.out.println("[TRIPLE_CORRECTION] tuple=" + tuple.stringListing()
            //     + " pairwise=" + String.format("%.6f", lowerbound)
            //     + " minimized=" + String.format("%.6f", tripleEnergy)
            //     + " correction=" + String.format("%.6f", correction)
            //     + " negative=" + (correction < 0));

            if (correction > 0) {
                correctionMatrix.setHigherOrder(tuple, correction);
            }
            else {
                // When triple energy is lower than pairwise bound, the matrix bound is too high.
                // This can happen when pairwise minimization doesn't find the global minimum
                // but triple minimization does. We should use the lower (more accurate) triple energy.
                System.err.println("Negative correction for "+tuple.stringListing() +
                                 " (correction=" + correction + "). Using triple energy as correction.");
                // Store a negative correction to lower the bound
                correctionMatrix.setHigherOrder(tuple, correction);
            }

            // NEW: Store triple DOF values in cache
            if (edu.duke.cs.osprey.ematrix.SubtreeDOFCache.ENABLE_TRIPLE_DOF_CACHE &&
                tripleDOFCache != null &&
                tuple.size() == 3 &&
                minimizedTuple.params != null) {

                try {
                    tripleDOFCache.storeTripleDOFs(
                        tuple,
                        minimizedTuple.params,
                        tripleEnergy,
                        minimizedTuple.pmol
                    );
                } catch (Exception e) {
                    System.err.println("[TRIPLE_DOF_CACHE_ERROR] " + e.getMessage());
                }
            }
        });
    }

    private RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for(int pos: positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }

    private void processPreminimization(ConfEnergyCalculator ecalc) {
        int maxMinimizations = 1;//parallelism.numThreads;
        List<MARKStarNode> topConfs = getTopConfs(maxMinimizations);
        // Need at least two confs to do any partial preminimization
        if (topConfs.size() < 2) {
            queue.addAll(topConfs);
            return;
        }
        RCTuple lowestBoundTuple = topConfs.get(0).toTuple();
        RCTuple overlap = findLargestOverlap(lowestBoundTuple, topConfs, 3);
        //Only continue if we have something to minimize
        for (MARKStarNode conf : topConfs) {
            RCTuple confTuple = conf.toTuple();
            if(minimizingEmat.getInternalEnergy(confTuple) == rigidEmat.getInternalEnergy(confTuple))
                continue;
            numPartialMinimizations++;
            minList.set(confTuple.size()-1,minList.get(confTuple.size()-1)+1);
            if (confTuple.size() > 2 && confTuple.size() < RCs.getNumPos ()){
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(minimizingEcalc, conf.toTuple());
                    return null;
                }, (econf) -> {
                });
            }
        }
        //minimizingEcalc.tasks.waitForFinish();
        ConfIndex index = new ConfIndex(RCs.getNumPos());
        if(overlap.size() > 3 && !correctionMatrix.hasHigherOrderTermFor(overlap)
                && minimizingEmat.getInternalEnergy(overlap) != rigidEmat.getInternalEnergy(overlap)) {
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(ecalc, overlap);
                    return null;
                }, (econf) -> {
                });
        }
        queue.addAll(topConfs);
    }

    private void computeTupleCorrection(ConfEnergyCalculator ecalc, RCTuple overlap) {
        if(correctionMatrix.hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = minimizingEmat.getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
        progress.reportPartialMinimization(1, epsilonBound);
        if(partiallyMinimizedLower > pairwiseLower)
        synchronized (correctionMatrix) {
            correctionMatrix.setHigherOrder(overlap, partiallyMinimizedLower - pairwiseLower);
        }
        progress.reportPartialMinimization(1, epsilonBound);
    }

    private List<MARKStarNode> getTopConfs(int numConfs) {
        List<MARKStarNode> topConfs = new ArrayList<>();
        while (topConfs.size() < numConfs&& !queue.isEmpty()) {
            MARKStarNode nextLowestConf = queue.poll();
            topConfs.add(nextLowestConf);
        }
        return topConfs;
    }


    private RCTuple findLargestOverlap(RCTuple conf, List<MARKStarNode> otherConfs, int minResidues) {
        RCTuple overlap = conf;
        for(MARKStarNode other: otherConfs) {
            overlap = overlap.intersect(other.toTuple());
            if(overlap.size() < minResidues)
                break;
        }
        return overlap;

    }

    protected void updateBound() {
        double curEpsilon = epsilonBound;
        Stopwatch time = new Stopwatch().start();
        epsilonBound = rootNode.computeEpsilonErrorBounds();
        time.stop();
        //System.out.println("Bound update time: "+time.getTime(2));
        debugEpsilon(curEpsilon);
        //System.out.println("Current epsilon:"+epsilonBound);
    }

    private boolean hasPrunedPair(ConfIndex confIndex, int nextPos, int nextRc) {

        // do we even have pruned pairs?
        PruningMatrix pmat = RCs.getPruneMat();
        if (pmat == null) {
            return false;
        }

        for (int i = 0; i < confIndex.numDefined; i++) {
            int pos = confIndex.definedPos[i];
            int rc = confIndex.definedRCs[i];
            assert (pos != nextPos || rc != nextRc);
            if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Try to tighten upper bound using PartialFixCache without full minimization
     * Returns the tightened upper bound, or original if PartialFixCache not available
     *
     * This method performs quick optimization (5 CCD iterations on M-set only) using
     * PartialFixCache to get a tighter upper bound without full minimization.
     */
    protected double tryQuickUpperBound(Node node) {
        // DEBUG OUTPUT DISABLED
        // System.out.println("[MARKStarBound DEBUG] tryQuickUpperBound() called, partialFixCache: " +
        //     (partialFixCache != null ? "EXISTS" : "NULL"));

        if (partialFixCache == null) {
            return node.getConfUpperBound();
        }

        // DEBUG OUTPUT DISABLED
        // System.out.println("[MARKStarBound DEBUG] Calling partialFixCache.minimizeWithPartialFixCache()");

        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();

            // Step 1: Create RCTuple from node assignments
            RCTuple conf = new RCTuple(node.assignments);

            // Step 2: Make the parametric molecule for this conformation
            ParametricMolecule pmol = context.ecalc.confSpace.makeMolecule(conf);

            // Step 3: Get initial DOFs (centered in voxel)
            cern.colt.matrix.DoubleMatrix1D initialDOFs =
                cern.colt.matrix.DoubleFactory1D.dense.make(pmol.dofs.size());
            pmol.dofBounds.getCenter(initialDOFs);

            // Step 4: Create residue interactions
            ResidueInteractions inters = context.ecalc.makeFragInters(conf);

            // Step 5: Access EnergyCalculator's context to get minimizer factory
            // The context is directly accessible from EnergyCalculator
            edu.duke.cs.osprey.energy.EnergyCalculator.Type.Context ecalcContext =
                context.ecalc.ecalc.context;

            // Step 6: Create energy function
            edu.duke.cs.osprey.energy.EnergyFunction efunc =
                ecalcContext.efuncs.make(inters, pmol.mol);

            // Step 7: Create objective function
            edu.duke.cs.osprey.minimization.ObjectiveFunction objectiveFunction =
                new edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction(pmol, efunc);

            // Step 8: Create minimizer
            edu.duke.cs.osprey.minimization.Minimizer minimizer =
                ecalcContext.minimizers.make(objectiveFunction);

            // Step 9: Call PartialFixCache to perform quick minimization
            edu.duke.cs.osprey.ematrix.PartialFixCache.MinimizationResult result =
                partialFixCache.minimizeWithPartialFixCache(
                    conf,
                    minimizer,
                    initialDOFs,
                    objectiveFunction,
                    pmol
                );

            // Step 10: Clean up resources
            efunc.close();
            minimizer.close();

            // Return the minimized energy as the tightened upper bound
            return result.energy;

        } catch (Exception e) {
            // If anything goes wrong, fall back to original upper bound
            // This ensures robustness - we never fail, just miss optimization opportunity
            if (debug) {
                System.err.println("WARNING: PartialFixCache tryQuickUpperBound failed for conf "
                    + SimpleConfSpace.formatConfRCs(node.assignments) + ": " + e.getMessage());
                e.printStackTrace();
            }
            return node.getConfUpperBound();
        }
    }
}
