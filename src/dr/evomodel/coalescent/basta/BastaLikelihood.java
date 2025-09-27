/*
 * BastaLikelihood.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.*;
import dr.evolution.tree.*;
import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TipStateAccessor;
import dr.evomodel.treelikelihood.AncestralStateTraitProvider;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.*;
import java.util.function.Function;
import java.util.logging.Logger;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.*;

/**
 * @author Guy Baele
 * @author Yucai Shao
 * @author Marc A. Suchard
 */

public class BastaLikelihood extends AbstractModelLikelihood implements
        TreeTraitProvider, AncestralStateTraitProvider, Citable, Profileable, Reportable, TipStateAccessor {

    private static final boolean COUNT_TOTAL_OPERATIONS = true;

    private final BastaLikelihoodDelegate likelihoodDelegate;

    private final Tree tree;
    private final PatternList patternList;
    private final SubstitutionModel substitutionModel;
    private final Parameter popSizeParameter;
    private final Parameter growthRateParameter;
    private final BranchRateModel branchRateModel;
    private final int stateCount;

    private final Helper treeTraits = new Helper();

    private final CoalescentIntervalTraversal treeTraversalDelegate;
    public final BestSignalsFromBigFastTreeIntervals treeIntervals;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown;

    private boolean populationSizesKnown;
    private boolean treeIntervalsKnown;
    private boolean transitionMatricesKnown;

    // Ancestral state reconstruction variables
    private final DataType dataType;
    private final String tag;
    private final CodeFormatter formatter;

    // State reconstruction settings
    private final boolean useMAP;
    private final boolean returnMarginalLogLikelihood;
    private final boolean conditionalProbabilitiesInLogSpace;

    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;
    private int[][] subIntervalStates;  // [interval][nodeNumber][pattern]
    private Map<Integer, List<Integer>> nodeIntervalMap = new HashMap<>();
    protected boolean areStatesRedrawn = false;
    protected boolean storedAreStatesRedrawn = false;
    protected double jointLogLikelihood;
    private double storedJointLogLikelihood;

    public BastaLikelihood(String name,
                           Tree treeModel,
                           PatternList patternList,
                           SubstitutionModel substitutionModel,
                           Parameter popSizeParameter,
                           Parameter growthRateParameter,
                           BranchRateModel branchRateModel,
                           BastaLikelihoodDelegate likelihoodDelegate,
                           int numberSubIntervals,
                           boolean useAmbiguities) {
        this(name, treeModel, patternList, substitutionModel, popSizeParameter,
                growthRateParameter, branchRateModel, likelihoodDelegate,
                numberSubIntervals, useAmbiguities,
                substitutionModel.getDataType(), "states", false, true, false);
    }

    public BastaLikelihood(String name,
                           Tree treeModel,
                           PatternList patternList,
                           SubstitutionModel substitutionModel,
                           Parameter popSizeParameter,
                           Parameter growthRateParameter,
                           BranchRateModel branchRateModel,
                           BastaLikelihoodDelegate likelihoodDelegate,
                           int numberSubIntervals,
                           boolean useAmbiguities,
                           DataType dataType,
                           String tag,
                           boolean useMAP,
                           boolean returnMarginalLogLikelihood,
                           boolean conditionalProbabilitiesInLogSpace) {

        super(name);

        assert likelihoodDelegate != null;
        assert treeModel != null;
        assert branchRateModel != null;
        assert patternList.getPatternCount() == 1;
        assert useAmbiguities;

        if (!(branchRateModel instanceof StrictClockBranchRates)) {
            throw new RuntimeException("Not yet implemented");
        }

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing BastaLikelihood with Ancestral State Reconstruction");

        this.patternList = patternList;

        this.likelihoodDelegate = likelihoodDelegate;
        addModel(likelihoodDelegate);

        this.tree = treeModel;

        this.branchRateModel = branchRateModel;
        addModel(branchRateModel);

        this.substitutionModel = substitutionModel;
        addModel(substitutionModel);

        this.popSizeParameter = popSizeParameter;
        addVariable(popSizeParameter);

        this.growthRateParameter = growthRateParameter;

        if (this.growthRateParameter != null) {
            addVariable(growthRateParameter);
            this.likelihoodDelegate.updateIsExponentialGrowth(true);
        }

        this.stateCount = substitutionModel.getDataType().getStateCount();

        if (tree instanceof TreeModel) {
            treeIntervals = new BestSignalsFromBigFastTreeIntervals((TreeModel) treeModel);
            addModel(treeIntervals);
        } else {
            throw new RuntimeException("Not yet implemented");
        }

        treeTraversalDelegate = new CoalescentIntervalTraversal(treeModel, treeIntervals, branchRateModel, numberSubIntervals);

        // Initialize ancestral state reconstruction settings
        this.dataType = dataType;
        this.tag = tag;
        this.useMAP = useMAP;
        this.returnMarginalLogLikelihood = returnMarginalLogLikelihood;
        this.conditionalProbabilitiesInLogSpace = conditionalProbabilitiesInLogSpace;

        // Initialize state storage arrays
        reconstructedStates = new int[treeModel.getNodeCount()][patternList.getPatternCount()];
        storedReconstructedStates = new int[treeModel.getNodeCount()][patternList.getPatternCount()];

        // Initialize sub-interval states
        initializeSubIntervalStates(treeIntervals, numberSubIntervals);

        boolean stripHiddenState = false;
        this.formatter = new CodeFormatter(dataType, stripHiddenState);

        // Add tree trait for accessing reconstructed states
        setupTraits();

        setTipData();

        likelihoodKnown = false;
        populationSizesKnown = false;
        treeIntervalsKnown = false;
        transitionMatricesKnown = false;
    }

    /**
     * Initialize storage for sub-interval states
     */

    private void initializeSubIntervalStates(BestSignalsFromBigFastTreeIntervals treeIntervals, int numberSubIntervals) {
        int totalSubIntervals = treeIntervals.getIntervalCount() * numberSubIntervals;
        int nodeCount = tree.getNodeCount();
        int patternCount = getPatternCount();

        // Three-dimensional array: [interval][nodeNumber][pattern]
        subIntervalStates = new int[nodeCount][patternCount];

        for (int i = 0; i < nodeCount; i++) {
            Arrays.fill(subIntervalStates[i], -1);
        }
    }

    /**
     * Set up tree traits for accessing reconstructed states
     */
    private void setupTraits() {
        TreeTrait<int[]> ancestralStateTrait = new TreeTrait.IA() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() {
                return int[].class;
            }

            public int[] getTrait(Tree tree, NodeRef node) {
                return getStatesForNode(tree, node);
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formattedState(getStatesForNode(tree, node), formatter);
            }
        };

        treeTraits.addTrait(ancestralStateTrait);
    }

    public CoalescentIntervalTraversal getTraversalDelegate() { return treeTraversalDelegate; }

    public SubstitutionModel getSubstitutionModel() { return substitutionModel; } // TODO generify for multiple models (e.g. epochs)

    public void setTipData() {

        int[] data = patternList.getPattern(0);

        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef node = tree.getExternalNode(i);

            int index = patternList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int datum = data[index];
            double[] partials = new double[stateCount];
            if (datum >= stateCount) {
                int[] states = new int[1];
                getTipStates(index, states);
                partials[states[0]] = 1.0;
            } else {
                partials[datum] = 1.0;
            }

            likelihoodDelegate.setPartials(node.getNumber(), partials);
        }
    }

    @Override
    public final Model getModel() {
        return this;
    }

    @Override @SuppressWarnings("Duplicates")
    public double getLogLikelihood() {
        if (COUNT_TOTAL_OPERATIONS) totalGetLogLikelihoodCount++;

        if (!likelihoodKnown) {
            long startTime;
            if (COUNT_TOTAL_OPERATIONS) {
                totalCalculateLikelihoodCount++;
                startTime = System.nanoTime();
            }

            logLikelihood = calculateLogLikelihood();

            if (COUNT_TOTAL_OPERATIONS) {
                long endTime = System.nanoTime();
                totalLikelihoodTime += (endTime - startTime) / 1000;
            }

            likelihoodKnown = true;
        }

        if (returnMarginalLogLikelihood) {
            return logLikelihood;
        }

        // Redraw states and return joint density for ancestral state reconstruction
        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }

        return jointLogLikelihood;
    }


    @Override
    public void setTipStates(int tipNum, int[] states) {
        double[] partials = new double[stateCount];
        partials[states[0]] = 1.0;
        likelihoodDelegate.setPartials(tipNum, partials);

        if (reconstructedStates != null && states != null && states.length > 0) {
            if (tipNum < reconstructedStates.length) {
                System.arraycopy(states, 0, reconstructedStates[tipNum], 0, Math.min(states.length, reconstructedStates[tipNum].length));
            }
        }

        likelihoodKnown = false;
        areStatesRedrawn = false;
    }

    @Override
    public void getTipStates(int tipNum, int[] states) {
        double[] partials = new double[stateCount];
        likelihoodDelegate.getPartials(tipNum, partials);
        for (int i = 0; i < stateCount; i++) {
            if (partials[i] > 0) {
                states[0] = i;
                break;
            }
        }
    }

    @Override
    public int getPatternCount() {
        return 1;
    }

    @Override
    public int getTipCount() {
        return tree.getExternalNodeCount();
    }

    @Override
    public void makeDirty() {
        if (COUNT_TOTAL_OPERATIONS) totalMakeDirtyCount++;

        likelihoodKnown = false;
        treeIntervalsKnown = false;
        populationSizesKnown = false;
        transitionMatricesKnown = false;
        areStatesRedrawn = false;

        likelihoodDelegate.makeDirty();
        updateAllNodes();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == popSizeParameter) {
            populationSizesKnown = false;
            likelihoodKnown = false;
        } else if (variable == growthRateParameter) {
            populationSizesKnown = false;
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeIntervals) {
            treeIntervalsKnown = false;
            transitionMatricesKnown = false;
            areStatesRedrawn = false;
            nodeIntervalMap.clear();
        } else if (model == branchRateModel) {
            treeIntervalsKnown = false; // TODO should not be necessary
            transitionMatricesKnown = false;
        } else if (model == substitutionModel) {
            transitionMatricesKnown = false;
        } else {
            throw new RuntimeException("Not yet implemented");
        }

        if (COUNT_TOTAL_OPERATIONS) totalModelChangedCount++;

        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        assert (likelihoodKnown) : "the likelihood should always be known at this point in the cycle";
        assert (populationSizesKnown);
        assert (treeIntervalsKnown);
        assert (transitionMatricesKnown);
        storedLogLikelihood = logLikelihood;

        // Store ancestral state reconstruction information
        if (areStatesRedrawn) {
            for (int i = 0; i < reconstructedStates.length; i++) {
                System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
            }
        }

        storedAreStatesRedrawn = areStatesRedrawn;
        storedJointLogLikelihood = jointLogLikelihood;
    }

    @Override

    protected final void restoreState() {
        //likelihoodDelegate.restoreState();
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
        populationSizesKnown = false;
        treeIntervalsKnown = false;
        transitionMatricesKnown = false;

        // Restore ancestral state reconstruction information
        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
        jointLogLikelihood = storedJointLogLikelihood;
    }

    @Override
    protected void acceptState() { } // nothing to do

    private double calculateLogLikelihood() {

        if (!transitionMatricesKnown) {
            // update eigen-decomposition
            likelihoodDelegate.updateEigenDecomposition(0, substitutionModel.getEigenDecomposition(), false);// TODO do conditionally and double-buffer
        }

        if (!populationSizesKnown) {
            // update population sizes
            likelihoodDelegate.updatePopulationSizes(0, popSizeParameter.getParameterValues(), false); // TODO do conditionally and double-buffer
        }

        if (this.growthRateParameter != null) {
            likelihoodDelegate.updateGrowthRates(0, growthRateParameter.getParameterValues(), false);
        }

        if (!treeIntervalsKnown) {
            // update operations on tree
            treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();
        }

        final List<BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<TransitionMatrixOperation> matrixOperations =
                transitionMatricesKnown ? NO_OPT :
                        treeTraversalDelegate.getMatrixOperations();
        final List<Integer> intervalStarts = treeTraversalDelegate.getIntervalStarts();

        if (!transitionMatricesKnown){
            likelihoodDelegate.flipTransitionMatrixBuffer(matrixOperations);
        }
        if (COUNT_TOTAL_OPERATIONS) {
            totalPropagationCount += branchOperations.size();
            totalMatrixUpdateCount += matrixOperations.size();
            totalIntervalReductionCount += treeTraversalDelegate.getCoalescentIntervalCount();
        }

        final NodeRef root = tree.getRoot();
        double logL = likelihoodDelegate.calculateLikelihood(branchOperations, matrixOperations,
                intervalStarts, root.getNumber(), this, transitionMatricesKnown);

        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        setAllNodesUpdated();

        treeIntervalsKnown = true;
        populationSizesKnown = true;
        transitionMatricesKnown = true;

        return logL;
    }

    public double[] getGradientLogDensity(StructuredCoalescentLikelihoodGradient wrt) {

        final List<BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<TransitionMatrixOperation> matrixOperations =
                transitionMatricesKnown ? NO_OPT :
                        treeTraversalDelegate.getMatrixOperations();
        final List<Integer> intervalStarts = treeTraversalDelegate.getIntervalStarts();

        final NodeRef root = tree.getRoot();

        calculateLogLikelihood(); // TODO Only execute if necessary

        double[] gradient = likelihoodDelegate.calculateGradient(branchOperations, matrixOperations, intervalStarts,
                root.getNumber(), wrt, this);

        return wrt.chainRule(gradient);
    }

    private void setAllNodesUpdated() {
        treeTraversalDelegate.setAllNodesUpdated();
    }

    /**
     * Set update flag for a node only
     */
    protected void updateNode(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS) totalRateUpdateSingleCount++;

        treeTraversalDelegate.updateNode(node);
        likelihoodKnown = false;
    }

    protected void updateAllNodes() {
        if (COUNT_TOTAL_OPERATIONS) totalRateUpdateAllCount++;

        treeTraversalDelegate.updateAllNodes();
        likelihoodKnown = false;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        double logL = getLogLikelihood();

        String delegateString = likelihoodDelegate.getReport();
        if (delegateString != null) {
            sb.append(delegateString);
        }

        sb.append(getClass().getName()).append("(").append(logL).append(")");

        if (COUNT_TOTAL_OPERATIONS)
            sb.append(
                    "\n  propagation operations = ").append(totalPropagationCount).append(
                    "\n  matrix updates = ").append(totalMatrixUpdateCount).append(
                    "\n  interval operations = ").append(totalIntervalReductionCount).append(
                    "\n  model changes = ").append(totalModelChangedCount).append(
                    "\n  make dirties = ").append(totalMakeDirtyCount).append(
                    "\n  calculate likelihoods = ").append(totalCalculateLikelihoodCount).append(
                    "\n  get likelihoods = ").append(totalGetLogLikelihoodCount).append(
                    "\n  all rate updates = ").append(totalRateUpdateAllCount).append(
                    "\n  partial rate updates = ").append(totalRateUpdateSingleCount).append(
                    "\n  average likelihood time = ").append(totalLikelihoodTime / totalCalculateLikelihoodCount);


        return sb.toString();
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public MutableTreeModel getTreeModel() {
        return null;
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    @Override
    public String formattedState(int[] state) {
        return null;
    }

    public void addTrait(TreeTrait trait) {
        treeTraits.addTrait(trait);
    }

    public void addTraits(TreeTrait[] traits) {
        treeTraits.addTraits(traits);
    }

    @Override
    public Citation.Category getCategory() { return Citation.Category.TREE_PRIORS; }

    @Override
    public String getDescription() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getDescription();
        } else {
            return null;
        }
    }

    public BastaLikelihoodDelegate getLikelihoodDelegate() {  return likelihoodDelegate; }

    @Override
    public List<Citation> getCitations() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getCitations();
        } else {
            return new ArrayList<>();
        }
    }

    @Override
    public long getTotalCalculationCount() {
        return likelihoodDelegate.getTotalCalculationCount();
    }

    public int[] getStatesForNode(Tree tree, NodeRef node) {
        if (tree != this.tree) {
            throw new RuntimeException("Can only reconstruct states on tree given to constructor");
        }

        if (!likelihoodKnown) {
            getLogLikelihood();
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNumber()];
    }

    /**
     * Maps nodes to their corresponding coalescent intervals with improved tracking
     */
    private void mapNodeToSubIntervals() {
        nodeIntervalMap.clear();
        CoalescentIntervalTraversal traversal = getTraversalDelegate();
        traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

        List<BranchIntervalOperation> branchIntervalOps = traversal.getBranchIntervalOperations();

        for (BranchIntervalOperation op : branchIntervalOps) {
            int nodeNumber = op.inputBuffer1 % tree.getNodeCount();
            if (!nodeIntervalMap.containsKey(nodeNumber)) {
                nodeIntervalMap.put(nodeNumber, new ArrayList<>());
            }
            nodeIntervalMap.get(nodeNumber).add(op.intervalNumber);

            if (op.inputBuffer2 > 0) {
                int nodeNumber2 = op.inputBuffer2 % tree.getNodeCount();
                if (!nodeIntervalMap.containsKey(nodeNumber2)) {
                    nodeIntervalMap.put(nodeNumber2, new ArrayList<>());
                }
                nodeIntervalMap.get(nodeNumber2).add(op.intervalNumber);
            }
        }

        // Sort interval lists for each node to ensure correct processing order
        for (List<Integer> intervals : nodeIntervalMap.values()) {
            Collections.sort(intervals);
        }
    }


    public void redrawAncestralStates(int traversalMethod) {
        jointLogLikelihood = 0;

        if (traversalMethod == 0) {
            mapNodeToSubIntervals();
            traverseSampleByCoalescentIntervals();
        } else if (traversalMethod == 1) {
            traverseSampleByNodes();
        } else {
            throw new IllegalArgumentException("Invalid traversal method: " + traversalMethod);
        }

        areStatesRedrawn = true;
    }


    public void redrawAncestralStates() {
        redrawAncestralStates(1); // Use original method by default
    }



    private void traverseSampleByCoalescentIntervals() {

        CoalescentIntervalTraversal traversal = getTraversalDelegate();
        traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();
        List<BranchIntervalOperation> ops = traversal.getBranchIntervalOperations();
        List<TransitionMatrixOperation> matrixOps = traversal.getMatrixOperations();
        List<Integer> intervalStarts = traversal.getIntervalStarts();
        int rootBuffer = 0;
        BranchIntervalOperation.initializeMap(tree, likelihoodDelegate.getMaxNumberOfCoalescentIntervals());
        Map<Integer, Integer> bufferToNodeMap = new HashMap<>();
        BranchIntervalOperation lastOp = ops.get(ops.size() - 1);

        for (BranchIntervalOperation op : ops) {
            int input1NodeNumber = op.inputBuffer1 % tree.getNodeCount();
            int input2NodeNumber = op.inputBuffer2 % tree.getNodeCount();
            int outputNodeNumber = op.outputBuffer % tree.getNodeCount();
            op.transform();
            if (!bufferToNodeMap.containsKey(op.inputBuffer1)) {
                bufferToNodeMap.put(op.inputBuffer1, input1NodeNumber);
            }

            if (op.inputBuffer2 >= 0) {

                if (!bufferToNodeMap.containsKey(op.inputBuffer2)) {
                    bufferToNodeMap.put(op.inputBuffer2, input2NodeNumber);
                }
            }

            if (!bufferToNodeMap.containsKey(op.outputBuffer)) {
                bufferToNodeMap.put(op.outputBuffer, outputNodeNumber);
            }

            if (op == lastOp) {
                rootBuffer = op.outputBuffer;
            }

        }

        Map<Integer, TransitionMatrixOperation> matrixOpMap = new HashMap<>();
        for (TransitionMatrixOperation op : matrixOps) {
            matrixOpMap.put(op.outputBuffer, op);
        }

        NodeRef rootNode = tree.getRoot();
        sampleRootState(rootNode, rootBuffer, true);

        for (int interval = intervalStarts.size() - 2; interval >= 0; interval--) {
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            // Process operations within this interval
            for (int i = start; i < end; i++) {
                BranchIntervalOperation op = ops.get(i);

                if (op.inputBuffer1 >= 0) {
                    processChildBuffer(op, op.inputBuffer1, op.inputMatrix1,
                            bufferToNodeMap, matrixOpMap);
                }

                if (op.inputBuffer2 >= 0) {
                    processChildBuffer(op, op.inputBuffer2, op.inputMatrix2,
                            bufferToNodeMap, matrixOpMap);
                }
            }
        }
    }

    /**
     * Process a child node using its buffer
     */
    private void processChildBuffer(BranchIntervalOperation op, int inputBuffer, int matrixBuffer,
                                    Map<Integer, Integer> bufferToNodeMap,
                                    Map<Integer, TransitionMatrixOperation> matrixOpMap) {

        Integer childNodeNumber = bufferToNodeMap.get(inputBuffer);

        NodeRef childNode = tree.getNode(childNodeNumber);


        double[] transitionMatrix = new double[stateCount * stateCount];
        if (matrixBuffer >= 0 && matrixOpMap.containsKey(matrixBuffer)) {
            TransitionMatrixOperation matrixOp = matrixOpMap.get(matrixBuffer);
            getMatrix(likelihoodDelegate, matrixOp.outputBuffer, transitionMatrix);
        } else {
            throw new RuntimeException("No transition matrix matches");
        }

        if (tree.isExternal(childNode) && inputBuffer < tree.getExternalNodeCount()) {
            double[] partials = new double[stateCount];
            likelihoodDelegate.getPartials(inputBuffer, partials);

            for (int j = 0; j < getPatternCount(); j++) {
                boolean isAmbiguous;
                int unambiguousState = -1;
                int statesWithProbability = 0;

                for (int k = 0; k < stateCount; k++) {
                    if (partials[k] > 0) {
                        statesWithProbability++;
                        unambiguousState = k;
                    }
                }

                isAmbiguous = (statesWithProbability != 1);

                if (isAmbiguous) {
                    sampleStateForSubInterval(childNode, childNodeNumber,
                            inputBuffer, transitionMatrix);
                    reconstructedStates[childNodeNumber][j] = subIntervalStates[childNodeNumber][j];
                } else {
                    subIntervalStates[childNodeNumber][j] = unambiguousState;
                    reconstructedStates[childNodeNumber][j] = unambiguousState;
                }
            }
            return;
        }

        sampleStateForSubInterval(childNode, childNodeNumber,
                inputBuffer, transitionMatrix);


        List<Integer> nodeIntervals = nodeIntervalMap.get(childNodeNumber);
        if (nodeIntervals != null) {
            int lastInterval = nodeIntervals.get(nodeIntervals.size() - 1);
            if (op.intervalNumber == lastInterval) {
                System.arraycopy(
                        subIntervalStates[childNodeNumber], 0,
                        reconstructedStates[childNodeNumber], 0,
                        getPatternCount()
                );
            }
        }
    }


    private void sampleRootState(NodeRef root, int rootBuffer, boolean copyToIntervals) {
        int nodeNum = root.getNumber();

        double[] partials = new double[stateCount];
        likelihoodDelegate.getPartials(rootBuffer, partials);

        double[] conditionalProbabilities = new double[stateCount];
        double[] frequencies = substitutionModel.getFrequencyModel().getFrequencies();

        for (int j = 0; j < getPatternCount(); j++) {
            for (int i = 0; i < stateCount; i++) {
                if (conditionalProbabilitiesInLogSpace) {
                    conditionalProbabilities[i] = Math.log(partials[i]) + Math.log(frequencies[i]);
                } else {
                    conditionalProbabilities[i] = partials[i] * frequencies[i];
                }
            }

            reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);

            if (copyToIntervals) {
                    subIntervalStates[nodeNum][j] = reconstructedStates[nodeNum][j];
            }

            if (!returnMarginalLogLikelihood) {
                jointLogLikelihood += Math.log(frequencies[reconstructedStates[nodeNum][j]]);
            }
        }
    }


    private void sampleStateForSubInterval(NodeRef node, int nodeNumber,
                                           int bufferIndex, double[] transitionMatrix) {
        int patternCount = getPatternCount();

        double[] nodeLikelihoods = new double[stateCount];
        likelihoodDelegate.getPartials(bufferIndex, nodeLikelihoods);

        int[] parentState = null;
        if (subIntervalStates[nodeNumber][0] != -1) {
            parentState = subIntervalStates[nodeNumber];
        }

        if (parentState == null) {
            NodeRef parent = tree.getParent(node);
            if (parent == null) {
                return;
            }

            int parentNumber = parent.getNumber();
            parentState = subIntervalStates[parentNumber];
        }

        double[] conditionalProbabilities = new double[stateCount];

        for (int j = 0; j < patternCount; j++) {
            int parentIndex = parentState[j];
            for (int i = 0; i < stateCount; i++) {
                if (conditionalProbabilitiesInLogSpace) {
                    conditionalProbabilities[i] = Math.log(transitionMatrix[i * stateCount + parentIndex]) +
                            Math.log(nodeLikelihoods[i]);
                } else {
                    conditionalProbabilities[i] = transitionMatrix[i * stateCount + parentIndex] *
                            nodeLikelihoods[i];
                }
            }

            subIntervalStates[nodeNumber][j] = drawChoice(conditionalProbabilities);
            if (!returnMarginalLogLikelihood) {
                double contrib = transitionMatrix[parentIndex + subIntervalStates[nodeNumber][j]];
                jointLogLikelihood += Math.log(contrib);
            }
        }
    }

    private void getMatrix(BastaLikelihoodDelegate delegate, int bufferIndex, double[] probabilities) {
        double[] matrices = new double[stateCount * stateCount];
        delegate.getTransitionMatrices(bufferIndex, matrices);
        System.arraycopy(
                matrices,
                0,
                probabilities, 0,
                stateCount * stateCount
        );
    }


    private void traverseSampleByNodes() {
        CoalescentIntervalTraversal traversal = getTraversalDelegate();
        traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

        List<BranchIntervalOperation> originalOps = traversal.getBranchIntervalOperations();

        int maxNumCoalescentIntervals = 0;
        BastaLikelihoodDelegate delegate = likelihoodDelegate;
        maxNumCoalescentIntervals = delegate.getMaxNumberOfCoalescentIntervals();


        BranchIntervalOperation.initializeMap(tree, maxNumCoalescentIntervals);
        Map<Integer, Integer> nodeToBufferMap = new HashMap<>();
        for (int i = 0; i < originalOps.size(); i++) {
            BranchIntervalOperation originalOp = originalOps.get(i);
            int input1NodeNumber = originalOp.inputBuffer1 % tree.getNodeCount();
            int input2NodeNumber = originalOp.inputBuffer2 % tree.getNodeCount();
            int outputNodeNumber = originalOp.outputBuffer % tree.getNodeCount();
            originalOp.transform();
            if (!nodeToBufferMap.containsKey(input1NodeNumber)) {
                    nodeToBufferMap.put(input1NodeNumber, originalOp.inputBuffer1);
            }

            if (originalOp.inputBuffer2 >= 0) {

                if (!nodeToBufferMap.containsKey(input2NodeNumber)) {
                    nodeToBufferMap.put(input2NodeNumber, originalOp.inputBuffer2);
                }
            }

            if (!nodeToBufferMap.containsKey(outputNodeNumber)) {
                nodeToBufferMap.put(outputNodeNumber, originalOp.outputBuffer);
            }

        }


        NodeRef rootNode = tree.getRoot();
        int rootBuffer = nodeToBufferMap.getOrDefault(rootNode.getNumber(), rootNode.getNumber());
        sampleRootState(rootNode, rootBuffer, false);
        traverseNodesPreOrderDirect(rootNode, nodeToBufferMap);
    }

    /**
     * Traverse the tree in pre-order with direct matrix computation
     */
    private void traverseNodesPreOrderDirect(NodeRef node, Map<Integer, Integer> nodeToBufferMap) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);

            if (tree.isExternal(child)) {
                int childNum = child.getNumber();
                int parentNum = node.getNumber();

                double[] partials = new double[stateCount];
                int childBuffer = nodeToBufferMap.getOrDefault(childNum, childNum);
                likelihoodDelegate.getPartials(childBuffer, partials);

                for (int j = 0; j < getPatternCount(); j++) {
                    // Check for unambiguous state
                    boolean isAmbiguous;
                    int unambiguousState = -1;
                    int statesWithProbability = 0;

                    for (int k = 0; k < stateCount; k++) {
                        if (partials[k] > 0) {
                            statesWithProbability++;
                            unambiguousState = k;
                        }
                    }

                    isAmbiguous = (statesWithProbability != 1);

                    if (isAmbiguous) {
                        int parentState = reconstructedStates[parentNum][j];
                        int parentIndex = parentState;
                        double branchLength = tree.getNodeHeight(node) - tree.getNodeHeight(child);
                        double[] transitionMatrix = new double[stateCount * stateCount];
                        substitutionModel.getTransitionProbabilities(branchLength, transitionMatrix);

                        double[] conditionalProbabilities = new double[stateCount];
                        for (int k = 0; k < stateCount; k++) {
                            if (conditionalProbabilitiesInLogSpace) {
                                conditionalProbabilities[k] = Math.log(transitionMatrix[k * stateCount + parentIndex]) +
                                        Math.log(partials[k] > 0 ? partials[k] : Double.MIN_VALUE);
                            } else {
                                conditionalProbabilities[k] = transitionMatrix[k * stateCount + parentIndex] * partials[k];
                            }
                        }
                        reconstructedStates[childNum][j] = drawChoice(conditionalProbabilities);

                    } else {
                        reconstructedStates[childNum][j] = unambiguousState;
                    }
                }
                continue;
            }

            double branchLength = tree.getNodeHeight(node) - tree.getNodeHeight(child);
            double[] transitionMatrix = new double[stateCount * stateCount];
            substitutionModel.getTransitionProbabilities(branchLength, transitionMatrix);

            int childBuffer = nodeToBufferMap.getOrDefault(child.getNumber(), child.getNumber());

            sampleChildNodeState(node, child, childBuffer, transitionMatrix);

            traverseNodesPreOrderDirect(child, nodeToBufferMap);
        }
    }



    private void sampleChildNodeState(NodeRef parent, NodeRef child, int bufferIndex, double[] transitionMatrix) {
        int parentNumber = parent.getNumber();
        int childNumber = child.getNumber();
        int patternCount = getPatternCount();

        double[] childPartials = new double[stateCount];
        likelihoodDelegate.getPartials(bufferIndex, childPartials);

        double[] conditionalProbabilities = new double[stateCount];

        for (int j = 0; j < patternCount; j++) {
            int parentState = reconstructedStates[parentNumber][j];
            int parentIndex = parentState * stateCount;

            for (int i = 0; i < stateCount; i++) {
                double transProb = transitionMatrix[parentIndex + i];
                double likelihood = childPartials[i];

                conditionalProbabilities[i] = conditionalProbabilitiesInLogSpace ?
                        Math.log(transProb) + Math.log(likelihood) :
                        transProb * likelihood;
            }

            reconstructedStates[childNumber][j] = drawChoice(conditionalProbabilities);

            if (!returnMarginalLogLikelihood) {
                double contrib = transitionMatrix[parentIndex + reconstructedStates[childNumber][j]];
                jointLogLikelihood += Math.log(contrib);
            }
        }
    }


    private int drawChoice(double[] measure) {
        if (useMAP) {
            // Use Maximum A Posteriori
            double max = measure[0];
            int choice = 0;
            for (int i = 1; i < measure.length; i++) {
                if (measure[i] > max) {
                    max = measure[i];
                    choice = i;
                }
            }
            return choice;
        } else {
            if (conditionalProbabilitiesInLogSpace) {
                return MathUtils.randomChoiceLogPDF(measure);
            }
            return MathUtils.randomChoicePDF(measure);
        }
    }

    private static String formattedState(int[] state, CodeFormatter formatter) {
        StringBuffer sb = new StringBuffer();
        sb.append("\"");
        formatter.reset();
        for (int i : state) {
            sb.append(formatter.getCodeString(i));
        }
        sb.append("\"");
        return sb.toString();
    }


    private class CodeFormatter {
        private final DataType dataType;
        private final Function<String, String> appender;
        private final Function<Integer, String> getter;
        private boolean first = true;

        CodeFormatter(DataType dataType, boolean stripHiddenState) {
            this.dataType = dataType;

            this.appender = (dataType instanceof GeneralDataType) ?
                    (codeString) -> codeString + " " : Function.identity();

            if (dataType instanceof HiddenCodons) {
                this.getter = (stripHiddenState) ?
                        ((HiddenCodons) dataType)::getTripletWithoutHiddenCode :
                        dataType::getTriplet;
            } else if (dataType instanceof HiddenDataType && stripHiddenState) {
                this.getter = ((HiddenDataType) dataType)::getCodeWithoutHiddenState;
            } else {
                this.getter = dataType::getCode;
            }
        }

        String getCodeString(int state) {
            String code = getter.apply(state);
            if (first) {
                first = false;
            } else {
                code = appender.apply(code);
            }
            return code;
        }

        void reset() { first = true; }
    }

    private final List<TransitionMatrixOperation> NO_OPT = new ArrayList<>();

    private int totalPropagationCount = 0;
    private int totalMatrixUpdateCount = 0;
    private int totalIntervalReductionCount = 0;
    private int totalGetLogLikelihoodCount = 0;
    private int totalModelChangedCount = 0;
    private int totalMakeDirtyCount = 0;
    private int totalCalculateLikelihoodCount = 0;
    private int totalRateUpdateAllCount = 0;
    private int totalRateUpdateSingleCount = 0;
    private int totalSizeUpdateAllCount = 0;

    private long totalLikelihoodTime = 0;

    public Parameter getPopSizes() {return popSizeParameter;}
}