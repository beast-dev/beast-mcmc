/*
 * AbstractBeagleGradientDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.HiddenCodons;
import dr.evolution.datatype.HiddenDataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.basta.BastaLikelihood;
import dr.evomodel.coalescent.basta.BastaLikelihoodDelegate;
import dr.evomodel.coalescent.basta.CoalescentIntervalTraversal;
import dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.markovjumps.CompleteHistoryAddOn;
import dr.evomodel.treelikelihood.MarkovJumpsTraitProvider;
import dr.inference.model.Model;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.*;
import java.util.function.Function;

/**
 * @author Yucai Shao
 * @author Marc Suchard
 */
public abstract class AbstractRealizedDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate {

    public static String NAME_SUFFIX = "_unified";

    public static class Fit extends AbstractRealizedDiscreteTraitDelegate {

        private final BeagleDataLikelihoodDelegate beagleDataLikelihoodDelegate;

        public Fit(String name,
                   TreeDataLikelihood likelihood,
                   BeagleDataLikelihoodDelegate likelihoodDelegate,
                   boolean useMap) {
            super(name, likelihood, likelihoodDelegate, likelihoodDelegate.getUseAmbiguities(), useMap);
            this.beagleDataLikelihoodDelegate = likelihoodDelegate;
        }

        @Override
        double[] getRootStateFrequencies() {
            return evolutionaryProcessDelegate.getRootStateFrequencies();
        }

        @Override
        int[] getTipStates(int node) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        double[] getPostOrderPartials(int node) {
            beagleDataLikelihoodDelegate.getPartials(node, tmpPostOrder);
            return tmpPostOrder;
        }

        @Override
        public double[] getTransitionMatrix(int node, int parent) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public OrderedEvents getOrderedInformation(int node, int parent, double[] probabilities) {


            double nodeTime = tree.getNodeHeight(tree.getNode(node));
            double parentTime = tree.getNodeHeight(tree.getNode(parent));

            double branchTime = parentTime - nodeTime;
            double branchRate = branchRateModel.getBranchRate(tree, tree.getNode(node));

            return new OrderedEvents(reconstructedStates[parent], reconstructedStates[node],
                    siteSpecificRateCategory, probabilities,
                    parentTime, nodeTime, branchTime, branchRate);
        }
    }

    public static class Bit extends AbstractRealizedDiscreteTraitDelegate {

        private final BastaLikelihood likelihood;
        private final BastaLikelihoodDelegate likelihoodDelegate;
        private final SubstitutionModel substitutionModel;
        private final boolean transpose;

        private EigenDecomposition eigenDecomposition;
        private Map<Integer, Integer> nodeMap;

        public Bit(String name,
                   BastaLikelihood likelihood,
                   boolean useMAP) {
            super(name, likelihood, likelihood, true, useMAP);

            this.likelihood = likelihood;
            this.likelihoodDelegate = likelihood.getLikelihoodDelegate();
//            this.branchRateModel = likelihood.getBranchRateModel();
            this.substitutionModel = likelihood.getSubstitutionModel();
            this.rootStateFrequencies = new double[stateCount];
            this.transpose = true; // TODO pass info
            Arrays.fill(rootStateFrequencies, 1.0);
        }

        private static final boolean TRANSPOSE_ONCE = true;

        @Override
        public OrderedEvents getOrderedInformation(int node, int parent, double[] unused) {

            double[] probabilities = getTransitionMatrix(tree.getNode(node), tree.getNode(parent), false);

            double nodeTime = tree.getNodeHeight(tree.getNode(node));
            double parentTime = tree.getNodeHeight(tree.getNode(parent));
            double branchTime = parentTime - nodeTime;
            double branchRate = branchRateModel.getBranchRate(tree, tree.getNode(node));

            return new OrderedEvents(reconstructedStates[node], reconstructedStates[parent],
                    siteSpecificRateCategory, probabilities, nodeTime, parentTime, branchTime, branchRate);
        }

        @Override
        public void setupStatistics() {

            super.setupStatistics();

            EigenDecomposition original = substitutionModel.getEigenDecomposition();
            if (TRANSPOSE_ONCE && transpose) {
                eigenDecomposition = original.transpose();
            } else {
                eigenDecomposition = original;
            }

            nodeMap = getNodeToBufferMap();

            if (MAS_DEBUG) {
                MathUtils.setSeed(666);
            }
        }

        private Map<Integer, Integer> getNodeToBufferMap() {

            CoalescentIntervalTraversal traversal = likelihood.getTraversalDelegate();
            traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

            List<ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation> originalOps = traversal.getBranchIntervalOperations();

            int maxNumCoalescentIntervals = likelihoodDelegate.getMaxNumberOfCoalescentIntervals();

            ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation.initializeMap(tree, maxNumCoalescentIntervals);
            Map<Integer, Integer> nodeToBufferMap = new HashMap<>();
            for (ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation originalOp : originalOps) {
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

            return nodeToBufferMap;
        }

        private int map(final int node) {
            return nodeMap.getOrDefault(node, node);
            // TODO map should be `int[]` and can be constructed during Traversal vectorization
        }

        @Override
        double[] getRootStateFrequencies() {
            return rootStateFrequencies;
        }

        @Override
        int[] getTipStates(int node) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        double[] getPostOrderPartials(int node) {
            likelihoodDelegate.getPartials(map(node), tmpPostOrder);
            return tmpPostOrder;
        }

        @Override
        public double[] getTransitionMatrix(int node, int parent) {
            return getTransitionMatrix(tree.getNode(node), tree.getNode(parent), true);
        }

        //        @Override
        double[] getTransitionMatrix(NodeRef node, NodeRef parent, boolean transpose) {

            double branchTime = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
            double branchLength = branchRateModel.getBranchRate(tree, node) * branchTime;
            double[] transitionMatrix = new double[stateCount * stateCount];
            if (transpose) {
                substitutionModel.getTransitionProbabilities(branchLength, transitionMatrix, eigenDecomposition);
            } else {
                substitutionModel.getTransitionProbabilities(branchLength, transitionMatrix);
            }
            if (!TRANSPOSE_ONCE) {
                EigenDecomposition.transposeInPlace(transitionMatrix, stateCount); // TODO save work by transposing eigen-decomposition once
            }
            return transitionMatrix;
        }

        private final double[] rootStateFrequencies;
    }

    final BranchRateModel branchRateModel;

    protected AbstractRealizedDiscreteTraitDelegate(String name,
                                                    ProcessAlongTree likelihood,
                                                    DiscreteProcessAlongTree likelihoodDelegate,
                                                    boolean useTipPartials,
                                                    boolean useMAP) {
        super(name, likelihood.getTree());
        this.tree = likelihood.getTree();
        this.likelihood = likelihood;
        this.likelihoodDelegate = likelihoodDelegate;
        this.useTipPartials = useTipPartials;
        this.branchRateModel = likelihood.getBranchRateModel();

        this.evolutionaryProcessDelegate = likelihoodDelegate.getEvolutionaryProcessDelegate();
        this.siteRateModel = likelihoodDelegate.getSiteRateModel();

        this.patternList = likelihoodDelegate.getPatternList();
        this.dataType = patternList.getDataType();
        this.patternCount = patternList.getPatternCount();
        this.stateCount = dataType.getStateCount();
        this.categoryCount = siteRateModel.getCategoryCount();

        likelihood.addModelListener(this);
        likelihood.addModelRestoreListener(this);

        boolean stripHiddenState = false;
        this.formatter = new CodeFormatter(dataType, stripHiddenState);

        this.ancestralStatesKnown = false;

        this.tmpPostOrder = new double[stateCount * patternCount * categoryCount];
        this.tmpConditional = new double[stateCount];
        this.tmpTransitionMatrix = new double[stateCount * stateCount * categoryCount];

        this.useMAP = useMAP;

        boolean logCompleteHistory = true;
        CompleteHistoryAddOn mj = new CompleteHistoryAddOn("name", tree,
                likelihoodDelegate.getEvolutionaryProcessDelegate().getBranchSubstitutionModel().getSubstitutionModels(),
                likelihoodDelegate.getSiteRateModel(), this,
                MarkovJumpsTraitProvider.ValueScaling.RAW, logCompleteHistory);

        registerAddOn(mj);
    }

    public void registerAddOn(RealizedDiscreteAddOn addOn) {
        addOn.constructTraits(treeTraitHelper);
        addOns.add(addOn);
    }

    @Override
    public boolean isVectorized() {
        return true;
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {

        TreeTrait<int[]> ancestralStateTrait = new TreeTrait.IA() {

            public String getTraitName() {
                return name;
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

        treeTraitHelper.addTrait(ancestralStateTrait);
    }

    public int[] getStatesForNode(Tree tree, NodeRef node) {

        if (tree != this.tree) {
            throw new RuntimeException("Can only reconstruct states on tree given to constructor");
        }

        simulationProcess.cacheSimulatedTraits(null);
        return reconstructedStates[node.getNumber()];
    }

    public int drawChoice(double[] measure) {
        if (useMAP) {
            // Use maximum a posteriori
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
            return MathUtils.randomChoicePDF(measure);
        }
    }

    @Override
    public void setupStatistics() {

        if (reconstructedStates == null) {
            reconstructedStates = new int[tree.getNodeCount()][patternCount];
        }

        if (siteSpecificRateCategory == null) {
            siteSpecificRateCategory = new int[patternCount];
        }

        // Simulate site-specific rate categories
        if (categoryCount == 1) {
            Arrays.fill(siteSpecificRateCategory, 0);
        } else {

            double[] rootPartials = getPostOrderPartials(tree.getRoot().getNumber());
            double[] posteriorWeightedCategory = tmpConditional;
            double[] priorWeightedCategory = siteRateModel.getCategoryProportions();

            for (int j = 0; j < patternCount; ++j) {
                for (int r = 0; r < categoryCount; r++) {
                    posteriorWeightedCategory[r] = 0;
                    for (int k = 0; k < stateCount; k++) {
                        posteriorWeightedCategory[r] += rootPartials[r * stateCount * patternCount +
                                j * stateCount + k];
                    }
                    posteriorWeightedCategory[r] *= priorWeightedCategory[r];
                }
                siteSpecificRateCategory[j] = drawChoice(posteriorWeightedCategory);
            }
        }
    }

    protected int[] siteSpecificRateCategory;

    abstract double[] getRootStateFrequencies();

    abstract int[] getTipStates(final int node);

    abstract double[] getPostOrderPartials(final int node);

    public abstract double[] getTransitionMatrix(final int node, final int parent);

//    public abstract int[] getStartingStates(int[] nodeStates,  int[] parentStates);

//    public abstract int[] getEndingStates(int[] nodeStates, )

    public abstract OrderedEvents getOrderedInformation(int node, int parent, double[] probabilities);

    public static class OrderedEvents {

        int[] startingStates;
        int[] endingStates;
        int[] categories;

        double[] probabilities;

        double startingTime;
        double endingTime;

        double branchTime;
        double branchRate;

        public OrderedEvents(int[] startingStates, int[] endingStates,
                             int[] categories, double[] probabilities,
                             double startingTime, double endingTime,
                             double branchTime, double branchRate) {
            this.startingStates = startingStates;
            this.endingStates = endingStates;
            this.categories = categories;
            this.probabilities = probabilities;

            if (startingTime == -0.0) {
                startingTime = 0.0;
            }

            this.startingTime = startingTime;
            this.endingTime = endingTime;
            this.branchTime = branchTime;
            this.branchRate = branchRate;

        }

//        public void normalize(StateHistory stateHistory) { }

        public double[] getTransitionMatrix() {
            return probabilities;
        }

        public double getBranchTime() {
            return branchTime;
        }

        public double getBranchRate() {
            return branchRate;
        }

        public int[] getCategories() {
            return categories;
        }

        public int[] getStartingStates() {
            return startingStates;
        }

        public int[] getEndingStates() {
            return endingStates;
        }

        public double getStartingTime() {
            return startingTime;
        }

        public double getEndingTime() {
            return endingTime;
        }
    }

    @Override
    protected void simulateRoot(final int rootNumber) {

        final double[] frequencies = getRootStateFrequencies();
        final double[] postOrderPartial = getPostOrderPartials(rootNumber);

        double[] conditional = tmpConditional;
        for (int j = 0; j < patternCount; j++) {

            int category = siteSpecificRateCategory[j];
            int offset = category * stateCount * patternCount + j * patternCount;

            for (int i = 0; i < stateCount; i++) {
                conditional[i] = postOrderPartial[offset + i] * frequencies[i];
            }

            reconstructedStates[rootNumber][j] = drawChoice(conditional);
        }
    }

    @Override
    protected void simulateNode(int parent, int node, int v2, int v3, int external) {

        if (external == 1) {
            if (useTipPartials) {
                simulateExternalNodeWithPartials(node, parent);
            } else {
                simulateExternalNodeWithStates(node, parent);
            }
        } else {
            simulateInternalNode(node, parent);
        }
    }

    @SuppressWarnings("unused")
    private void simulateExternalNodeWithStates(final int node, final int parent) {

        int[] tipStates = getTipStates(node);
//
//        getTipStates(nodeNum, reconstructedStates[nodeNum]);
//
//        // Check for ambiguity codes and sample them
//
//        for (int j = 0; j < patternCount; j++) {
//
//            final int thisState = reconstructedStates[nodeNum][j];
//
//            if (dataType.isAmbiguousState(thisState)) {
//
//                final int parentIndex = parentState[j] * stateCount;
//                int category = rateCategory == null ? 0 : rateCategory[j];
//                int matrixIndex = category * stateCount * stateCount;
//
//                getMatrix(nodeNum, probabilities);
//                System.arraycopy(probabilities, parentIndex + matrixIndex, conditionalProbabilities, 0, stateCount);
//
//                if (useAmbiguities && !dataType.isUnknownState(thisState)) { // Not completely unknown
//                    boolean[] stateSet = dataType.getStateSet(thisState);
//
//                    for (int k = 0; k < stateCount; k++) {
//                        if (!stateSet[k]) {
//                            conditionalProbabilities[k] = 0.0;
//                        }
//                    }
//                }
//
//                if (conditionalProbabilitiesInLogSpace) {
//                    for (int k = 0; k < stateCount; k++) {
//                        conditionalProbabilities[k] = Math.log(conditionalProbabilities[k]);
//                    }
//                }
//                reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
//            }


        throw new RuntimeException("Not yet implemented");
    }

    private void simulateExternalNodeWithPartials(final int node, final int parent) {

        double[] partialLikelihood = getPostOrderPartials(node);
        double[] conditionalProbabilities = tmpConditional;
        double[] probabilities = null;
        int[] parentStates = reconstructedStates[parent];

        for (int j = 0; j < patternCount; j++) {

            int parentState = parentStates[j];
            int category = siteSpecificRateCategory[j];

            int matrixOffset = category * stateCount * stateCount;
            int partialOffset = category * stateCount * patternCount + j * stateCount;

            // Check for unambiguous state
            int unambiguousState = -1;
            int statesWithProbability = 0;

            for (int i = 0; i < stateCount; ++i) {
                if (partialLikelihood[partialOffset + i] > 0.0) {
                    ++statesWithProbability;
                    unambiguousState = i;
                }
            }

            boolean isAmbiguous = (statesWithProbability != 1);

            final int nodeState;
            if (isAmbiguous) {

                if (probabilities == null) {
                    probabilities = getTransitionMatrix(node, parent);
                }

                for (int i = 0; i < stateCount; ++i) {
                    conditionalProbabilities[i] = partialLikelihood[partialOffset + i]
                            * probabilities[matrixOffset + parentState * stateCount + i];
                }

                nodeState = drawChoice(conditionalProbabilities);
            } else {
                nodeState = unambiguousState;
            }

            reconstructedStates[node][j] = nodeState;

            if (MAS_DEBUG) {
                String id = parent + "->" + node + " ";
                if (isAmbiguous) {
                    System.err.println("new: " + id + new WrappedVector.Raw(conditionalProbabilities));
                } else {
                    System.err.println("new: " + id + unambiguousState);
                }
            }
        }

        hook(node, parent, probabilities);
    }

    protected void hook(final int node, final int parent, double[] probabilities) {

        for (RealizedDiscreteAddOn addOn : addOns) {
            addOn.hookCalculation(node, parent, probabilities);
        }
    }

//    private final MarkovJumpsOnTreeWrapper mj;

    private void simulateInternalNode(final int node, final int parent) {

        double[] partialLikelihood = getPostOrderPartials(node);
        double[] conditionalProbabilities = tmpConditional;
        double[] probabilities = getTransitionMatrix(node, parent);
        int[] parentStates = reconstructedStates[parent];

        for (int j = 0; j < patternCount; ++j) {

            int parentState = parentStates[j];
            int category = siteSpecificRateCategory[j];

            int matrixOffset = category * stateCount * stateCount;
            int partialOffset = category * stateCount * patternCount + j * stateCount;

            for (int i = 0; i < stateCount; ++i) {
                conditionalProbabilities[i] = partialLikelihood[partialOffset + i]
                        * probabilities[matrixOffset + parentState * stateCount + i];
            }

            if (MAS_DEBUG) {
                String id = parent + "->" + node + " ";
                System.err.println("new: " + id + new WrappedVector.Raw(conditionalProbabilities));
                System.err.println("new: " + id + new WrappedVector.Raw(probabilities));
                if (MAS_KILL) {
                    if (count > 10) {
                        System.exit(-1);
                    }
                }
                ++count;
            }

            reconstructedStates[node][j] =  drawChoice(conditionalProbabilities);
        }

        hook(node, parent, probabilities);
    }

    public interface RealizedDiscreteAddOn {

        void constructTraits(Helper treeTraitHelper);

        void hookCalculation(final int node, final int parent, double[] probabilities);
    }

    private final List<RealizedDiscreteAddOn> addOns = new ArrayList<>();

    private static final boolean MAS_KILL = false;
    private static final boolean MAS_DEBUG = false;
    private int count = 0;

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        ancestralStatesKnown = false;
    }

    @Override
    public void modelRestored(Model model) {
        ancestralStatesKnown = false;
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {

        int k = 0;
        for (ProcessOnTreeDelegate.NodeOperation op : nodeOperations) {

            operations[k    ] = op.getNodeNumber(); // Parent sample
            operations[k + 1] = op.getLeftChild();  // Node sample
            operations[k + 4] = op.getLeftChild() < tree.getExternalNodeCount() ? 1 : 0;    // Is node external?

            k += 5;
        }

        return nodeOperations.size();
    }

    private static String formattedState(int[] state, CodeFormatter formatter) {
        StringBuilder sb = new StringBuilder();
        sb.append("\"");
        formatter.reset();
        for (int i : state) {
            sb.append(formatter.getCodeString(i));
        }
        sb.append("\"");
        return sb.toString();
    }

    private static class CodeFormatter {

        private final Function<String, String> appender;
        private final Function<Integer, String> getter;
        private boolean first = true;

        CodeFormatter(DataType dataType, boolean stripHiddenState) {

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

    @Override
    public int getSingleOperationSize() {
        return Beagle.OPERATION_TUPLE_SIZE;
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    protected int[][] reconstructedStates;

    protected final Tree tree;
    protected final ProcessAlongTree likelihood;
    protected final DiscreteProcessAlongTree likelihoodDelegate;
    protected final EvolutionaryProcessDelegate evolutionaryProcessDelegate;
    protected final SiteRateModel siteRateModel;
    protected final PatternList patternList;
    protected final DataType dataType;
    protected final CodeFormatter formatter;

    protected final int patternCount;
    protected final int stateCount;
    protected final int categoryCount;

    protected boolean ancestralStatesKnown;

    protected double[] tmpPostOrder;
    protected double[] tmpConditional;
    protected double[] tmpTransitionMatrix;

    private final boolean useTipPartials;
    private final boolean useMAP;
}
