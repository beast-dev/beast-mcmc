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
import dr.inference.model.Model;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;

/**
 * @author Yucai Shao
 * @author Marc Suchard
 */
public abstract class AbstractRealizedDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate {

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
        double[] getTransitionMatrix(int node, int parent) {
            throw new RuntimeException("Not yet implemented");
        }

//        @Override
//        double[] getPostOrderPartials(final NodeRef node) {
//            beagleDataLikelihoodDelegate.getPartials(node.getNumber(), tmpPostOrder);
//            return tmpPostOrder;
//        }

//        @Override
//        double[] getTransitionMatrix(NodeRef node, NodeRef parent) {
//            throw new RuntimeException("TODO");
//        }

//        @Override
//        int getActiveNodeIndex(int nodeNumber) {
//            return -1;
//        }
//
//        @Override
//        int getActiveMatrixIndex(int nodeNumber) {
//            return -1;
//        }
    }

    public static class Bit extends AbstractRealizedDiscreteTraitDelegate {

        private final BastaLikelihood likelihood;
        private final BastaLikelihoodDelegate likelihoodDelegate;
        private final BranchRateModel branchRateModel;
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
            this.branchRateModel = likelihood.getBranchRateModel();
            this.substitutionModel = likelihood.getSubstitutionModel();
            this.rootStateFrequencies = new double[stateCount];
            this.transpose = true; // TODO pass info
            Arrays.fill(rootStateFrequencies, 1.0);
        }

        private static final boolean TRANSPOSE_ONCE = false; // TODO Make work

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

//        private int map(final NodeRef node) {
//            return nodeMap.getOrDefault(node.getNumber(), node.getNumber());
//        }

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

//        @Override
//        double[] getPostOrderPartials(final NodeRef node) {
//            likelihoodDelegate.getPartials(map(node), tmpPostOrder);
//            return tmpPostOrder;
//        }

        @Override
        double[] getPostOrderPartials(int node) {
            likelihoodDelegate.getPartials(map(node), tmpPostOrder);
            return tmpPostOrder;
        }

        @Override
        double[] getTransitionMatrix(int node, int parent) {
            return getTransitionMatrix(tree.getNode(node), tree.getNode(parent));
        }

//        @Override
        double[] getTransitionMatrix(NodeRef node, NodeRef parent) {

            double branchTime = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
            double branchLength = branchRateModel.getBranchRate(tree, node) * branchTime;
            double[] transitionMatrix = new double[stateCount * stateCount];
            substitutionModel.getTransitionProbabilities(branchLength, transitionMatrix, eigenDecomposition);
            if (!TRANSPOSE_ONCE) {
                EigenDecomposition.transposeInPlace(transitionMatrix, stateCount); // TODO save work by transposing eigen-decomposition once
            }
            return transitionMatrix;
        }

//        @Override
//        int getActiveNodeIndex(int nodeNumber) {
//            return -1;
//        }
//
//        @Override
//        int getActiveMatrixIndex(int nodeNumber) {
//            return -1;
//        }

        private final double[] rootStateFrequencies;
    }

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

    private int[] siteSpecificRateCategory;

    abstract double[] getRootStateFrequencies();

    abstract int[] getTipStates(final int node);

    abstract double[] getPostOrderPartials(final int node);

    abstract double[] getTransitionMatrix(final int node, final int parent);

//    abstract double[] getPostOrderPartials(final NodeRef node);

//    abstract double[] getTransitionMatrix(final NodeRef node, final NodeRef parent);

//    @Override
//    protected void simulateRoot(final NodeRef root) {
//
//        final double[] frequencies = getRootStateFrequencies();
//        final double[] postOrderPartial = getPostOrderPartials(root);
//
//        double[] conditional = tmpConditional;
//        for (int j = 0; j < patternCount; j++) {
//
//            int category = siteSpecificRateCategory[j];
//            int offset = category * stateCount * patternCount + j * patternCount;
//
//            for (int i = 0; i < stateCount; i++) {
//                conditional[i] = postOrderPartial[offset + i] * frequencies[i];
//            }
//
//            reconstructedStates[root.getNumber()][j] = drawChoice(conditional);
//        }
//    }

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

//    @Override
//    protected void simulateNode(NodeOperation op) {
//
//        final NodeRef parent = tree.getNode(op.getNodeNumber());
//        final NodeRef node = tree.getNode(op.getLeftChild());
//
//        if (tree.isExternal(node)) {
//            if (useTipPartials) {
//                simulateExternalNodeWithPartials(node, parent);
//            } else {
//                simulateExternalNodeWithStates(node, parent);
//            }
//        } else {
//            simulateInternalNode(node, parent);
//        }
//    }

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

            if (isAmbiguous) {

                if (probabilities == null) {
                    probabilities = getTransitionMatrix(node, parent);
                }

                for (int i = 0; i < stateCount; ++i) {
                    conditionalProbabilities[i] = partialLikelihood[partialOffset + i]
                            * probabilities[matrixOffset + parentState * stateCount + i];
                }

                reconstructedStates[node][j] = drawChoice(conditionalProbabilities);
            } else {
                reconstructedStates[node][j] = unambiguousState;
            }

            if (MAS_DEBUG) {
                String id = parent + "->" + node + " ";
                if (isAmbiguous) {
                    System.err.println("new: " + id + new WrappedVector.Raw(conditionalProbabilities));
                } else {
                    System.err.println("new: " + id + unambiguousState);
                }
            }
        }
    }

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
                System.err.println("old: " + id + new WrappedVector.Raw(probabilities));
                if (MAS_KILL) {
                    if (count > 10) {
                        System.exit(-1);
                    }
                }
                ++count;
            }

            reconstructedStates[node][j] =  drawChoice(conditionalProbabilities);
        }
//        hookCalculation(tree, parent, node, parentState, state, probabilities, rateCategory);
    }

    private static final boolean MAS_KILL = false;
    private static final boolean MAS_DEBUG = false;
    private int count = 0;

//    private void simulateExternalNodeWithPartials(final NodeRef node, final NodeRef parent) {
//    }
//
//    private void simulateExternalNodeWithStates(final NodeRef node, final NodeRef parent) {
//    }

//    private void simulateInternalNode(final NodeRef node, final NodeRef parent) {
//
//        double[] conditionalProbabilities = tmpConditional;
//        double[] partialLikelihood = getPostOrderPartials(node);
//        double[] probabilities = getTransitionMatrix(node, parent);
//        int[] parentStates = reconstructedStates[parent.getNumber()];
//
//        for (int j = 0; j < patternCount; ++j) {
//
//            int parentState = parentStates[j];
//            int childIndex = j * stateCount;
//
//            int category = siteSpecificRateCategory[j];
//            int matrixOffset = category * stateCount * stateCount;
//            int partialOffset = category * stateCount * patternCount;
//
//            for (int i = 0; i < stateCount; i++) {
//                    conditionalProbabilities[i] = partialLikelihood[partialOffset + childIndex + i]
//                            * probabilities[matrixOffset + parentState * stateCount + i]; // need a transpose here!!!
//            }
//
//            reconstructedStates[node.getNumber()][j] =  drawChoice(conditionalProbabilities);
//        }
////        hookCalculation(tree, parent, node, parentState, state, probabilities, rateCategory);
//    }

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
//            operations[k + 2] = getActiveNodeIndex(op.getLeftChild());   // Node post-order partial
//            operations[k + 3] = getActiveMatrixIndex(op.getLeftChild()); // Node branch info
            operations[k + 4] = op.getLeftChild() < tree.getExternalNodeCount() ? 1 : 0;    // Is node external?

//            operations[k    ] = op.getNodeNumber(); // Parent sample
//            operations[k + 1] = op.getLeftChild();  // Node sample
//            operations[k + 2] = op.getLeftChild() < tree.getExternalNodeCount() ? 1 : 0;    // Is node external?

            k += 5;
        }

        return nodeOperations.size();
    }

//    private static final int OPERATION_LENGTH = 3;

//    abstract int getActiveNodeIndex(int nodeNumber);
//
//    abstract int getActiveMatrixIndex(int nodeNumber);

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
