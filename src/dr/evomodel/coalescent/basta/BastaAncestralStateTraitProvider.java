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
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.HiddenCodons;
import dr.evolution.datatype.HiddenDataType;
import dr.evolution.tree.*;
import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.AncestralStateTraitProvider;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;
import java.util.function.Function;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.BranchIntervalOperation;

/**
 * @author Yucai Shao
 * @author Marc A. Suchard
 */

public class BastaAncestralStateTraitProvider extends AbstractModel implements AncestralStateTraitProvider, Citable {

    private final BastaLikelihood bastaLikelihood;
    private  BastaLikelihoodDelegate likelihoodDelegate;

    private  Tree tree;
    private  PatternList patternList;
    private  SubstitutionModel substitutionModel;
    private  BranchRateModel branchRateModel;
    private  int stateCount;

    private final TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper();

    private  CoalescentIntervalTraversal treeTraversalDelegate;
    public  BestSignalsFromBigFastTreeIntervals treeIntervals;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown;

    private boolean treeIntervalsKnown;
    private boolean transitionMatricesKnown;

    // Ancestral state reconstruction variables
    private  DataType dataType;
    private final String tag;
    private final CodeFormatter formatter;

    // State reconstruction settings
    private  boolean useMAP;
    private  boolean returnMarginalLogLikelihood;
    private  boolean conditionalProbabilitiesInLogSpace;
    private  boolean useOriginalDrawChoice;

    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;
//    private int[][] subIntervalStates;  // [interval][nodeNumber][pattern]
    private Map<Integer, List<Integer>> nodeIntervalMap = new HashMap<>();
    protected boolean areStatesRedrawn = false;
    protected boolean storedAreStatesRedrawn = false;
    protected double jointLogLikelihood;
    private double storedJointLogLikelihood;

    public BastaAncestralStateTraitProvider(String name,
                                            BastaLikelihood bastaLikelihood,
                                            String tag,
                                            boolean useMAP,
                                            boolean returnMarginalLogLikelihood,
                                            boolean conditionalProbabilitiesInLogSpace) {

        super(name);
        this.bastaLikelihood = bastaLikelihood;
        this.tag = tag;

        boolean stripHiddenState = false;
        this.formatter = new CodeFormatter(dataType, stripHiddenState);
    }

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

    public final Model getModel() {
        return this;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

        // Store ancestral state reconstruction information
        if (areStatesRedrawn) {
            for (int i = 0; i < reconstructedStates.length; i++) {
                System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
            }
        }

        storedAreStatesRedrawn = areStatesRedrawn;
    }

    @Override
    protected final void restoreState() {
        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
    }

    @Override
    protected void acceptState() { } // nothing to do

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public MutableTreeModel getTreeModel() {
        return null;
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return null;
    }

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
        return null;
    }

    @Override
    public List<Citation> getCitations() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getCitations();
        } else {
            return new ArrayList<>();
        }
    }

    public int[] getStatesForNode(Tree tree, NodeRef node) {

        if (tree != this.tree) {
            throw new RuntimeException("Can only reconstruct states on tree given to constructor");
        }

        bastaLikelihood.getLogLikelihood();

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNumber()];
    }


//    private void mapNodeToSubIntervals() {
//        nodeIntervalMap.clear();
//        CoalescentIntervalTraversal traversal = getTraversalDelegate();
//        traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();
//
//        List<BranchIntervalOperation> branchIntervalOps = traversal.getBranchIntervalOperations();
//
//        for (BranchIntervalOperation op : branchIntervalOps) {
//            int nodeNumber = op.inputBuffer1 % tree.getNodeCount();
//            if (!nodeIntervalMap.containsKey(nodeNumber)) {
//                nodeIntervalMap.put(nodeNumber, new ArrayList<>());
//            }
//            nodeIntervalMap.get(nodeNumber).add(op.intervalNumber);
//
//            if (op.inputBuffer2 > 0) {
//                int nodeNumber2 = op.inputBuffer2 % tree.getNodeCount();
//                if (!nodeIntervalMap.containsKey(nodeNumber2)) {
//                    nodeIntervalMap.put(nodeNumber2, new ArrayList<>());
//                }
//                nodeIntervalMap.get(nodeNumber2).add(op.intervalNumber);
//            }
//        }
//
//        // Sort interval lists for each node to ensure correct processing order
//        for (List<Integer> intervals : nodeIntervalMap.values()) {
//            Collections.sort(intervals);
//        }
//    }

    private void redrawAncestralStates(SamplingTraversalMethod traversalMethod) {

        if (traversalMethod == SamplingTraversalMethod.PRE_ORDER) {
            traverseSampleByNodes();
        } else {
            throw new IllegalArgumentException("Invalid traversal method: " + traversalMethod);
        }

        areStatesRedrawn = true;
    }

    public void redrawAncestralStates() {
        redrawAncestralStates(SamplingTraversalMethod.PRE_ORDER); // Use original method by default
    }

    enum SamplingTraversalMethod {
        LEVEL_ORDER,
        PRE_ORDER
    }

//    private void traverseSampleByCoalescentIntervals() {
//
//        CoalescentIntervalTraversal traversal = getTraversalDelegate();
//        traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();
//        List<BranchIntervalOperation> ops = traversal.getBranchIntervalOperations();
//        List<TransitionMatrixOperation> matrixOps = traversal.getMatrixOperations();
//        List<Integer> intervalStarts = traversal.getIntervalStarts();
//        int rootBuffer = 0;
//        BranchIntervalOperation.initializeMap(tree, likelihoodDelegate.getMaxNumberOfCoalescentIntervals());
//        Map<Integer, Integer> bufferToNodeMap = new HashMap<>();
//        BranchIntervalOperation lastOp = ops.get(ops.size() - 1);
//
//        for (BranchIntervalOperation op : ops) {
//            int input1NodeNumber = op.inputBuffer1 % tree.getNodeCount();
//            int input2NodeNumber = op.inputBuffer2 % tree.getNodeCount();
//            int outputNodeNumber = op.outputBuffer % tree.getNodeCount();
//            op.transform();
//            if (!bufferToNodeMap.containsKey(op.inputBuffer1)) {
//                bufferToNodeMap.put(op.inputBuffer1, input1NodeNumber);
//            }
//
//            if (op.inputBuffer2 >= 0) {
//
//                if (!bufferToNodeMap.containsKey(op.inputBuffer2)) {
//                    bufferToNodeMap.put(op.inputBuffer2, input2NodeNumber);
//                }
//            }
//
//            if (!bufferToNodeMap.containsKey(op.outputBuffer)) {
//                bufferToNodeMap.put(op.outputBuffer, outputNodeNumber);
//            }
//
//            if (op == lastOp) {
//                rootBuffer = op.outputBuffer;
//            }
//
//        }
//
//        Map<Integer, TransitionMatrixOperation> matrixOpMap = new HashMap<>();
//        for (TransitionMatrixOperation op : matrixOps) {
//            matrixOpMap.put(op.outputBuffer, op);
//        }
//
//        NodeRef rootNode = tree.getRoot();
//        sampleRootState(rootNode, rootBuffer, true);
//
//        for (int interval = intervalStarts.size() - 2; interval >= 0; interval--) {
//            int start = intervalStarts.get(interval);
//            int end = intervalStarts.get(interval + 1);
//
//            // Process operations within this interval
//            for (int i = start; i < end; i++) {
//                BranchIntervalOperation op = ops.get(i);
//
//                if (op.inputBuffer1 >= 0) {
//                    processChildBuffer(op, op.inputBuffer1, op.inputMatrix1,
//                            bufferToNodeMap, matrixOpMap);
//                }
//
//                if (op.inputBuffer2 >= 0) {
//                    processChildBuffer(op, op.inputBuffer2, op.inputMatrix2,
//                            bufferToNodeMap, matrixOpMap);
//                }
//            }
//        }
//    }


//    private void processChildBuffer(BranchIntervalOperation op, int inputBuffer, int matrixBuffer,
//                                    Map<Integer, Integer> bufferToNodeMap,
//                                    Map<Integer, TransitionMatrixOperation> matrixOpMap) {
//
//        Integer childNodeNumber = bufferToNodeMap.get(inputBuffer);
//
//        NodeRef childNode = tree.getNode(childNodeNumber);
//
//
//        double[] transitionMatrix = new double[stateCount * stateCount];
//        if (matrixBuffer >= 0 && matrixOpMap.containsKey(matrixBuffer)) {
//            TransitionMatrixOperation matrixOp = matrixOpMap.get(matrixBuffer);
//            getMatrix(likelihoodDelegate, matrixOp.outputBuffer, transitionMatrix);
//        } else {
//            throw new RuntimeException("No transition matrix matches");
//        }
//
//        if (tree.isExternal(childNode) && inputBuffer < tree.getExternalNodeCount()) {
//            double[] partials = new double[stateCount];
//            likelihoodDelegate.getPartials(inputBuffer, partials);
//
//            for (int j = 0; j < getPatternCount(); j++) {
//                boolean isAmbiguous;
//                int unambiguousState = -1;
//                int statesWithProbability = 0;
//
//                for (int k = 0; k < stateCount; k++) {
//                    if (partials[k] > 0) {
//                        statesWithProbability++;
//                        unambiguousState = k;
//                    }
//                }
//
//                isAmbiguous = (statesWithProbability != 1);
//
//                if (isAmbiguous) {
//                    sampleStateForSubInterval(childNode, childNodeNumber,
//                            inputBuffer, transitionMatrix);
//                    reconstructedStates[childNodeNumber][j] = subIntervalStates[childNodeNumber][j];
//                } else {
//                    subIntervalStates[childNodeNumber][j] = unambiguousState;
//                    reconstructedStates[childNodeNumber][j] = unambiguousState;
//                }
//            }
//            return;
//        }
//
//        sampleStateForSubInterval(childNode, childNodeNumber,
//                inputBuffer, transitionMatrix);
//
//
//        List<Integer> nodeIntervals = nodeIntervalMap.get(childNodeNumber);
//        if (nodeIntervals != null) {
//            int lastInterval = nodeIntervals.get(nodeIntervals.size() - 1);
//            if (op.intervalNumber == lastInterval) {
//                System.arraycopy(
//                        subIntervalStates[childNodeNumber], 0,
//                        reconstructedStates[childNodeNumber], 0,
//                        getPatternCount()
//                );
//            }
//        }
//    }

    private void sampleRootState(NodeRef root, int rootBuffer, boolean copyToIntervals) {
        int nodeNum = root.getNumber();

        double[] partials = new double[stateCount];
        likelihoodDelegate.getPartials(rootBuffer, partials);

        double[] conditionalProbabilities = new double[stateCount];
        double[] frequencies = substitutionModel.getFrequencyModel().getFrequencies();

        for (int j = 0; j < bastaLikelihood.getPatternCount(); j++) {
            for (int i = 0; i < stateCount; i++) {
                if (conditionalProbabilitiesInLogSpace) {
                    conditionalProbabilities[i] = Math.log(partials[i]) + Math.log(frequencies[i]);
                } else {
                    conditionalProbabilities[i] = partials[i] * frequencies[i];
                }
            }

            reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);

//            if (copyToIntervals) {
//                    subIntervalStates[nodeNum][j] = reconstructedStates[nodeNum][j];
//            }

            if (!returnMarginalLogLikelihood) {
                jointLogLikelihood += Math.log(frequencies[reconstructedStates[nodeNum][j]]);
            }
        }
    }


//    private void sampleStateForSubInterval(NodeRef node, int nodeNumber,
//                                           int bufferIndex, double[] transitionMatrix) {
//
//        int patternCount = bastaLikelihood.getPatternCount();
//
//        double[] nodeLikelihoods = new double[stateCount];
//        likelihoodDelegate.getPartials(bufferIndex, nodeLikelihoods);
//
//        int[] parentState = null;
//        if (subIntervalStates[nodeNumber][0] != -1) {
//            parentState = subIntervalStates[nodeNumber];
//        }
//
//        if (parentState == null) {
//            NodeRef parent = tree.getParent(node);
//            if (parent == null) {
//                return;
//            }
//
//            int parentNumber = parent.getNumber();
//            parentState = subIntervalStates[parentNumber];
//        }
//
//        double[] conditionalProbabilities = new double[stateCount];
//
//        for (int j = 0; j < patternCount; j++) {
//            int parentIndex = parentState[j];
//            for (int i = 0; i < stateCount; i++) {
//                if (conditionalProbabilitiesInLogSpace) {
//                    conditionalProbabilities[i] = Math.log(transitionMatrix[i * stateCount + parentIndex]) +
//                            Math.log(nodeLikelihoods[i]);
//                } else {
//                    conditionalProbabilities[i] = transitionMatrix[i * stateCount + parentIndex] *
//                            nodeLikelihoods[i];
//                }
//            }
//
//            subIntervalStates[nodeNumber][j] = drawChoice(conditionalProbabilities);
//            if (!returnMarginalLogLikelihood) {
//                double contrib = transitionMatrix[subIntervalStates[nodeNumber][j] * stateCount + parentIndex];
//                jointLogLikelihood += Math.log(contrib);
//            }
//        }
//    }

//    private void getMatrix(BastaLikelihoodDelegate delegate, int bufferIndex, double[] probabilities) {
//        double[] matrices = new double[stateCount * stateCount];
//        delegate.getTransitionMatrices(bufferIndex, matrices);
//        System.arraycopy(
//                matrices,
//                0,
//                probabilities, 0,
//                stateCount * stateCount
//        );
//    }


    // TODO this mapping should be computed during post-order update, not here
    private Map<Integer, Integer> getNodeToBufferMap() {

        CoalescentIntervalTraversal traversal = bastaLikelihood.getTraversalDelegate();
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

        return nodeToBufferMap;
    }

    private void traverseSampleByNodes() {

        Map<Integer, Integer> nodeToBufferMap = getNodeToBufferMap();

        NodeRef rootNode = tree.getRoot();
        int rootBuffer = nodeToBufferMap.getOrDefault(rootNode.getNumber(), rootNode.getNumber());

        sampleRootState(rootNode, rootBuffer, false);

        traverseNodesPreOrderDirect(rootNode, nodeToBufferMap);
    }

    private void traverseNodesPreOrderDirect(NodeRef node, Map<Integer, Integer> nodeToBufferMap) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);

            if (tree.isExternal(child)) {
                int childNum = child.getNumber();
                int parentNum = node.getNumber();

                double[] partials = new double[stateCount];
                int childBuffer = nodeToBufferMap.getOrDefault(childNum, childNum);
                likelihoodDelegate.getPartials(childBuffer, partials);

                for (int j = 0; j < bastaLikelihood.getPatternCount(); j++) {
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
        int patternCount = bastaLikelihood.getPatternCount();

        double[] childPartials = new double[stateCount];
        likelihoodDelegate.getPartials(bufferIndex, childPartials);

        double[] conditionalProbabilities = new double[stateCount];

        for (int j = 0; j < patternCount; j++) {
            int parentState = reconstructedStates[parentNumber][j];
            int parentIndex = parentState;

            for (int i = 0; i < stateCount; i++) {
                double transProb = transitionMatrix[i * stateCount + parentIndex];
                double likelihood = childPartials[i];

                conditionalProbabilities[i] = conditionalProbabilitiesInLogSpace ?
                        Math.log(transProb) + Math.log(likelihood) :
                        transProb * likelihood;
            }

            reconstructedStates[childNumber][j] = drawChoice(conditionalProbabilities);

            if (!returnMarginalLogLikelihood) {
                double contrib = transitionMatrix[reconstructedStates[childNumber][j] * stateCount + parentIndex];
                jointLogLikelihood += Math.log(contrib);
            }
        }
    }

    private int drawChoice(double[] measure) {
        if (useOriginalDrawChoice) {
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
        } else {
            double sum = 0.0;
            for (int i = 0; i < measure.length; i++) {
                sum += measure[i];
            }

            if (sum < 0.0000000001) {
                //System.out.println("Warning: Probability sum is extremely small: " + sum);
                double max = measure[0];
                int choice = 0;
                for (int i = 1; i < measure.length; i++) {
                    if (measure[i] > max) {
                        max = measure[i];
                        choice = i;
                    }
                }
                return choice;
            }

            double[] normalizedMeasure = new double[measure.length];
            for (int i = 0; i < measure.length; i++) {
                normalizedMeasure[i] = measure[i] / sum;
            }

            double x = MathUtils.nextDouble();
            sum = 0.0;
            for (int i = 0; i < normalizedMeasure.length; i++) {
                sum += normalizedMeasure[i];
                if (x < sum) {
                    return i;
                }
            }

            return normalizedMeasure.length - 1;
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

    private static class CodeFormatter {

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
}