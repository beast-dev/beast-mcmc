/*
 * AncestralStateBeagleTreeLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treelikelihood;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.Map;
import java.util.Set;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */

public class AncestralStateBeagleTreeLikelihood extends BeagleTreeLikelihood implements TreeTraitProvider, AncestralStateTraitProvider {

//    public AncestralStateBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
//                                              BranchSubstitutionModel branchSubstitutionModel, SiteRateModel siteRateModel,
//                                              BranchRateModel branchRateModel, boolean useAmbiguities,
//                                              PartialsRescalingScheme scalingScheme,
//                                              DataType dataType,
//                                              String tag,
//                                              SubstitutionModel substModel) {
//        this(patternList, treeModel, branchSubstitutionModel, siteRateModel, branchRateModel, useAmbiguities, scalingScheme,
//                dataType, tag, substModel, false, true);
//    }

    public AncestralStateBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                              BranchModel branchModel,
                                              SiteRateModel siteRateModel,
                                              BranchRateModel branchRateModel,
                                              TipStatesModel tipStatesModel,
                                              boolean useAmbiguities,
                                              PartialsRescalingScheme scalingScheme,
                                              boolean delayRescalingUntilUnderflow,
                                              Map<Set<String>, Parameter> partialsRestrictions,
                                              final DataType dataType,
                                              final String tag,
//                                              SubstitutionModel substModel,
                                              boolean useMAP,
                                              boolean returnML) {

        super(patternList, treeModel, branchModel, siteRateModel, branchRateModel, tipStatesModel, useAmbiguities, scalingScheme, delayRescalingUntilUnderflow,
                partialsRestrictions);

        this.dataType = dataType;
//        this.tag = tag;

        probabilities = new double[stateCount * stateCount * categoryCount];
        partials = new double[stateCount * patternCount * categoryCount];
//        rootPartials = new double[stateCount*patternCount];
//        cumulativeScaleBuffers = new int[nodeCount][];
//        scaleBufferIndex = getScaleBufferCount() - 1;

        // Save tip states locally so these do not need to be transfers back

        if (useAmbiguities()) {
            tipPartials = new double[tipCount][];
        } else {
            tipStates = new int[tipCount][];
        }

        for (int i = 0; i < tipCount; i++) {
            // Find the id of tip i in the patternList
            String id = treeModel.getTaxonId(i);
            int index = patternList.getTaxonIndex(id);
            if (useAmbiguities()) {
                tipPartials[i] = getPartials(patternList, index);
            } else {
                tipStates[i] = getStates(patternList, index);
            }
        }

        reconstructedStates = new int[treeModel.getNodeCount()][patternCount];
        storedReconstructedStates = new int[treeModel.getNodeCount()][patternCount];

        this.useMAP = useMAP;
        this.returnMarginalLogLikelihood = returnML;

        treeTraits.addTrait(new TreeTrait.IA() {
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
                return formattedState(getStatesForNode(tree, node), dataType);
            }
        });

    }

    private double[] getPartials(PatternList patternList, int sequenceIndex) {
        double[] partials = new double[patternCount * stateCount];

        boolean[] stateSet;

        int v = 0;
        for (int i = 0; i < patternCount; i++) {


            if (patternList instanceof UncertainSiteList) {
                ((UncertainSiteList) patternList).fillPartials(sequenceIndex, i, partials, v);
                v += stateCount;
                // TODO Add this functionality to SimpleSiteList to avoid if statement here
            } else {

                int state = patternList.getPatternState(sequenceIndex, i);
                stateSet = dataType.getStateSet(state);

                for (int j = 0; j < stateCount; j++) {
                    if (stateSet[j]) {
                        partials[v] = 1.0;
                    } else {
                        partials[v] = 0.0;
                    }
                    v++;
                }
            }
        }  // TODO Note code duplication with BTL, refactor when debugged

        return partials;
    }

    private int[] getStates(PatternList patternList,
                            int sequenceIndex) {

        int[] states = new int[patternCount];
        for (int i = 0; i < patternCount; i++) {
            states[i] = patternList.getPatternState(sequenceIndex, i);
        }
        return states;
    }

    public BranchModel getBranchModel() {
        return branchModel;
    }

    protected Helper treeTraits = new Helper();

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);
        fireModelChanged(model);
    }

    public int[] getStatesForNode(Tree tree, NodeRef node) {
        if (tree != treeModel) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
            calculateLogLikelihood();
            likelihoodKnown = true;
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNumber()];
    }

    @Override
    protected int getScaleBufferCount() {
        return internalNodeCount + 2;
    }

    private int drawChoice(double[] measure) {
        if (useMAP) {
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

    public void makeDirty() {
        super.makeDirty();
        areStatesRedrawn = false;
    }

    public void redrawAncestralStates() {
        // Sample states
        jointLogLikelihood = 0;
        traverseSample(treeModel, treeModel.getRoot(), null, null);
        areStatesRedrawn = true;
    }

//    protected double calculateLogLikelihood() {
//        areStatesRedrawn = false;
//        return super.calculateLogLikelihood();
//    }

    protected double calculateLogLikelihood() {
        areStatesRedrawn = false;
        double marginalLogLikelihood = super.calculateLogLikelihood();
        if (returnMarginalLogLikelihood) {
            return marginalLogLikelihood;
        }
        // redraw states and return joint density of drawn states
        redrawAncestralStates();
        return jointLogLikelihood;
    }

    public String formattedState(int[] state) {
        return formattedState(state, dataType);
    }

    private static String formattedState(int[] state, DataType dataType) {
        StringBuffer sb = new StringBuffer();
        sb.append("\"");
        if (dataType instanceof GeneralDataType) {
            boolean first = true;
            for (int i : state) {
                if (!first) {
                    sb.append(" ");
                } else {
                    first = false;
                }

                sb.append(dataType.getCode(i));
            }

        } else {
            for (int i : state) {
                if (dataType.getClass().equals(Codons.class)) {
                    sb.append(dataType.getTriplet(i));
                } else {
                    sb.append(dataType.getChar(i));
                }
            }
        }
        sb.append("\"");
        return sb.toString();
    }

    protected void getMatrix(int branchIndex, double[] probabilities) {
        beagle.getTransitionMatrix(substitutionModelDelegate.getMatrixIndex(branchIndex), probabilities);
        // NB: It may be faster to compute matrices in BEAST via substitutionModel
    }

    public void setTipStates(int tipNum, int[] states) {
        System.arraycopy(states, 0, tipStates[tipNum], 0, states.length);
        beagle.setTipStates(tipNum, states);
        makeDirty();
    }

//    public void getTipPartials(int tipNum, double[] partials) {
//        System.arraycopy(tipPartials[tipNum], 0, partials, 0, partials.length);
//    }

    public void getTipStates(int tipNum, int[] states) {
        // Saved locally to reduce BEAGLE library access
        System.arraycopy(tipStates[tipNum], 0, states, 0, states.length);
    }

//    public int traverseCollectScaleBuffers(TreeModel tree, NodeRef node) {
//
//        if (true) // Currently do nothing
//            return 0;
//
//            return 0;
//
//        int nodeNum = node.getNumber();
//
//        NodeRef child0 = tree.getChild(node,0);
//        NodeRef child1 = tree.getChild(node,1);
//
//        int len0 = traverseCollectScaleBuffers(tree,child0);
//        int len1 = traverseCollectScaleBuffers(tree,child1);
//        int thisLen = len0 + len1 + 1;
//        int offset = 0;
//
//        int[] scaleBuffer = new int[thisLen];
//        if (len0 > 0) {
//            System.arraycopy(cumulativeScaleBuffers[child0.getNumber()],0,scaleBuffer,offset,len0);
//            offset += len0;
//        }
//        if (len1 > 0) {
//            System.arraycopy(cumulativeScaleBuffers[child1.getNumber()],0,scaleBuffer,offset,len1);
//            offset += len1;
//        }
//        scaleBuffer[offset] = scaleBufferHelper.getOffsetIndex(nodeNum - tipCount);
//        cumulativeScaleBuffers[nodeNum] = scaleBuffer;
//
//        return thisLen;
//    }

    public void storeState() {

        super.storeState();

        if (areStatesRedrawn) {
            for (int i = 0; i < reconstructedStates.length; i++) {
                System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
            }
        }
        // TODO MAS: I do not understand why these are NOT necessary

        storedAreStatesRedrawn = areStatesRedrawn;
        storedJointLogLikelihood = jointLogLikelihood;
    }

    public void restoreState() {

        super.restoreState();

        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
        jointLogLikelihood = storedJointLogLikelihood;
    }

    public void traverseSample(TreeModel tree, NodeRef node, int[] parentState, int[] rateCategory) {

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!tree.isExternal(node)) {

            if (parent == null) {

                // This is the root node
                getPartials(nodeNum, partials);

                boolean sampleCategory = categoryCount > 1;
                double[] posteriorWeightedCategory = null;
                double[] priorWeightedCategory = null;

                if (sampleCategory) {
                    rateCategory = new int[patternCount];
                    posteriorWeightedCategory = new double[categoryCount];
                    priorWeightedCategory = siteRateModel.getCategoryProportions();
                }

                for (int j = 0; j < patternCount; j++) {

                    // Sample across-site-rate-variation, if it exists
                    if (sampleCategory) {
                        for (int r = 0; r < categoryCount; r++) {
                            posteriorWeightedCategory[r] = 0;
                            for (int k = 0; k < stateCount; k++) {
                                posteriorWeightedCategory[r] += partials[r * stateCount * patternCount +
                                        j * stateCount + k];
                            }
                            posteriorWeightedCategory[r] *= priorWeightedCategory[r];
                        }
                        rateCategory[j] = drawChoice(posteriorWeightedCategory);
                    }

                    // Sample root character state
                    int partialsIndex = (rateCategory == null ? 0 : rateCategory[j]) * stateCount * patternCount;
                    System.arraycopy(partials, partialsIndex + j * stateCount, conditionalProbabilities, 0, stateCount);

                    double[] frequencies = substitutionModelDelegate.getRootStateFrequencies(); // TODO May have more than one set of frequencies
                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] *= frequencies[i];
                    }
                    try {
                        state[j] = drawChoice(conditionalProbabilities);
                    } catch (Error e) {
                        System.err.println(e.toString());
                        System.err.println("Please report error to Marc");
                        state[j] = 0;
                    }
                    reconstructedStates[nodeNum][j] = state[j];

                    if (!returnMarginalLogLikelihood) {
                        jointLogLikelihood += Math.log(frequencies[state[j]]);
                    }
                }

                if (sampleCategory) {
                    if (this.rateCategory == null) {
                        this.rateCategory = new int[patternCount];
                    }
                    System.arraycopy(rateCategory, 0, this.rateCategory, 0, patternCount);
                }

            } else {

                // This is an internal node, but not the root
                double[] partialLikelihood = new double[stateCount * patternCount * categoryCount];
                getPartials(nodeNum, partialLikelihood);

                // Sibon says that this actually works now
//                if (categoryCount > 1)
//                    throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");

                getMatrix(nodeNum, probabilities);

                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    int category = rateCategory == null ? 0 : rateCategory[j];
                    int matrixIndex = category * stateCount * stateCount;
                    int partialIndex = category * stateCount * patternCount;

                    for (int i = 0; i < stateCount; i++)
                        conditionalProbabilities[i] = partialLikelihood[partialIndex + childIndex + i]
                                * probabilities[matrixIndex + parentIndex + i];

                    state[j] = drawChoice(conditionalProbabilities);
                    reconstructedStates[nodeNum][j] = state[j];

                    if (!returnMarginalLogLikelihood) {
                        double contrib = probabilities[parentIndex + state[j]];
                        jointLogLikelihood += Math.log(contrib);
                    }
                }

                hookCalculation(tree, parent, node, parentState, state, probabilities, rateCategory);
            }

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            traverseSample(tree, child1, state, rateCategory);

            NodeRef child2 = tree.getChild(node, 1);
            traverseSample(tree, child2, state, rateCategory);
        } else {

            // This is an external leaf

            if (useAmbiguities()) {

                getMatrix(nodeNum, probabilities);
                double[] partials = tipPartials[nodeNum];

                for (int j = 0; j < patternCount; j++) {
                    final int parentIndex = parentState[j] * stateCount;
                    int category = rateCategory == null ? 0 : rateCategory[j];
                    int matrixIndex = category * stateCount * stateCount;

                    System.arraycopy(probabilities, parentIndex + matrixIndex, conditionalProbabilities, 0, stateCount);
                    for (int k = 0; k < stateCount; ++k) {
                        conditionalProbabilities[k] *= partials[j * stateCount + k];
                    }
                    reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);

                    if (!returnMarginalLogLikelihood) {
                        double contrib = probabilities[parentIndex + reconstructedStates[nodeNum][j]];
                        jointLogLikelihood += Math.log(contrib);
                    }
                }

            } else {

                getTipStates(nodeNum, reconstructedStates[nodeNum]);

                // Check for ambiguity codes and sample them

                for (int j = 0; j < patternCount; j++) {

                    final int thisState = reconstructedStates[nodeNum][j];

                    if (dataType.isAmbiguousState(thisState)) {

                        final int parentIndex = parentState[j] * stateCount;
                        int category = rateCategory == null ? 0 : rateCategory[j];
                        int matrixIndex = category * stateCount * stateCount;

                        getMatrix(nodeNum, probabilities);
                        System.arraycopy(probabilities, parentIndex + matrixIndex, conditionalProbabilities, 0, stateCount);

                        if (useAmbiguities && !dataType.isUnknownState(thisState)) { // Not completely unknown
                            boolean[] stateSet = dataType.getStateSet(thisState);

                            for (int k = 0; k < stateCount; k++) {
                                if (!stateSet[k]) {
                                    conditionalProbabilities[k] = 0.0;
                                }
                            }
                        }
                        reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
                    }

                    if (!returnMarginalLogLikelihood) {
                        final int parentIndex = parentState[j] * stateCount;
                        getMatrix(nodeNum, probabilities);
                        if (!returnMarginalLogLikelihood) {
                            double contrib = probabilities[parentIndex + reconstructedStates[nodeNum][j]];
                            jointLogLikelihood += Math.log(contrib);
                        }
                    }
                }
            }

            hookCalculation(tree, parent, node, parentState, reconstructedStates[nodeNum], null, rateCategory);
        }
    }

    protected void hookCalculation(Tree tree, NodeRef parentNode, NodeRef childNode,
                                   int[] parentStates, int[] childStates,
                                   double[] probabilities, int[] rateCategory) {
        // Do nothing
    }

    private final DataType dataType;
    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;

    //    private final String tag;
    protected boolean areStatesRedrawn = false;
    protected boolean storedAreStatesRedrawn = false;

    private boolean useMAP = false;
    private boolean returnMarginalLogLikelihood = true;

    private double jointLogLikelihood;
    private double storedJointLogLikelihood;

    private int[][] tipStates;
    private double[][] tipPartials;

    private double[] probabilities;
    private double[] partials;

    protected int[] rateCategory = null;
//    private double[] rootPartials;
//    private int[][] cumulativeScaleBuffers;
//    private int scaleBufferIndex;

}
