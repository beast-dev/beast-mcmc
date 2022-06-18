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

import dr.evolution.datatype.*;
import dr.evolution.tree.*;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RateEpochBranchRateModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.markovjumps.MarkovJumpsCore;
import dr.math.MathUtils;

import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import java.util.*;

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

    public AncestralStateBeagleTreeLikelihood(PatternList patternList, MutableTreeModel treeModel,
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
                                              boolean returnML,
                                              boolean conditionalProbabilitiesInLogSpace) {

        super(patternList, treeModel, branchModel, siteRateModel, branchRateModel, tipStatesModel, useAmbiguities, scalingScheme, delayRescalingUntilUnderflow,
                partialsRestrictions);
        this.conditionalProbabilitiesInLogSpace = conditionalProbabilitiesInLogSpace;
        this.dataType = dataType;
//        this.tag = tag;

        probabilities = new double[stateCount * stateCount * categoryCount];
        partials = new double[stateCount * patternCount * categoryCount];
        
        probabilitiesAlongBranch = new ArrayList<double[]>();
        probabilitiesConvolved = new ArrayList<double[]>();
        combinedWeights = new ArrayList<Double>();
        combinedMatrixOrder = new ArrayList<Integer>();
        combinedRates = new ArrayList<Double>();
        
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

        boolean stripHiddenState = false; // TODO Pass as option
        this.formatter = new CodeFormatter(dataType, stripHiddenState);

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
                return formattedState(getStatesForNode(tree, node), formatter);
            }
        });

    }
    public AncestralStateBeagleTreeLikelihood(PatternList patternList, MutableTreeModel treeModel,
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
                                              boolean returnML){
        this(patternList, treeModel, branchModel, siteRateModel, branchRateModel, tipStatesModel, useAmbiguities,
                scalingScheme, delayRescalingUntilUnderflow, partialsRestrictions, dataType, tag, useMAP, returnML, false);
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

    protected int drawChoice(double[] measure) {
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
            if(conditionalProbabilitiesInLogSpace){
                return MathUtils.randomChoiceLogPDF(measure);
            }
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
        return formattedState(state, formatter);
    }

    private static String formattedState(int[] state, CodeFormatter formatter) {
        StringBuffer sb = new StringBuffer();
        sb.append("\"");
        formatter.reset();
        for (int i : state) {
            sb.append(formatter.getCodeString(i));
        }
//        if (dataType instanceof GeneralDataType) {
//            boolean first = true;
//            for (int i : state) {
//                if (!first) {
//                    sb.append(" ");
//                } else {
//                    first = false;
//                }
//
//                sb.append(dataType.getCode(i));
//            }
//
//        } else {
//            for (int i : state) {
//                if (dataType instanceof Codons) {
//                    sb.append(dataType.getTriplet(i));
//                } else {
//                    sb.append(dataType.getChar(i));
//                }
//            }
//        }
        sb.append("\"");
        return sb.toString();
    }

    protected void getMatrix(int branchIndex, double[] probabilities) {
        beagle.getTransitionMatrix(substitutionModelDelegate.getMatrixIndex(branchIndex), probabilities);
        // NB: It may be faster to compute matrices in BEAST via substitutionModel
    }
    
    protected void getTransitionProbabilityMatrices(Tree tree, NodeRef childNode) {
    
        NodeRef parentNode = tree.getParent(childNode);

        final double parentTime = tree.getNodeHeight(parentNode);
        final double childTime = tree.getNodeHeight(childNode);
        final double branchTime = parentTime - childTime;
        final double branchRate = branchRateModel.getBranchRate(tree, childNode);
        
        // fetch the matrix and/or rate epochal layout along the branch
        // combinedWeights contains the (absolute instead of relative) duration of each piece
        combinedWeights.clear();
        combinedMatrixOrder.clear();
        combinedRates.clear();
        
        BranchModel.Mapping matrixMapping = branchModel.getBranchModelMapping(childNode);
        int[] matrixOrder = matrixMapping.getOrder();
        double[] matrixWeights = matrixMapping.getWeights();
        int nmatrices = matrixOrder.length;
        
        if (!(branchModel instanceof EpochBranchModel)) {
            for (int j = 0; j < nmatrices; j++) {
                matrixWeights[j] *= branchTime;
            }
        }
        
        BranchRateModel.Mapping rateMapping = branchRateModel.getBranchRateModelMapping(tree, childNode);
        double[] rates = rateMapping.getRates();
        double[] rateWeights = rateMapping.getWeights();
        int nrates = rates.length;
        
        if (!(branchRateModel instanceof RateEpochBranchRateModel)) {
            for (int j = 0; j < nrates; j++) {
                rateWeights[j] *= branchTime;
            }
        }
        
        // generate cumulative sum vector for matrices
        double[] matrixWeightsCumsum = new double[nmatrices + 1];
        for (int j = 0; j < nmatrices - 1; j++) {
            matrixWeightsCumsum[j + 1] = matrixWeightsCumsum[j] + matrixWeights[j];
        }
        
        // generate cumulative sum vector for rates
        double[] rateWeightsCumsum = new double[nrates + 1];
        for (int j = 0; j < nrates - 1; j++) {
            rateWeightsCumsum[j + 1] = rateWeightsCumsum[j] + rateWeights[j];
        }
        
        // make sure that the last time in the matrix and rate weights cumsum vectors match (which should equal the branch duration)
        matrixWeightsCumsum[nmatrices] = branchTime;
        rateWeightsCumsum[nrates] = branchTime;

        if (nmatrices == 1 && nrates == 1) {
            
            combinedWeights.add(matrixWeights[0]);
            combinedMatrixOrder.add(matrixOrder[0]);
            combinedRates.add(branchRate);
            
        } else if (nmatrices > 1 && nrates == 1) {
            
            for (double w : matrixWeights) {
                combinedWeights.add(w);
            }
            for (int i : matrixOrder) {
                combinedMatrixOrder.add(i);
            }
            for (int j = 0; j < nmatrices; j++) {
                combinedRates.add(branchRate);
            }
            
        } else if (nmatrices == 1 && nrates > 1) {
            
            for (double w : rateWeights) {
                combinedWeights.add(w);
            }
            for (int j = 0; j < nrates; j++) {
                combinedMatrixOrder.add(matrixOrder[0]);
            }
            for (double r : rates) {
                combinedRates.add(r);
            }
            
        } else {
            
            int matrixId = 0;
            int rateId = 0;
            double lastCumsum = 0.0;
            
            while (matrixId < nmatrices || rateId < nrates) {
                combinedMatrixOrder.add(matrixOrder[matrixId]);
                combinedRates.add(rates[rateId]);
                    
                if (matrixWeightsCumsum[matrixId + 1] < rateWeightsCumsum[rateId + 1]) {
                    combinedWeights.add(matrixWeightsCumsum[matrixId + 1] - lastCumsum);
                    lastCumsum = matrixWeightsCumsum[matrixId + 1];
                    matrixId++;
                } else if (matrixWeightsCumsum[matrixId + 1] > rateWeightsCumsum[rateId + 1]) {
                    combinedWeights.add(rateWeightsCumsum[rateId + 1] - lastCumsum);
                    lastCumsum = rateWeightsCumsum[rateId + 1];
                    rateId++;
                } else {
                    combinedWeights.add(matrixWeightsCumsum[matrixId + 1] - lastCumsum);
                    lastCumsum = matrixWeightsCumsum[matrixId + 1];
                    matrixId++;
                    rateId++;
                }
            }
        }
        
        int npieces = combinedWeights.size();
        
        // compute transition probability matrix for each piece
        probabilitiesAlongBranch.clear();
        for (int j = 0; j < npieces; j++) {
            probabilitiesAlongBranch.add(new double[stateCount * stateCount * categoryCount]);
        
            for (int k = 0; k < categoryCount; k++) {
                final double edgeLength = combinedWeights.get(j) * combinedRates.get(j) * siteRateModel.getRateForCategory(k);
                double[] probabilitiesTmp = new double[stateCount * stateCount];
                substitutionModelDelegate.getSubstitutionModel(combinedMatrixOrder.get(j)).getTransitionProbabilities(edgeLength, probabilitiesTmp);
                System.arraycopy(probabilitiesTmp, 0, probabilitiesAlongBranch.get(j), k * stateCount * stateCount, stateCount * stateCount);
            }
        }
        
        // convolve transition probability matrix
        probabilitiesConvolved.clear();        
        for (int j = npieces - 2; j >= 0; j--) {
            probabilitiesConvolved.add(0, new double[stateCount * stateCount * categoryCount]);
            
            for (int k = 0; k < categoryCount; k++) {
                double[] probabilitiesTmp0 = new double[stateCount * stateCount];
                double[] probabilitiesTmp1 = new double[stateCount * stateCount];
                double[] probabilitiesTmp2 = new double[stateCount * stateCount];
                
                System.arraycopy(probabilitiesAlongBranch.get(j), k * stateCount * stateCount, probabilitiesTmp1, 0, stateCount * stateCount);                
                if (j == npieces - 2) {
                    System.arraycopy(probabilitiesAlongBranch.get(j + 1), k * stateCount * stateCount, probabilitiesTmp2, 0, stateCount * stateCount);
                } else {
                    System.arraycopy(probabilitiesConvolved.get(1), k * stateCount * stateCount, probabilitiesTmp2, 0, stateCount * stateCount);
                }
                
                MarkovJumpsCore.matrixMultiply(probabilitiesTmp1, probabilitiesTmp2, stateCount, probabilitiesTmp0);
                System.arraycopy(probabilitiesTmp0, 0, probabilitiesConvolved.get(0), k * stateCount * stateCount, stateCount * stateCount);
            }
        }
    }
    
    protected void getTransitionProbabilityMatrix(Tree tree, NodeRef childNode, double[] probabilities, boolean update) {
        
        if (update) {
            getTransitionProbabilityMatrices(tree, childNode);
        }
        
        if (probabilitiesConvolved.size() > 0) {
            System.arraycopy(probabilitiesConvolved.get(0), 0, probabilities, 0, categoryCount * stateCount * stateCount);
        } else {
            System.arraycopy(probabilitiesAlongBranch.get(0), 0, probabilities, 0, categoryCount * stateCount * stateCount);
        }
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

    public void traverseSample(Tree tree, NodeRef node, int[] parentState, int[] rateCategory) {

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
                    int partialsIndex = (rateCategory == null ? 0 : rateCategory[j]) * stateCount * patternCount + j * stateCount;


                    double[] frequencies = substitutionModelDelegate.getRootStateFrequencies(); // TODO May have more than one set of frequencies

                    for (int i = 0; i < stateCount; i++) {
                         if (conditionalProbabilitiesInLogSpace) {
                             conditionalProbabilities[i] = Math.log(partials[partialsIndex + i]) + Math.log(frequencies[i]);
                         } else {
                             conditionalProbabilities[i] = partials[partialsIndex + i] * frequencies[i];
                         }
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

//                 getMatrix(nodeNum, probabilities);
                getTransitionProbabilityMatrix(tree, node, probabilities, true);

                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    int category = rateCategory == null ? 0 : rateCategory[j];
                    int matrixIndex = category * stateCount * stateCount;
                    int partialIndex = category * stateCount * patternCount;

                    for (int i = 0; i < stateCount; i++) {
                        if (conditionalProbabilitiesInLogSpace) {
                            conditionalProbabilities[i] = Math.log(partialLikelihood[partialIndex + childIndex + i])
                                    + Math.log(probabilities[matrixIndex + parentIndex + i]);
                        } else {
                            conditionalProbabilities[i] = partialLikelihood[partialIndex + childIndex + i]
                                    * probabilities[matrixIndex + parentIndex + i];
                        }

                    }

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
            getTransitionProbabilityMatrix(tree, node, probabilities, true);

            if (useAmbiguities()) {

//                 getMatrix(nodeNum, probabilities);
                double[] partials = tipPartials[nodeNum];

                for (int j = 0; j < patternCount; j++) {
                    final int parentIndex = parentState[j] * stateCount;
                    int category = rateCategory == null ? 0 : rateCategory[j];
                    int matrixIndex = category * stateCount * stateCount;

                    int probabilityIndex = parentIndex + matrixIndex;
                    for (int k = 0; k < stateCount; k++) {
                        if(conditionalProbabilitiesInLogSpace){
                            conditionalProbabilities[k] = Math.log(probabilities[probabilityIndex + k])+ Math.log(partials[j * stateCount + k]);
                        }else{
                            conditionalProbabilities[k] = probabilities[probabilityIndex + k]* partials[j * stateCount + k];
                        }
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

//                         getMatrix(nodeNum, probabilities);
                        System.arraycopy(probabilities, parentIndex + matrixIndex, conditionalProbabilities, 0, stateCount);

                        if (useAmbiguities && !dataType.isUnknownState(thisState)) { // Not completely unknown
                            boolean[] stateSet = dataType.getStateSet(thisState);

                            for (int k = 0; k < stateCount; k++) {
                                if (!stateSet[k]) {
                                    conditionalProbabilities[k] = 0.0;
                                }
                            }
                        }

                        if (conditionalProbabilitiesInLogSpace) {
                            for (int k = 0; k < stateCount; k++) {
                                conditionalProbabilities[k] = Math.log(conditionalProbabilities[k]);
                            }
                        }
                        reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
                    }

                    if (!returnMarginalLogLikelihood) {
                        final int parentIndex = parentState[j] * stateCount;
//                         getMatrix(nodeNum, probabilities);
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
    protected List<double[]> probabilitiesAlongBranch;
    protected List<double[]> probabilitiesConvolved;
    protected List<Double> combinedWeights = new ArrayList<Double>();
    protected List<Integer> combinedMatrixOrder = new ArrayList<Integer>();
    protected List<Double> combinedRates = new ArrayList<Double>();
    private double[] partials;

    protected int[] rateCategory = null;
    private final boolean conditionalProbabilitiesInLogSpace;
//    private double[] rootPartials;
//    private int[][] cumulativeScaleBuffers;
//    private int scaleBufferIndex;

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

     private final CodeFormatter formatter;
}
