/*
 * AncestralStateTreeLikelihood.java
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

package dr.oldevomodel.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.math.MathUtils;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
@Deprecated // Switching to BEAGLE
public class AncestralStateTreeLikelihood extends TreeLikelihood implements TreeTraitProvider {
    public static final String STATES_KEY = "states";

//    private boolean useExtraReconstructedStates = false;


    /**
     * Constructor.
     * Now also takes a DataType so that ancestral states are printed using data codes
     *
     * @param patternList     -
     * @param treeModel       -
     * @param siteModel       -
     * @param branchRateModel -
     * @param useAmbiguities  -
     * @param storePartials   -
     * @param dataType        - need to provide the data-type, so that corrent data characters can be returned
     * @param tag             - string label for reconstruction characters in tree log
     * @param forceRescaling  -
     * @param useMAP          - perform maximum aposteriori reconstruction
     * @param returnML        - report integrate likelihood of tip data
     */
    public AncestralStateTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                        SiteModel siteModel, BranchRateModel branchRateModel,
                                        boolean useAmbiguities, boolean storePartials,
                                        final DataType dataType,
                                        final String tag,
                                        boolean forceRescaling,
                                        boolean useMAP,
                                        boolean returnML) {
        super(patternList, treeModel, siteModel, branchRateModel, null, useAmbiguities, false, storePartials,
                false, forceRescaling);
        this.dataType = dataType;
        this.tag = tag;

        reconstructedStates = new int[treeModel.getNodeCount()][patternCount];
        storedReconstructedStates = new int[treeModel.getNodeCount()][patternCount];

        this.useMAP = useMAP;
        this.returnMarginalLogLikelihood = returnML;
      
        treeTraits.addTrait(STATES_KEY, new TreeTrait.IA() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public int[] getTrait(Tree tree, NodeRef node) {
                return getStatesForNode(tree,node);
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formattedState(getStatesForNode(tree,node), dataType);
            }
        });

        if (useAmbiguities) {
            Logger.getLogger("dr.evomodel.treelikelihood").info("Ancestral reconstruction using ambiguities is currently "+
            "not support without BEAGLE");
            System.exit(-1);
        }

    }

    public AncestralStateTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                        SiteModel siteModel, BranchRateModel branchRateModel,
                                        boolean useAmbiguities, boolean storePartials,
                                        DataType dataType,
                                        String tag,
                                        boolean forceRescaling) {
        this(patternList, treeModel, siteModel, branchRateModel, useAmbiguities,
                storePartials, dataType, tag, forceRescaling, false, true);
    }

    public void storeState() {

        super.storeState();

        for (int i = 0; i < reconstructedStates.length; i++) {
            System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
        }

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

    public DataType getDataType() {
        return dataType;
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


    public void redrawAncestralStates() {
        jointLogLikelihood = 0;
        traverseSample(treeModel, treeModel.getRoot(), null);
        areStatesRedrawn = true;
    }

//    private boolean checkConditioning = true;


    protected void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);
        fireModelChanged(model);

    }

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

    protected TreeTraitProvider.Helper treeTraits = new Helper();

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
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
                sb.append(dataType.getChar(i));
            }
        }
        sb.append("\"");
        return sb.toString();
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

    /**
     * Traverse (pre-order) the tree sampling the internal node states.
     *
     * @param tree        - TreeModel on which to perform sampling
     * @param node        - current node
     * @param parentState - character state of the parent node to 'node'
     */
    public void traverseSample(TreeModel tree, NodeRef node, int[] parentState) {

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        // This function assumes that all partial likelihoods have already been calculated
        // If the node is internal, then sample its state given the state of its parent (pre-order traversal).

        double[] conditionalProbabilities = new double[stateCount];
        int[] state = new int[patternCount];

        if (!tree.isExternal(node)) {

            if (parent == null) {

                double[] rootPartials = getRootPartials();

                // This is the root node
                for (int j = 0; j < patternCount; j++) {

                    System.arraycopy(rootPartials, j * stateCount, conditionalProbabilities, 0, stateCount);
                    double[] frequencies = frequencyModel.getFrequencies();
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

                    //System.out.println("Pr(j) = " + frequencies[state[j]]);
                    jointLogLikelihood += Math.log(frequencies[state[j]]);
                }

            } else {

                // This is an internal node, but not the root
                double[] partialLikelihood = new double[stateCount * patternCount];
                if (categoryCount > 1)
                    throw new RuntimeException("Reconstruction not implemented for multiple categories yet.");

                likelihoodCore.getPartials(nodeNum, partialLikelihood);

//				final double branchRate = branchRateModel.getBranchRate(tree, node);
//
//				            // Get the operational time of the branch
//				final double branchTime = branchRate * ( tree.getNodeHeight(parent) - tree.getNodeHeight(node) );
//
//				for (int i = 0; i < categoryCount; i++) {
//
//				                siteModel.getTransitionProbabilitiesForCategory(i, branchTime, probabilities);
//
//				}
//




                ((AbstractLikelihoodCore) likelihoodCore).getNodeMatrix(nodeNum, 0, probabilities);


                for (int j = 0; j < patternCount; j++) {

                    int parentIndex = parentState[j] * stateCount;
                    int childIndex = j * stateCount;

                    for (int i = 0; i < stateCount; i++) {
                        conditionalProbabilities[i] = partialLikelihood[childIndex + i] * probabilities[parentIndex + i];
                    }

                    state[j] = drawChoice(conditionalProbabilities);
                    reconstructedStates[nodeNum][j] = state[j];

                    double contrib = probabilities[parentIndex + state[j]];
                    //System.out.println("Pr(" + parentState[j] + ", " + state[j] +  ") = " + contrib);
                    jointLogLikelihood += Math.log(contrib);
                }
            }

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            traverseSample(tree, child1, state);

            NodeRef child2 = tree.getChild(node, 1);
            traverseSample(tree, child2, state);
        } else {

            // This is an external leaf

            ((AbstractLikelihoodCore) likelihoodCore).getNodeStates(nodeNum, reconstructedStates[nodeNum]);

            // Check for ambiguity codes and sample them

            for (int j = 0; j < patternCount; j++) {

                final int thisState = reconstructedStates[nodeNum][j];
                final int parentIndex = parentState[j] * stateCount;
                ((AbstractLikelihoodCore) likelihoodCore).getNodeMatrix(nodeNum, 0, probabilities);
                if (dataType.isAmbiguousState(thisState)) {

                    System.arraycopy(probabilities, parentIndex, conditionalProbabilities, 0, stateCount);
                    reconstructedStates[nodeNum][j] = drawChoice(conditionalProbabilities);
                }

                double contrib = probabilities[parentIndex + reconstructedStates[nodeNum][j]];
                //System.out.println("Pr(" + parentState[j] + ", " + reconstructedStates[nodeNum][j] +  ") = " + contrib);
                jointLogLikelihood += Math.log(contrib);
            }
        }
    }

    private DataType dataType;
    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;

    private String tag;
    private boolean areStatesRedrawn = false;
    private boolean storedAreStatesRedrawn = false;

    private boolean useMAP = false;
    private boolean returnMarginalLogLikelihood = true;

    private double jointLogLikelihood;
    private double storedJointLogLikelihood;


}
