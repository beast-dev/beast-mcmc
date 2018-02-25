/*
 * DataSimulationDelegate.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import beagle.InstanceDetails;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.tree.*;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.*;
import dr.inference.model.Model;

import java.util.List;

import static dr.evolution.tree.TreeTrait.DA.factory;

/**
 * AbstractDiscreteTraitDelegate - interface for a plugin delegate for data simulation on a tree.
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class AbstractDiscreteTraitDelegate extends ProcessSimulationDelegate.AbstractDelegate
        implements TreeTrait.TraitInfo<double[]> {

    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final Beagle beagle;
    private final EvolutionaryProcessDelegate evolutionaryProcessDelegate;
    private final SiteRateModel siteRateModel;
    private final int preOrderPartialOffset;

    public AbstractDiscreteTraitDelegate(String name,
                                  Tree tree,
                                  BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree);
        this.likelihoodDelegate = likelihoodDelegate;
        this.beagle = likelihoodDelegate.getBeagleInstance();
        assert(this.likelihoodDelegate.isUsePreOrder()); /// TODO: reinitialize beagle instance if usePreOrder = false
        evolutionaryProcessDelegate = this.likelihoodDelegate.getEvolutionaryProcessDelegate();
        siteRateModel = this.likelihoodDelegate.getSiteRateModel();

        nodeCount = tree.getNodeCount();
        tipCount = tree.getExternalNodeCount();
        internalNodeCount = nodeCount - tipCount;

        patternCount = likelihoodDelegate.getPatternList().getPatternCount();
        stateCount = likelihoodDelegate.getPatternList().getDataType().getStateCount();
        categoryCount = siteRateModel.getCategoryCount();

        // one partials buffer for each tip and two for each internal node (for store restore)
        partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);
        // mirror the preOrderpartialBufferHelper in BeagleDataLikelihoodDelegate
        // every node gets two partisl (store and restore) except for the root node
        preOrderBufferHelper = new BufferIndexHelper(nodeCount, 1);
        // put preOrder partials in the same order as defined in BufferIndexHelper right after postOrder partials
        preOrderPartialOffset = partialBufferHelper.getBufferCount();

        patternList = likelihoodDelegate.getPatternList();
    }

    @Override
    public void simulate(final int[] operations, final int operationCount,
                         final int rootNodeNumber) {
        //This function updates preOrder Partials for all nodes
        if(DEBUG){
            System.err.println("Setting Root preOrder partial.");
        }

        this.simulateRoot(rootNodeNumber); //set up pre-order partials at root node first

        if(DEBUG){
            System.err.println("Now update preOrder partials at all other nodes");
        }
        int[] beagleoperations = new int[operationCount * Beagle.OPERATION_TUPLE_SIZE];
        int k = 0; int j = 0;
        for(int i = 0; i < operationCount; ++i){
            beagleoperations[k++] = getPreOrderPartialIndex(operations[j++]);
            beagleoperations[k++] = Beagle.NONE;
            beagleoperations[k++] = Beagle.NONE;
            beagleoperations[k++] = getPreOrderPartialIndex(operations[j++]);
            beagleoperations[k++] = evolutionaryProcessDelegate.getMatrixIndex(operations[j++]);
            beagleoperations[k++] = partialBufferHelper.getOffsetIndex(operations[j++]);
            beagleoperations[k++] = evolutionaryProcessDelegate.getMatrixIndex(operations[j++]);
        }
        beagle.updatePrePartials(beagleoperations, operationCount, Beagle.NONE);  // Update all nodes with no rescaling

        //super.simulate(operations, operationCount, rootNodeNumber); // TODO Should override this to compute pre-order partials
        //TODO:error control
    }

    @Override
    public void setupStatistics() {
        throw new RuntimeException("Not used (?) with BEAGLE");
    }

    @Override
    protected void simulateRoot(int rootNumber) {
        //This function sets preOrderPartials at Root for now.
        if (DEBUG) {
            System.err.println("Simulate root node " + rootNumber);
        }
        int[] rootPreIndices = {getPreOrderPartialIndex(rootNumber)};
        int[] rootFreqIndices = {0}; /// as in BeagleDataLikelihoodDelegate.calculateLikelihood()
        InstanceDetails instanceDetails = beagle.getDetails();
        beagle.setRootPrePartials(rootPreIndices, rootFreqIndices, 1);
        //TODO: find the right error message for control
    }

    @Override
    protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
        throw new RuntimeException("Not used with BEAGLE");
    }

    @Override
    protected void constructTraits(Helper treeTraitHelper) {
        treeTraitHelper.addTrait(factory(this));
    }

    public static String getName(String name) {
        return "derivative." + name;
    }

    @Override
    public String getTraitName() {
        return getName(name);
    }

    @Override
    public TreeTrait.Intent getTraitIntent() {
        return TreeTrait.Intent.NODE;
    }

    @Override
    public Class getTraitClass() {
        return double[].class;
    }

    @Override
    public double[] getTrait(Tree tree, NodeRef node) {
        assert (tree == this.tree);
        assert (node == null); // Implies: get trait for all nodes at same time

        //update all preOrder partials first
        simulationProcess.cacheSimulatedTraits(node);

        final double[] postOrderPartial = new double[stateCount * patternCount * categoryCount];
        final double[] preOrderPartial = new double[stateCount * patternCount * categoryCount];
        final double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();
        final double[] rootPostOrderPartials = new double[stateCount * patternCount * categoryCount];
        final double[] gradient = new double[patternList.getStateCount()];

        //create a matrix for fetching the infinitesimal matrix Q
        double [] Q = new  double[stateCount * stateCount];
        double[] clikelihood = new double[categoryCount * patternCount];  // likelihood for each category doesn't come in free.

        beagle.getPartials(partialBufferHelper.getOffsetIndex(tree.getRoot().getNumber()), Beagle.NONE, rootPostOrderPartials);

        double [] grand_denominator = new double[patternCount];
        double [] grand_numerator = new double[patternCount];

        int l, j, k, s, t, m, u;
        for(m = 0; m < patternCount; ++m){
            grand_denominator[m] = 0;
            grand_numerator[m] = 0;
        }

        beagle.getPartials(partialBufferHelper.getOffsetIndex(node.getNumber()), Beagle.NONE, postOrderPartial);
        beagle.getPartials(getPreOrderPartialIndex(node.getNumber()), Beagle.NONE, preOrderPartial);

        evolutionaryProcessDelegate.getSubstitutionModel(node.getNumber()).getInfinitesimalMatrix(Q);  //store the Q matrix

        double[] tmpNumerator = new double[patternCount * categoryCount];
        l = 0; j = 0;
        for(s = 0; s < categoryCount; s++){
            for(m = 0; m < patternCount; m++){
                double clikelihood_tmp = 0;
                for( k = 0; k < stateCount; k++){
                    clikelihood_tmp += frequencies[k] * rootPostOrderPartials[l++];
                }
                clikelihood[j++] = clikelihood_tmp;
            }
        }

        //now calculate weights
        t = 0;
        double denominator = 0;
        double numerator = 0;
        double tmp = 0;
        double[] weights = siteRateModel.getCategoryProportions();
        for(s = 0; s < categoryCount; s++){
            double rs = siteRateModel.getRateForCategory(s);
            double ws = weights[s];
            u = 0;
            for(m = 0; m < patternCount; m++){
                l = 0;
                numerator = 0;
                denominator = 0;
                for(k = 0; k < stateCount; k++){
                    tmp = 0;
                    for(j = 0; j < stateCount; j++){
                        tmp += Q[l++] * postOrderPartial[u + j];
                    }
                    numerator += tmp * preOrderPartial[u + k];
                    denominator += postOrderPartial[u + k] * preOrderPartial[u + k];
                }
                u += stateCount;
                tmpNumerator[t] = ws * rs * numerator / denominator * clikelihood[t];
                grand_numerator[m] += tmpNumerator[t];
                grand_denominator[m] += ws * clikelihood[t];
                t++;
            }
        }
        for( m = 0; m < patternList.getStateCount(); m++){
            int sitePatternIndex = ((SitePatterns) patternList).getPatternIndex(m);
            gradient[m] = grand_numerator[sitePatternIndex] / grand_denominator[sitePatternIndex];
        }
        // TODO See TipGradientViaFullConditionalDelegate.getTrait() as an example of using post- and pre-order partials together
        return gradient;
    }

    @Override
    public boolean isTraitLoggable() {
        return false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        // TODO
    }

    @Override
    public void modelRestored(Model model) {

    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        int k = 0;
        for(int i = 0; i < nodeOperations.size(); ++i){
            NodeOperation tmpNodeOperation = nodeOperations.get(i);
            //nodeNumber = ParentNodeNumber, leftChild = nodeNumber, rightChild = siblingNodeNumber
            operations[k++] = tmpNodeOperation.getLeftChild();
            operations[k++] = tmpNodeOperation.getNodeNumber();
            operations[k++] = tmpNodeOperation.getLeftChild();
            operations[k++] = tmpNodeOperation.getRightChild();
            operations[k++] = tmpNodeOperation.getRightChild();
        }
        flipBufferIndices(nodeOperations);
        return nodeOperations.size();
    }

    //This function defines the mapping from nodeNum to its preorder partial index
    final private int getPreOrderPartialIndex(final int i){
        NodeRef nodeRef = tree.getNode(i);
        if(tree.isRoot(nodeRef)){
            return preOrderPartialOffset + 0;
        }else{
            return preOrderPartialOffset + preOrderBufferHelper.getOffsetIndex(i + 1);
        }
    }

    final private void flipBufferIndices(List<NodeOperation> nodeOperations) {
        // Flip all the buffers to be written to first...
        for (int nodeNum = tipCount; nodeNum < nodeCount; ++nodeNum){
            this.partialBufferHelper.flipOffset(nodeNum);
        }

        // Now root = 0 (stored as a tip in postOrder), all other nodeNum += 1
        for (int nodeNum = 1; nodeNum < nodeCount; ++nodeNum){
            this.preOrderBufferHelper.flipOffset(nodeNum);
        }
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private final int nodeCount;
    private final int tipCount;
    private final int internalNodeCount;


    private final int patternCount;
    private final int stateCount;
    private final int categoryCount;

    private final BufferIndexHelper partialBufferHelper;
    private final BufferIndexHelper preOrderBufferHelper;
    private final PatternList patternList;

    private static final boolean DEBUG = true;
}

