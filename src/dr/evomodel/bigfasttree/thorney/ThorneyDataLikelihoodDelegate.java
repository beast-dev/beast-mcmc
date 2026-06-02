/*
 * ThorneyDataLikelihoodDelegate.java
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

package dr.evomodel.bigfasttree.thorney;

import java.util.Arrays;
import java.util.List;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.TreeTraversal.TraversalType;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;


//mutationbranch map is akin to patterns list
// branch length likelihood delegate is like a substitution model.
// 

public class ThorneyDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate {
    final static String NAME="ThorneyDataLikelihoodDelegate";



    public ThorneyDataLikelihoodDelegate(Tree tree, MutationBranchMap mutationMap, BranchLengthLikelihoodDelegate branchLengthLikelihoodDelegate) {
        super(NAME);
        this.tree = tree;

        this.mutationMap = mutationMap;

        this.branchLengthLikelihoodDelegate = branchLengthLikelihoodDelegate;
        totalCalculationCount=0;

        updateNode = new boolean[tree.getNodeCount()];
        Arrays.fill(updateNode, true);

        storedUpdateNode  = new boolean[tree.getNodeCount()];
        
        branchLengths = new double[tree.getNodeCount()];
        storedBranchLengths = new double[tree.getNodeCount()];

        branchLogL = new double[tree.getNodeCount()];
        storedBranchLogL = new double[tree.getNodeCount()];
    }

    @Override
    public TraversalType getOptimalTraversalType() {
        // I don't think this matters at this point
        return TreeTraversal.TraversalType.POST_ORDER;
    }

    @Override
    public long getTotalCalculationCount() {
        return totalCalculationCount;
    }

    @Override
    public String getReport() {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'getReport'");
    }

    @Override
    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        //check the substitution model
    }

    @Override
    public void storeState() {
        // TODO Auto-generated method stub

        System.arraycopy(branchLogL, 0, storedBranchLogL, 0, branchLogL.length);
        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, storedBranchLengths.length);    }

    @Override
    public void restoreState() {
        double[] tmp = storedBranchLogL;
        storedBranchLogL=branchLogL;
        branchLogL=tmp;
    
        double[] tmp1 = storedBranchLengths;
        storedBranchLengths = branchLengths;
        branchLengths = tmp1;
    }

    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations,
            int rootNodeNumber) throws LikelihoodException {
                // summing this way accounts for numerical stability but feels like over kill
                for (BranchOperation branchOperation : branchOperations) {
                    NodeRef node = tree.getNode(branchOperation.getBranchNumber());
                    updateNode(node);
                    branchLengths[node.getNumber()] = branchOperation.getBranchLength();
                }
                NodeRef root = tree.getNode(rootNodeNumber);
                int rootChild1 = tree.getChild(root, 0).getNumber();
                int rootChild2 = tree.getChild(root, 1).getNumber();
                updateNode[rootChild1] = updateNode[rootChild1] || updateNode[rootChild2] ;

                totalCalculationCount+=1;

                return calculateLogLikelihood(root, rootNodeNumber, rootChild1,rootChild2);


    }

    public double calculateLogLikelihood(NodeRef node, int root, int rootChild1, int rootChild2){
        int nodeIndex = node.getNumber();

        if (updateNode[nodeIndex]) {
            double logL;
            if(nodeIndex==root){
                logL=0;
            }else{
                MutationList mutations = mutationMap.getMutations(node); //TODO make branchMutations
                double branchLength = branchLengths[nodeIndex];
//                    if (nodeIndex == rootChild1) {
//                        // sum the branches on both sides of the root
//                        NodeRef node2 = tree.getNode(rootChild2);
//                        time += tree.getBranchLength(node2);
//                        mutations += branchLengthProvider.getBranchLength(tree, node2);
//                    }
//                gamma.setScale(1.0);
//                branchLogL[i] = gamma.logPdf(x);
// get mutations/ get branch length / have substitution model handy and return the likelihood of seeing all this mutations 
                logL = branchLengthLikelihoodDelegate.getLogLikelihood(mutations, branchLength);  // Substitution model

            }
            for (int i = 0; i < tree.getChildCount(node); i++) {
                logL += this.calculateLogLikelihood(tree.getChild(node, i),root,rootChild1,rootChild2);
            }
            branchLogL[nodeIndex] = logL;
            updateNode[nodeIndex] = false;
        }
        return branchLogL[nodeIndex];
    }
    /**
     * Set update flag for node and remove it's old contribution to the likelihood.
     * Also handle the root and children so that the 1 branch between children is marked as updated.
     * @param node
     */
    protected void updateNode(NodeRef node) {

        updateNode[node.getNumber()] = true;
        NodeRef parent = tree.getParent(node);
        if (parent != null && !updateNode[parent.getNumber()]) {
            updateNode(parent);
        }
    }


    @Override
    public int getTraitCount() {
        // TODO Auto-generated method stub
        return 1;

    }

    @Override
    public int getTraitDim() {
        // TODO Auto-generated method stub 
        // pattern count in 
        return tree.getNodeCount(); // multiplied substitution model;
        // throw new UnsupportedOperationException("Unimplemented method 'getTraitDim'");
    }

    @Override
    public RateRescalingScheme getRateRescalingScheme() {
        // TODO Auto-generated method stub
        return RateRescalingScheme.NONE;
    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
                // Do nothing like beagle
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        // TODO Auto-generated method stub
        throw new RuntimeException("Not yet implemented");    }

    @Override
    public void setComputePostOrderStatisticsOnly(boolean computePostOrderStatistics) {
                // Do nothing like beagle

    }

    @Override
    public boolean providesPostOrderStatisticsOnly() { return false;    }

    public int getPartitionCat(){
        // not meaningful right now, need to update
        return 0;
    };

    public double[] getSiteLogLikelihoods(){
        throw new RuntimeException("getSiteLogLikelihoods() not implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) { // substitution model here
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'handleModelChangedEvent'");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, ChangeType type) {
        // TODO Auto-generated method stub
        throw new UnsupportedOperationException("Unimplemented method 'handleVariableChangedEvent'");
    }

    @Override
    protected void acceptState() {
        // nada
    }

    public PreOrderSettings getPreOrderSettings() {
        return null;
    }

    @Override
    public boolean getPreferGPU() {
        return true;
    }

    public boolean getUseAmbiguities() {
        return true;
    }

    public PartialsRescalingScheme getRescalingScheme() {
        return null;
    }

    public boolean getDelayRescalingUntilUnderflow() {
        return true;
    }
    // get the mutation map
    protected MutationBranchMap getMutationMap() {
        return mutationMap;
    }

    // get the branch likelihood delegate
    protected BranchLengthLikelihoodDelegate getBranchLengthLikelihoodDelegate() {
        return branchLengthLikelihoodDelegate;
    }

    private long totalCalculationCount;

    private final MutationBranchMap mutationMap;

    private final Tree tree;
    private BranchLengthLikelihoodDelegate branchLengthLikelihoodDelegate;

        /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;
    protected boolean[] storedUpdateNode;

    private double[] branchLengths;
    private double[] storedBranchLengths;

    private double[] branchLogL;
    private double[] storedBranchLogL;
}
