/*
 * AbstractTreeLikelihood.java
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

package dr.evomodel.bigFastTree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;
import org.apache.commons.math.special.Gamma;
import org.apache.commons.math.util.FastMath;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * ApproximatePoissonTreeLikelihood - a tree likelihood which uses an ML tree as expected number
 * of substitutions and assumes that these are drawn from a gamma (with mean == variance).
 *
 * This is similar to work by Jeff Thorne and colleagues:
 *
 *
 * And more recently Didelot et al (2018) BioRxiv
 *
 * @author Andrew Rambaut
 * @author JT McCrone
 */

public class ApproximatePoissonTreeLikelihood extends AbstractModelLikelihood implements Reportable {
    public ApproximatePoissonTreeLikelihood(String name, Tree dataTree, int sequenceLength, TreeModel treeModel, BranchRateModel branchRateModel) {

        super(name);

        this.sequenceLength = sequenceLength;

        this.treeModel = treeModel;
        addModel(treeModel);
        cachedRoot = treeModel.getRoot().getNumber();
        cachedRootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        cachedRootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();

        this.branchRateModel = branchRateModel;
        addModel(branchRateModel);


        HashMap<BitSet, NodeRef> dataTreeMap = getBitSetNodeMap(dataTree,dataTree);

        HashMap<BitSet, NodeRef> treeModelMap = getBitSetNodeMap(dataTree,treeModel);
        // reverse map to be <NodeRef, BitSet>
        HashMap <NodeRef, BitSet> treeModelNodeMap = new HashMap<>();
        for (Map.Entry<BitSet, NodeRef> entry: treeModelMap.entrySet()){
            treeModelNodeMap.put(entry.getValue(), entry.getKey());
        }
        // An array where the entry points to the node in the data tree that maps to the node in the time tree
        // Nodes that result from resolving a polytomy in the data tree will point to the polytomy node.
        this.nodeInDataTree = new int[this.treeModel.getNodeCount()];


        setUpNodeMap(treeModelNodeMap,dataTreeMap,treeModel.getRoot(),null);


        updateNode = new boolean[treeModel.getNodeCount()];
        branchLengths = new double[dataTree.getNodeCount()];

        for (int i = 0; i < dataTree.getNodeCount(); i++) {
                double x = dataTree.getBranchLength(dataTree.getNode(i)) * sequenceLength;
                branchLengths[i] = Math.round(x);
        }
        for (int i = 0; i <this.treeModel.getNodeCount() ; i++) {
            updateNode[i] = true;
        }

//        distanceMatrix = new double[treeModel.getExternalNodeCount()][];
//        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
//            distanceMatrix[i] = new double[treeModel.getExternalNodeCount()];
//            for (int j = i + 1; j < treeModel.getExternalNodeCount(); j++) {
//                distanceMatrix[i][j] = distanceMatrix[j][i] = TreeUtils.getPathLength(dataTree, dataTree.getNode(i), dataTree.getNode(j));
//            }
//        }

        branchLogL = new double[treeModel.getNodeCount()];
        storedBranchLogL = new double[treeModel.getNodeCount()];
        likelihoodKnown = false;
    }

    /**
     * A private recursive method that sets the map from the treemodel to the data tree. If the data tree is resolved it
     * will be a 1 to 1 map. If there are polytomies in the data tree all inserted nodes will map to the polytomy node.
      * @param treeModelNodeMap
     * @param dataTreeMap
     * @param node
     * @param parentsClade
     */
    private void setUpNodeMap(HashMap <NodeRef, BitSet> treeModelNodeMap, HashMap<BitSet, NodeRef> dataTreeMap,NodeRef node, BitSet parentsClade){
        BitSet clade = treeModelNodeMap.get(node);
        int j = node.getNumber();
        if(!dataTreeMap.containsKey(clade)){
           clade=parentsClade;
        }
        nodeInDataTree[j] = dataTreeMap.get(clade).getNumber();
        for (int i = 0; i <treeModel.getChildCount(node) ; i++) {
            NodeRef child = treeModel.getChild(node,i);
            setUpNodeMap(treeModelNodeMap,dataTreeMap,child,clade);
        }

    }

    private double getNumberOfMutations(int i) {
        int dataNode = nodeInDataTree[i];
        int dataParentNode = nodeInDataTree[treeModel.getParent(treeModel.getNode(i)).getNumber()];
        if(dataNode==dataParentNode){
            return 0d;
        }else{
            return branchLengths[dataNode];
        }

    }
    /**
     * Gets a HashMap of clade bitsets to nodes in tree. This is useful for comparing the topology of trees
     * @param referenceTree  the tree that will be used to define taxa and tip numbers
     * @param tree the tree for which clades are being defined
     * @return A HashMap with a BitSet of descendent taxa as the key and a node as value
     */
    private HashMap<BitSet, NodeRef> getBitSetNodeMap(Tree referenceTree,Tree  tree) {
        HashMap<BitSet, NodeRef> map = new HashMap<>();
        addBits(referenceTree,tree,tree.getRoot(),map);
        return map;
    }

    /**
     *  A private recursive function used by getBitSetNodeMap
     *  This is modeled after the addClades in CladeSet and getClades in compatibility statistic
     * @param referenceTree  the tree that will be used to define taxa and tip numbers
     * @param tree the tree for which clades are being defined
     * @param node current node
     * @param map map that is being appended to
     */
    private BitSet addBits(Tree referenceTree, Tree tree, NodeRef node, HashMap map) {
        BitSet bits = new BitSet();
        if (tree.isExternal(node)) {
            String taxonId = tree.getNodeTaxon(node).getId();
            bits.set(referenceTree.getTaxonIndex(taxonId));

        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef node1 = tree.getChild(node, i);
                bits.or(addBits(referenceTree,tree, node1, map));
            }
        }

        map.put(bits, node);
        return bits;
    }

    /**
     * Set update flag for node and remove it's old contribution to the likelihood.
     * Also handle the root and children so that the 1 branch between children is marked as updated.
     * @param node
     */
    protected void updateNode(NodeRef node) {

        updateNode[node.getNumber()] = true;
        NodeRef parent = treeModel.getParent(node);
        if (parent != null && !updateNode[parent.getNumber()]) {
            updateNode(parent);
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its direct children
     */
    protected void updateNodeAndChildren(NodeRef node) {

        updateNode(node);

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNode(child);
        }

//        likelihoodKnown = false;
    }



    /**
     * Set update flag for a node and all its descendents
     */
    protected void updateNodeAndDescendents(NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNodeAndDescendents(child);
        }

        likelihoodKnown = false;
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            updateNode[i] = true;
        }
        likelihoodKnown = false;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    /**
     * Handles model changed events from the submodels.
     */
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        fireModelChanged();

        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    NodeRef node= ((TreeChangedEvent) object).getNode();
                    updateNodeAndChildren(node);
                    if(treeModel.getRoot()== node){
                        // If the root is changed then the old root's children must be updated so that they are
                        // no longer counted as 1 branch in the likelihood
                        // This currently updates more often than needed.
                        NodeRef oldRootNode = treeModel.getNode(cachedRoot);
                        if(oldRootNode!=node){
                            updateNodeAndChildren(oldRootNode);

                        }
                    }

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }

        } else if (model == branchRateModel) {
            if (index == -1) {
                updateAllNodes();
            } else {
                updateNode(treeModel.getNode(index));
            }

        } else {

            throw new RuntimeException("Unknown componentChangedEvent");
        }
    }

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

        storedCachedRoot = cachedRoot;
        storedCachedRootChild1 = cachedRootChild1;
        storedCachedRootChild2 = cachedRootChild2;


        System.arraycopy(branchLogL, 0, storedBranchLogL, 0, branchLogL.length);
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        cachedRoot = storedCachedRoot;
        cachedRootChild1 = storedCachedRootChild1;
        cachedRootChild2 = storedCachedRootChild2;
        double[] tmp = storedBranchLogL;
        storedBranchLogL = branchLogL;
        branchLogL = tmp;
    }

    protected void acceptState() {
        cachedRoot = treeModel.getRoot().getNumber();
        cachedRootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        cachedRootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    private double calculateLogLikelihood(NodeRef node,int root, int rootChild1, int rootChild2) {
        int nodeIndex = node.getNumber();

        if (updateNode[nodeIndex]) {
            double logL;
            if(nodeIndex==root || nodeIndex==rootChild2){
                logL=0;
            }else{
                double expected = treeModel.getBranchLength(node) * branchRateModel.getBranchRate(treeModel, node);

                double x = getNumberOfMutations(nodeIndex);

                if (nodeIndex == rootChild1) {
                    // sum the branches on both sides of the root
                    NodeRef node2 = treeModel.getNode(rootChild2);
                    expected += treeModel.getBranchLength(node2) * branchRateModel.getBranchRate(treeModel, node2);
                    x += getNumberOfMutations(rootChild2);
                }
                double mean = expected * sequenceLength;

//                gamma.setScale(1.0);
//                branchLogL[i] = gamma.logPdf(x);
                logL = SaddlePointExpansion.logPoissonProbability(mean, (int) x); //SaddlePointExpansion.logBinomialProbability((int)x, sequenceLength, expected, 1.0D - expected);
            }

            for (int i = 0; i < treeModel.getChildCount(node); i++) {
                logL += calculateLogLikelihood(treeModel.getChild(node, i),root,rootChild1,rootChild2);
            }
            branchLogL[nodeIndex] = logL;
            updateNode[nodeIndex] = false;
        }
        return branchLogL[nodeIndex];
    }

    private double calculateLogLikelihoodLinearTraversal() {

//        makeDirty();
//Could make this faster by only adding the changed values
        int root = treeModel.getRoot().getNumber();
        int rootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        int rootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();


        double logL = 0.0;
        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            if ( i != root && i != rootChild2) {
                NodeRef node = treeModel.getNode(i);
                // skip the root and the second child of the root (this is added to the first child)

                double expected = treeModel.getBranchLength(node) * branchRateModel.getBranchRate(treeModel, node);
                double x = getNumberOfMutations(i);

                if (i == rootChild1) {
                    // sum the branches on both sides of the root
                    NodeRef node2 = treeModel.getNode(rootChild2);
                    expected += treeModel.getBranchLength(node2) * branchRateModel.getBranchRate(treeModel, node2);
                    x += getNumberOfMutations(rootChild2);
                }
                double mean = expected * sequenceLength;

//                gamma.setScale(1.0);
//                branchLogL[i] = gamma.logPdf(x);

                logL += SaddlePointExpansion.logPoissonProbability(mean,(int)x); //SaddlePointExpansion.logBinomialProbability((int)x, sequenceLength, expected, 1.0D - expected);
            }
        }

        return logL;
    }

    private double calculateLogLikelihood() {


//Could make this faster by only adding the changed values
        int root = treeModel.getRoot().getNumber();
        int rootChild1 = treeModel.getChild(treeModel.getRoot(), 0).getNumber();
        int rootChild2 = treeModel.getChild(treeModel.getRoot(), 1).getNumber();
        updateNode[rootChild1] = updateNode[rootChild1] || updateNode[rootChild2] ;
        assert updateNode[root];


        return calculateLogLikelihood(treeModel.getRoot(), root, rootChild1, rootChild2);

    }


    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            if(testing){
                makeDirty();
                assert logLikelihood==calculateLogLikelihood();
                double linearLL = calculateLogLikelihoodLinearTraversal();
                assert linearLL !=logLikelihood;
                assert Math.abs(linearLL-logLikelihood)<1E-10;
            }
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        likelihoodKnown = false;
        updateAllNodes();
    }

    public String getReport() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }


    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the tree
     */
    private final TreeModel treeModel;

    private final BranchRateModel branchRateModel;

    private final int sequenceLength;

    private final double[] branchLengths;
    //private final double[][] distanceMatrix;

    /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;
    private boolean storedLikelihoodKnown = false;

    private double[] branchLogL;
    private double[] storedBranchLogL;
    private final int[] nodeInDataTree;
    private int cachedRoot;
    private int storedCachedRoot;
    private int cachedRootChild1;
    private int cachedRootChild2;
    private int storedCachedRootChild1;
    private int storedCachedRootChild2;
    private final boolean testing =false;
}

// Grabbed some stuff from Commons Maths as it is not public
// This code is under the Apache License 2.0

final class SaddlePointExpansion {
    private static final double HALF_LOG_2_PI = 0.5D * FastMath.log(6.283185307179586D);
    private static final double[] EXACT_STIRLING_ERRORS = new double[]{0.0D, 0.15342640972002736D, 0.08106146679532726D, 0.05481412105191765D, 0.0413406959554093D, 0.03316287351993629D, 0.02767792568499834D, 0.023746163656297496D, 0.020790672103765093D, 0.018488450532673187D, 0.016644691189821193D, 0.015134973221917378D, 0.013876128823070748D, 0.012810465242920227D, 0.01189670994589177D, 0.011104559758206917D, 0.010411265261972096D, 0.009799416126158804D, 0.009255462182712733D, 0.008768700134139386D, 0.00833056343336287D, 0.00793411456431402D, 0.007573675487951841D, 0.007244554301320383D, 0.00694284010720953D, 0.006665247032707682D, 0.006408994188004207D, 0.006171712263039458D, 0.0059513701127588475D, 0.0057462165130101155D, 0.005554733551962801D};

    private SaddlePointExpansion() {
    }

    static double getStirlingError(double z) {
        double ret;
        double z2;
        if (z < 15.0D) {
            z2 = 2.0D * z;
            if (FastMath.floor(z2) == z2) {
                ret = EXACT_STIRLING_ERRORS[(int)z2];
            } else {
                ret = Gamma.logGamma(z + 1.0D) - (z + 0.5D) * FastMath.log(z) + z - HALF_LOG_2_PI;
            }
        } else {
            z2 = z * z;
            ret = (0.08333333333333333D - (0.002777777777777778D - (7.936507936507937E-4D - (5.952380952380953E-4D - 8.417508417508417E-4D / z2) / z2) / z2) / z2) / z;
        }

        return ret;
    }

    static double getDeviancePart(double x, double mu) {
        double ret;
        if (FastMath.abs(x - mu) < 0.1D * (x + mu)) {
            double d = x - mu;
            double v = d / (x + mu);
            double s1 = v * d;
            double s = 0.0D / 0.0;
            double ej = 2.0D * x * v;
            v *= v;

            for(int j = 1; s1 != s; ++j) {
                s = s1;
                ej *= v;
                s1 += ej / (double)(j * 2 + 1);
            }

            ret = s1;
        } else {
            ret = x * FastMath.log(x / mu) + mu - x;
        }

        return ret;
    }

    static double logBinomialProbability(int x, int n, double p, double q) {
        double ret;
        if (x == 0) {
            if (p < 0.1D) {
                ret = -getDeviancePart((double)n, (double)n * q) - (double)n * p;
            } else {
                ret = (double)n * FastMath.log(q);
            }
        } else if (x == n) {
            if (q < 0.1D) {
                ret = -getDeviancePart((double)n, (double)n * p) - (double)n * q;
            } else {
                ret = (double)n * FastMath.log(p);
            }
        } else {
            ret = getStirlingError((double)n) - getStirlingError((double)x) - getStirlingError((double)(n - x)) - getDeviancePart((double)x, (double)n * p) - getDeviancePart((double)(n - x), (double)n * q);
            double f = 6.283185307179586D * (double)x * (double)(n - x) / (double)n;
            ret += -0.5D * FastMath.log(f);
        }

        return ret;
    }

    static public double logPoissonProbability(double mean,int x) {
        double ret;
        if (x >= 0 && x != 2147483647) {
            if (x == 0) {
                ret = FastMath.exp(-mean);
            } else {
                ret = FastMath.exp(-getStirlingError((double)x) - getDeviancePart((double)x, mean)) / FastMath.sqrt(6.283185307179586D * (double)x);
            }
        } else {
            ret = 0.0D;
        }

        return Math.log(ret);
    }
}
