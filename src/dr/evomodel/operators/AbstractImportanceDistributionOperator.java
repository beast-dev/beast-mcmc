/*
 * AbstractImportanceDistributionOperator.java
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

/**
 *
 */
package dr.evomodel.operators;

import dr.evolution.tree.Clade;
import dr.evolution.tree.MutableTree.InvalidTreeException;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.AbstractCladeImportanceDistribution;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.operators.*;
import dr.math.MathUtils;

import java.util.*;

/**
 * @author Sebastian Hoehna
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public abstract class AbstractImportanceDistributionOperator extends
        SimpleMCMCOperator implements GeneralOperator {

    private long transitions = 0;

    private OperatorSchedule schedule;

    protected TreeModel tree;

    protected AbstractCladeImportanceDistribution probabilityEstimater;

    private int sampleEvery;

    private int samples;

    private int sampleCount;

    private Queue<NodeRef> internalNodes;

    private Map<Integer, NodeRef> externalNodes;

    private boolean burnin = false;

    /**
     *
     */
    public AbstractImportanceDistributionOperator(TreeModel tree, double weight) {
        super();

        this.tree = tree;

        setWeight(weight);
        this.samples = 10000;
        this.sampleEvery = 10;

        init();
    }

    /**
     *
     */
    public AbstractImportanceDistributionOperator(TreeModel tree,
                                                  double weight, int samples, int sampleEvery) {
        super();

        this.tree = tree;
        setWeight(weight);
        this.samples = samples;
        this.sampleEvery = sampleEvery;

        init();
    }

    private void init() {
        schedule = getOperatorSchedule(tree);
        sampleCount = 0;
        internalNodes = new LinkedList<NodeRef>();
        externalNodes = new HashMap<Integer, NodeRef>();
        fillExternalNodes(tree.getRoot());
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.AbstractImportanceSampler#doOperation()
      */
    public double doOperation() {
        // dummy method
        return 0.0;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.AbstractImportanceSampler#doOperation()
      */
    public double doOperation(Likelihood likelihood) {
        if (!burnin) {
            if (sampleCount < samples * sampleEvery) {
                sampleCount++;
                if (sampleCount % sampleEvery == 0) {
                    probabilityEstimater.addTree(tree);
                }
                setAcceptCount(0);
                setRejectCount(0);
                setTransitions(0);

                return doUnguidedOperation();

            } else {
                return doImportanceDistributionOperation(likelihood);
            }
        } else {

            return doUnguidedOperation();

        }
    }

    protected double doImportanceDistributionOperation(Likelihood likelihood) {
        final NodeRef root = tree.getRoot();
        BitSet all = new BitSet();
        all.set(0, (tree.getNodeCount() + 1) / 2);
        Clade rootClade = new Clade(all, tree.getNodeHeight(root));

        internalNodes.clear();
        fillInternalNodes(root);
        // remove the root
        internalNodes.poll();

        externalNodes.clear();
        fillExternalNodes(root);

        double prob;
        double back = probabilityEstimater.getTreeProbability(tree);

        try {
            tree.beginTreeEdit();

            List<Clade> originalClades = new ArrayList<Clade>();
            extractClades(tree, tree.getRoot(), originalClades, null);

            double[] originalNodeHeights = getAbsoluteNodeHeights(originalClades);
            Arrays.sort(originalNodeHeights);
            back += getChanceForNodeHeights(originalNodeHeights);

            prob = createTree(root, rootClade);

            assignDummyHeights(root);

//			assignCladeHeights(tree.getRoot(), originalClades, null);
//			double[] originalNodeHeights = getAbsoluteNodeHeights(originalClades);
//			Arrays.sort(originalNodeHeights);
//			prob += setMissingNodeHeights(tree.getChild(tree.getRoot(),0));
//			prob += setMissingNodeHeights(tree.getChild(tree.getRoot(),1));

            prob += setNodeHeights(originalNodeHeights);

//			List<Clade> newClades = new ArrayList<Clade>();
//			extractClades(tree, tree.getRoot(), newClades, null);

            tree.endTreeEdit();

            tree.checkTreeIsValid();
        } catch (InvalidTreeException e) {
            throw new RuntimeException(e.getMessage());
        }

        tree.pushTreeChangedEvent(root);

        return back - prob;
    }

    private void assignDummyHeights(NodeRef node) {
        double rootHeight = tree.getNodeHeight(node) * tree.getInternalNodeCount();
        tree.setNodeHeight(node, rootHeight);

        int childcount = tree.getChildCount(node);
        for (int i = 0; i < childcount; i++) {
            NodeRef child = tree.getChild(node, i);
            if (!tree.isExternal(child)) {
                assignDummyHeights(child, rootHeight / 2.0);
            }
        }

    }

    private void assignDummyHeights(NodeRef node, double height) {
        assert (!tree.isExternal(node));
        tree.setNodeHeight(node, height);

        int childcount = tree.getChildCount(node);
        for (int i = 0; i < childcount; i++) {
            NodeRef child = tree.getChild(node, i);
            if (!tree.isExternal(child)) {
                assignDummyHeights(child, height / 2.0);
            }
        }
    }

    private double createTree(NodeRef node, Clade c)
            throws InvalidTreeException {
        double prob = 0.0;
        if (c.getSize() == 2) {
            // this clade only contains two tips
            // the split between them is trivial

            int leftTipIndex = c.getBits().nextSetBit(0);
            int rightTipIndex = c.getBits().nextSetBit(leftTipIndex + 1);
            NodeRef leftTip = externalNodes.get(leftTipIndex);
            NodeRef rightTip = externalNodes.get(rightTipIndex);

            removeChildren(node);
            NodeRef leftParent = tree.getParent(leftTip);
            if (leftParent != null)
                tree.removeChild(leftParent, leftTip);
            NodeRef rightParent = tree.getParent(rightTip);
            if (rightParent != null)
                tree.removeChild(rightParent, rightTip);
            tree.addChild(node, leftTip);
            tree.addChild(node, rightTip);
        } else {
            Clade[] clades = new Clade[2];
            prob = splitClade(c, clades);
            NodeRef leftChild, rightChild;

            if (clades[0].getSize() == 1) {
                int tipIndex = clades[0].getBits().nextSetBit(0);
                leftChild = externalNodes.get(tipIndex);
            } else {
                leftChild = internalNodes.poll();
                // TODO set the node height for the new node
                tree.setNodeHeight(leftChild, tree.getNodeHeight(node) * 0.5);
                prob += createTree(leftChild, clades[0]);
            }

            if (clades[1].getSize() == 1) {
                int tipIndex = clades[1].getBits().nextSetBit(0);
                rightChild = externalNodes.get(tipIndex);
            } else {
                rightChild = internalNodes.poll();
                // TODO set the node height for the new node
                tree.setNodeHeight(rightChild, tree.getNodeHeight(node) * 0.5);
                prob += createTree(rightChild, clades[1]);
            }

            removeChildren(node);
            NodeRef leftParent = tree.getParent(leftChild);
            if (leftParent != null)
                tree.removeChild(leftParent, leftChild);
            NodeRef rightParent = tree.getParent(rightChild);
            if (rightParent != null)
                tree.removeChild(rightParent, rightChild);
            tree.addChild(node, leftChild);
            tree.addChild(node, rightChild);
        }

        return prob;
    }

    /**
     * @param parent
     * @warning assumes strictly bifurcating trees
     */
    private void removeChildren(NodeRef parent) {
        // assumes strictly bifurcating trees

        NodeRef child = tree.getChild(parent, 0);
        if (child != null) {
            tree.removeChild(parent, child);
        }

        child = tree.getChild(parent, 1);
        if (child != null) {
            tree.removeChild(parent, child);
        }
    }

    private double splitClade(Clade c, Clade[] children) {
        return probabilityEstimater.splitClade(c, children);
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param tree   - the tree from which the clades are extracted
     * @param node   - the starting node. All clades below starting at this branch
     *               are added
     * @param clades - the list in which the clades are stored
     * @param bits   - a bit set to which the current bits of the clades are added
     */
    private void extractClades(Tree tree, NodeRef node,
                               List<Clade> clades, BitSet bits) {

        // create a new bit set for this clade
        BitSet bits2 = new BitSet();

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            int index = node.getNumber();
            bits2.set(index);

        } else {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                extractClades(tree, child, clades, bits2);
            }
            // add my bit set to the list
            clades.add(new Clade(bits2, tree.getNodeHeight(node)));
        }

        // add my bit set to the bit set I was given
        // this is needed for adding all children clades together
        if (bits != null) {
            bits.or(bits2);
        }
    }

    /**
     * Creates a list with all clades of the tree
     *
     * @param node   - the starting node. All clades below starting at this branch
     *               are added
     * @param clades - the list in which the clades are stored
     */
    private void assignCladeHeights(NodeRef node,
                                    HashMap<Clade, Double> clades, BitSet bits) {

        // create a new bit set for this clade
        BitSet bits2 = new BitSet();

        // check if the node is external
        if (tree.isExternal(node)) {

            // if so, the only taxon in the clade is I
            int index = node.getNumber();
            bits2.set(index);

        } else {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                assignCladeHeights(child, clades, bits2);
            }
            Clade c = new Clade(bits2, tree.getNodeHeight(node));
            if (clades.containsKey(c)) {
                tree.setNodeHeight(node, clades.get(c));
                clades.remove(c);
            }
        }

        // add my bit set to the bit set I was given
        // this is needed for adding all children clades together
        if (bits != null) {
            bits.or(bits2);
        }
    }

    private double[] getRelativeNodeHeights(Tree tree) {
        int count = tree.getInternalNodeCount();
        double[] nodeHeights = new double[count];

        for (int i = 0; i < count; i++) {
            NodeRef node = tree.getInternalNode(i);
            NodeRef parent = tree.getParent(node);
            nodeHeights[i] = tree.getNodeHeight(node)
                    / tree.getNodeHeight(parent);
        }

        return nodeHeights;
    }

    private double[] getAbsoluteNodeHeights(Tree tree) {
        int count = tree.getInternalNodeCount();
        double[] nodeHeights = new double[count];

        for (int i = 0; i < count; i++) {
            NodeRef node = tree.getInternalNode(i);
            nodeHeights[i] = tree.getNodeHeight(node);
        }

        return nodeHeights;
    }

    private double[] getAbsoluteNodeHeights(List<Clade> clades) {
        double[] nodeHeights = new double[clades.size()];

        int count = 0;
        for (Clade c : clades) {
            nodeHeights[count] = c.getHeight();
            count++;
        }

        return nodeHeights;
    }

    private double getChanceForNodeHeights(double[] nodeHeights) {
        return getChanceOfPermuation(nodeHeights);
        // return getChanceOfUniformNodeHeights(tree.getRoot());
//		return probabilityEstimater.getChanceForNodeHeights(tree, likelihood, prior);
    }

    private double getChanceOfUniformNodeHeights(NodeRef parent) {
        double prob = 0.0;

        NodeRef leftChild = tree.getChild(parent, 0);
        NodeRef rightChild = tree.getChild(parent, 1);

        if (!tree.isExternal(leftChild)) {
            prob += Math.log(1.0 / tree.getNodeHeight(parent));
            prob += getChanceOfUniformNodeHeights(leftChild);
        }
        if (!tree.isExternal(rightChild)) {
            prob += Math.log(1.0 / tree.getNodeHeight(parent));
            prob += getChanceOfUniformNodeHeights(rightChild);
        }

        return prob;
    }

    private double getChanceOfPermuation(double[] nodeHeights) {
        List<NodeRef> nodes = new LinkedList<NodeRef>();

        NodeRef root = tree.getRoot();

        NodeRef leftChild = tree.getChild(root, 0);
        NodeRef rightChild = tree.getChild(root, 1);

        if (!tree.isExternal(leftChild)) {
            nodes.add(leftChild);
        }
        if (!tree.isExternal(rightChild)) {
            nodes.add(rightChild);
        }

        int pointer = nodeHeights.length - 2;
        double prob = 0.0;
        while (!nodes.isEmpty()) {
            int index = getHighestNode(nodes);
            prob += Math.log(1.0 / nodes.size());
            NodeRef n = nodes.remove(index);
            tree.setNodeHeight(n, nodeHeights[pointer]);
            pointer--;
            leftChild = tree.getChild(n, 0);
            rightChild = tree.getChild(n, 1);

            if (!tree.isExternal(leftChild)) {
                nodes.add(leftChild);
            }
            if (!tree.isExternal(rightChild)) {
                nodes.add(rightChild);
            }
        }

        return prob;
    }

    private int getHighestNode(List<NodeRef> nodes) {
        double maxHeight = 0;
        int index = 0;

        for (int i = 0; i < nodes.size(); i++) {
            NodeRef n = nodes.get(i);
            if (tree.getNodeHeight(n) > maxHeight) {
                maxHeight = tree.getNodeHeight(n);
                index = i;
            }
        }

        return index;
    }

    private double setNodeHeights(double[] nodeHeights) {
        // return setUniformNodeHeights(tree.getRoot());
        return assignPermutedNodeHeights(nodeHeights);
//		return probabilityEstimater.setNodeHeights(tree, likelihood, prior);
    }

    private double setUniformNodeHeights(NodeRef parent) {
        double prob = 0.0;

        NodeRef leftChild = tree.getChild(parent, 0);
        NodeRef rightChild = tree.getChild(parent, 1);

        if (!tree.isExternal(leftChild)) {
            double max = tree.getNodeHeight(parent);
            double height = max * MathUtils.nextDouble();
            tree.setNodeHeight(leftChild, height);
            prob += Math.log(1.0 / max);
            prob += setUniformNodeHeights(leftChild);
        }
        if (!tree.isExternal(rightChild)) {
            double max = tree.getNodeHeight(parent);
            double height = max * MathUtils.nextDouble();
            tree.setNodeHeight(rightChild, height);
            prob += Math.log(1.0 / max);
            prob += setUniformNodeHeights(rightChild);
        }

        return prob;
    }

    private double assignPermutedNodeHeights(double[] nodeHeights) {
        List<NodeRef> nodes = new LinkedList<NodeRef>();

        NodeRef root = tree.getRoot();

        NodeRef leftChild = tree.getChild(root, 0);
        NodeRef rightChild = tree.getChild(root, 1);

        if (!tree.isExternal(leftChild)) {
            nodes.add(leftChild);
        }
        if (!tree.isExternal(rightChild)) {
            nodes.add(rightChild);
        }

        int pointer = nodeHeights.length - 2;
        double prob = 0.0;
        while (!nodes.isEmpty()) {
            int index = MathUtils.nextInt(nodes.size());
            prob += Math.log(1.0 / nodes.size());
            NodeRef n = nodes.remove(index);
            tree.setNodeHeight(n, nodeHeights[pointer]);
            pointer--;
            leftChild = tree.getChild(n, 0);
            rightChild = tree.getChild(n, 1);

            if (!tree.isExternal(leftChild)) {
                nodes.add(leftChild);
            }
            if (!tree.isExternal(rightChild)) {
                nodes.add(rightChild);
            }
        }

        return prob;
    }

    private double setMissingNodeHeights(NodeRef node) {
        double prob = 0.0;

        // check if the node is external
        if (!tree.isExternal(node)) {

            // otherwise, call all children and add its taxon together to one
            // clade
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                setMissingNodeHeights(child);
            }
            double min = getMinNodeHeight(node);
            double max = getMaxNodeHeight(node);
            if (max <= min) {
                max = tree.getNodeHeight(tree.getRoot());
            }
            prob += Math.log(1.0 / (max - min));
            double height = min + MathUtils.nextDouble() * (max - min);
            tree.setNodeHeight(node, height);
        }

        return prob;
    }

    private double getMinNodeHeight(NodeRef node) {
        double min = Double.MAX_VALUE;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            double height = tree.getNodeHeight(child);
            if (height < min) {
                min = height;
            }
        }
        return min;
    }

    private double getMaxNodeHeight(NodeRef node) {
        return tree.getNodeHeight(tree.getParent(node));
    }

    private void fillInternalNodes(NodeRef node) {
        if (!tree.isExternal(node)) {
            internalNodes.add(node);
            int childCount = tree.getChildCount(node);
            for (int i = 0; i < childCount; i++) {
                fillInternalNodes(tree.getChild(node, i));
            }
        }
    }

    private void fillExternalNodes(NodeRef node) {
        if (!tree.isExternal(node)) {
            int childCount = tree.getChildCount(node);
            for (int i = 0; i < childCount; i++) {
                fillExternalNodes(tree.getChild(node, i));
            }
        } else {
            Integer i = node.getNumber();
            externalNodes.put(i, node);
        }
    }

    private OperatorSchedule getOperatorSchedule(TreeModel treeModel) {

        ExchangeOperator narrowExchange = new ExchangeOperator(
                ExchangeOperator.NARROW, treeModel, 10);
        ExchangeOperator wideExchange = new ExchangeOperator(
                ExchangeOperator.WIDE, treeModel, 3);
        SubtreeSlideOperator subtreeSlide = new SubtreeSlideOperator(treeModel,
                10.0, 1.0, true, false, false, false, CoercionMode.COERCION_ON);
        NNI nni = new NNI(treeModel, 10.0);
        WilsonBalding wilsonBalding = new WilsonBalding(treeModel, 3.0);
        FNPR fnpr = new FNPR(treeModel, 5.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(narrowExchange);
        schedule.addOperator(wideExchange);
        schedule.addOperator(subtreeSlide);
        schedule.addOperator(nni);
        schedule.addOperator(wilsonBalding);
        schedule.addOperator(fnpr);

        return schedule;
    }

    protected double doUnguidedOperation() {
        int index = schedule.getNextOperatorIndex();
        SimpleMCMCOperator operator = (SimpleMCMCOperator) schedule
                .getOperator(index);

        return operator.doOperation();
    }

    /**
     * @return the number of transitions since last call to reset().
     */
    public long getTransitions() {
        return transitions;
    }

    /**
     * Set the number of transitions since last call to reset(). This is used to
     * restore the state of the operator
     */
    public void setTransitions(int transitions) {
        this.transitions = transitions;
    }

    public double getTransistionProbability() {
        long accepted = getAcceptCount();
        long rejected = getRejectCount();
        long transition = getTransitions();
        return (double) transition / (double) (accepted + rejected);
    }

    public void reset() {
        super.reset();
        transitions = 0;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.50;
    }

    public double getMaximumAcceptanceLevel() {
        return 1.0;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.75;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 1.0;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.SimpleMCMCOperator#getOperatorName()
      */
    @Override
    public abstract String getOperatorName();

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
      */
    public abstract String getPerformanceSuggestion();

}
