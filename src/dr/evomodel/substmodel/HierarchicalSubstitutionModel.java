/*
 * HierarchicalSubstitutionModel.java
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


package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;

/**
 * Hierarchical Rate Substitution Model
 *
 * @author Jeff Thorne
 * @author Tae-Kun Seo
 * @author Sangchul Choi
 * @author Xiang Ji
 */
public class HierarchicalSubstitutionModel extends AbstractModel implements SubstitutionModel {

    private TreeParameterModel hierarchicalRates;
    private TreeModel tree;
    private FrequencyModel frequency;
    private DataType dataType;
    private final int stateCount;
    private final BranchRateModel branchRateModel;

    private final double[] probNoEventsOnBranch;
    private final double[] probNoEventsFromRoot;
    private final double[] H;
    private final double[] pi;
    private final double[] HSums;
    private final int[] MRCAs;
    private final int[] nodeIndexForState;

    private final List<Integer> postOrderTraversal;
    private final List<Integer> preOrderTraversal;




    public HierarchicalSubstitutionModel(String name, TreeParameterModel hierarchicalRates, TreeModel tree,
                                         FrequencyModel frequency) {
        super(name);
        this.hierarchicalRates = hierarchicalRates;
        this.tree = tree;
        this.frequency = frequency;
        this.dataType = frequency.getDataType();
        this.stateCount = tree.getTaxonCount();
        this.branchRateModel = new StrictClockBranchRates(new Parameter.Default(tree.getNodeCount(), 1.0));
        this.postOrderTraversal = getPostOrderTraversalOrder(tree);
        this.preOrderTraversal = postOrderTraversal.reversed();

        this.probNoEventsOnBranch = new double[tree.getInternalNodeCount()];
        this.probNoEventsFromRoot = new double[tree.getNodeCount()];
        this.H = new double[tree.getNodeCount()];
        this.HSums = new double[tree.getNodeCount()];
        this.pi = new double[tree.getNodeCount()];
        this.MRCAs = new int[tree.getExternalNodeCount() * (tree.getExternalNodeCount() - 1) / 2];
        this.nodeIndexForState = new int[tree.getExternalNodeCount()];
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            nodeIndexForState[i] = getState(i);
        }

        getMRCAs();
        addModel(tree);
        addModel(hierarchicalRates);
        addModel(frequency);

    }

    private void getCachedValues(double distance) {

        for (int i = 0; i < postOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(postOrderTraversal.get(i));
            if (tree.isExternal(node)) {
                pi[node.getNumber()] = frequency.getFrequency(node.getNumber());
            } else {
                double sum = 0;
                for (int j = 0; j < tree.getChildCount(node); j++) {
                    sum += pi[tree.getChild(node, j).getNumber()];
                }
                pi[node.getNumber()] = sum;
            }
        }

        for (int i = 0; i < preOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(preOrderTraversal.get(i));

            if (!tree.isExternal(node)) {
                final int index = preOrderTraversal.get(i) - tree.getExternalNodeCount();
                final double poissonRate = hierarchicalRates.getNodeValue(tree, node);
                probNoEventsOnBranch[index] = Math.exp(-poissonRate * distance);
                if (tree.isRoot(node)) {
                    probNoEventsFromRoot[node.getNumber()] = 1.0;
                } else {
                    probNoEventsFromRoot[node.getNumber()] = probNoEventsFromRoot[tree.getParent(node).getNumber()] * probNoEventsOnBranch[tree.getParent(node).getNumber() - tree.getExternalNodeCount()];
                }

                H[preOrderTraversal.get(i)] = probNoEventsFromRoot[node.getNumber()] * (1.0 - probNoEventsOnBranch[index]) / pi[node.getNumber()];
            } else {
                probNoEventsFromRoot[node.getNumber()] = probNoEventsFromRoot[tree.getParent(node).getNumber()] * probNoEventsOnBranch[tree.getParent(node).getNumber() - tree.getExternalNodeCount()];
                H[preOrderTraversal.get(i)] = probNoEventsFromRoot[node.getNumber()] / pi[node.getNumber()];
            }
        }

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(i + tree.getExternalNodeCount());
            final double poissonRate = hierarchicalRates.getNodeValue(tree, node);
            probNoEventsOnBranch[i] = Math.exp(-poissonRate * distance);


        }
    }

    private int getState(int nodeNumber) {
        NodeRef node = tree.getNode(nodeNumber);
        if (tree.isExternal(node)) {
            return dataType.getState(tree.getTaxonId(nodeNumber));
        }
        throw new RuntimeException("Should not be called for internal node on character tree.");
    }



    @Override
    public void getTransitionProbabilities(double distance, double[] matrix) {
        getCachedValues(distance);

        for (int i = 0; i < preOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(preOrderTraversal.get(i));

            if (tree.isRoot(node)) {
                HSums[node.getNumber()] = H[node.getNumber()];
            } else {
                HSums[node.getNumber()] = HSums[tree.getParent(node).getNumber()] + H[node.getNumber()];
            }

        }

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                final double piJ = pi[j];
                if (i == j) {
                    matrix[i * stateCount + j] = piJ * HSums[nodeIndexForState[i]];
                } else {
                    matrix[i * stateCount + j] = piJ * HSums[MRCAs[getPairNodeIndex(nodeIndexForState[i], nodeIndexForState[j])]];
                }
            }
        }
    }

    private List<Integer> getPostOrderTraversalOrder(TreeModel tree) {
        List<Integer> order = new ArrayList<>();
        postOrderTraversalNode(tree, tree.getRoot(), order);
        return order;
    }

    private void postOrderTraversalNode(TreeModel tree, NodeRef node, List<Integer> order) {
        int nodeNumber = node.getNumber();

        if (!tree.isExternal(node)) {
            NodeRef leftChild = tree.getChild(node, 0);
            NodeRef rightChild = tree.getChild(node, 1);

            postOrderTraversalNode(tree, leftChild, order);
            postOrderTraversalNode(tree, rightChild, order);

        }
        order.add(nodeNumber);

    }

    private void getMRCAs() {

        List<Integer>[] descendantNodeNumbers = new List[tree.getNodeCount()];
        for (int i = 0; i < postOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(postOrderTraversal.get(i));
            if (tree.isExternal(node)) {
                descendantNodeNumbers[node.getNumber()] = new ArrayList<>();
                descendantNodeNumbers[node.getNumber()].add(node.getNumber());
            } else {
                List<Integer> leftDescendantNodeNumbers = descendantNodeNumbers[tree.getChild(node, 0).getNumber()];
                List<Integer> rightDescendantNodeNumbers = descendantNodeNumbers[tree.getChild(node, 1).getNumber()];
                for (int j = 0; j < leftDescendantNodeNumbers.size(); j++) {
                    final int leftNodeNumber = leftDescendantNodeNumbers.get(j);
                    for (int k = 0; k < rightDescendantNodeNumbers.size(); k++) {
                        final int rightNodeNumber = rightDescendantNodeNumbers.get(k);
                        MRCAs[getPairNodeIndex(leftNodeNumber, rightNodeNumber)] = node.getNumber();
                    }
                }
                descendantNodeNumbers[node.getNumber()] = new ArrayList<>(leftDescendantNodeNumbers.size() + rightDescendantNodeNumbers.size());
                descendantNodeNumbers[node.getNumber()].addAll(leftDescendantNodeNumbers);
                descendantNodeNumbers[node.getNumber()].addAll(rightDescendantNodeNumbers);
            }
        }
    }

    private int getPairNodeIndex(int leftNodeNumber, int rightNodeNumber) {
        final int larger = leftNodeNumber > rightNodeNumber ? leftNodeNumber : rightNodeNumber;
        final int smaller = leftNodeNumber < rightNodeNumber ? leftNodeNumber : rightNodeNumber;
        return larger - smaller - 1 + (2 * stateCount - smaller - 1) * smaller / 2;
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        throw new RuntimeException("Should not be called");
    }

    @Override
    public FrequencyModel getFrequencyModel() {
        return frequency;
    }

    @Override
    public void getInfinitesimalMatrix(double[] matrix) {
        for (int i = 0; i < stateCount; i++) {
            double offDiagonalSum = 0;
            for (int j = 0; j < stateCount; j++) {
                if (i != j) {

                }
            }
        }
    }

    @Override
    public DataType getDataType() {
        return dataType;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }
}
