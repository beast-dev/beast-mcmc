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
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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

    public static final int BUNDLE_VERSION = 1;
    private static final int BUNDLE_VERSION_SHIFT = 65536;

    private final HierarchicalRateProvider hierarchicalRates;
    private final TreeModel tree;
    private final FrequencyModel frequency;
    private final DataType dataType;
    private final int stateCount;
    private final int externalNodeCount;

    private final double[] probNoEventsOnBranch;
    private final double[] probNoEventsFromRoot;
    private final double[] H;
    private final double[] pi;
    private final double[] HSums;
    private final int[] MRCAs;
    private final int[] stateOfNode;
    private final int[] nodeOfState;
    private final List<Integer>[] descendantTips;

    private final List<Integer> postOrderTraversal;
    private final List<Integer> preOrderTraversal;

    private final boolean binary;
    private final int[] bundleNumber;
    private final int[] internalNodeOfBundleRank;

    private EigenDecomposition eigenDecomposition;
    private boolean eigenDecompositionKnown = false;

    // The model as BEAGLE's hierarchical representation reads it from an eigen buffer (see
    // BEAGLE_FLAG_HIERARCHICAL_REPRESENTATION in beagle.h): eigenValues[0] = version word,
    // eigenValues[1 .. S-1] = lambda per internal node in bundle numbering, eigenVectors[0 .. 3(S-1)) =
    // post-order triples as in getBundleHeader; the inverse-eigenvector array is unused.
    private double[] beagleModelTriples;
    private double[] beagleModelLambdas;
    private double[] beagleModelUnused;
    private boolean beagleModelKnown = false;

    public HierarchicalSubstitutionModel(String name, HierarchicalRateProvider hierarchicalRateProvider,
                                         FrequencyModel frequency) {
        super(name);
        this.hierarchicalRates = hierarchicalRateProvider;
        this.tree = hierarchicalRateProvider.getTreeModel();
        this.frequency = frequency;
        this.dataType = frequency.getDataType();
        this.stateCount = tree.getExternalNodeCount();
        this.externalNodeCount = tree.getExternalNodeCount();

        if (dataType.getStateCount() != stateCount) {
            throw new IllegalArgumentException("Character tree has " + stateCount + " tips but the data type has "
                    + dataType.getStateCount() + " states");
        }

        this.postOrderTraversal = getPostOrderTraversalOrder(tree);
        List<Integer> reversedOrder = new ArrayList<>(postOrderTraversal);
        Collections.reverse(reversedOrder);
        this.preOrderTraversal = reversedOrder;

        this.probNoEventsOnBranch = new double[tree.getInternalNodeCount()];
        this.probNoEventsFromRoot = new double[tree.getNodeCount()];
        this.H = new double[tree.getNodeCount()];
        this.HSums = new double[tree.getNodeCount()];
        this.pi = new double[tree.getNodeCount()];
        this.MRCAs = new int[stateCount * (stateCount - 1) / 2];

        this.stateOfNode = new int[externalNodeCount];
        this.nodeOfState = new int[stateCount];
        Arrays.fill(nodeOfState, -1);
        for (int i = 0; i < externalNodeCount; i++) {
            final int state = getState(i);
            if (state < 0 || state >= stateCount || nodeOfState[state] != -1) {
                throw new IllegalArgumentException("Character tree tip " + tree.getTaxonId(i)
                        + " is not a distinct state of " + dataType.getDescription());
            }
            stateOfNode[i] = state;
            nodeOfState[state] = i;
        }

        this.descendantTips = new List[tree.getNodeCount()];
        getMRCAs();

        boolean allBinary = true;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            if (tree.getChildCount(tree.getInternalNode(i)) != 2) {
                allBinary = false;
            }
        }
        this.binary = allBinary;

        this.bundleNumber = new int[tree.getNodeCount()];
        this.internalNodeOfBundleRank = new int[tree.getInternalNodeCount()];
        int rank = 0;
        for (int nodeNumber : postOrderTraversal) {
            if (nodeNumber < externalNodeCount) {
                bundleNumber[nodeNumber] = stateOfNode[nodeNumber];
            } else {
                bundleNumber[nodeNumber] = stateCount + rank;
                internalNodeOfBundleRank[rank] = nodeNumber;
                rank++;
            }
        }

        addModel(tree);
        hierarchicalRates.addToModelListener(this);
        addModel(frequency);
    }

    private void computePi() {
        for (int i = 0; i < postOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(postOrderTraversal.get(i));
            if (tree.isExternal(node)) {
                pi[node.getNumber()] = frequency.getFrequency(stateOfNode[node.getNumber()]);
            } else {
                double sum = 0;
                for (int j = 0; j < tree.getChildCount(node); j++) {
                    sum += pi[tree.getChild(node, j).getNumber()];
                }
                pi[node.getNumber()] = sum;
            }
        }
    }

    private void getCachedValues(double distance) {

        computePi();

        for (int i = 0; i < preOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(preOrderTraversal.get(i));

            if (!tree.isExternal(node)) {
                final int index = node.getNumber() - externalNodeCount;
                final double poissonRate = hierarchicalRates.getRate(tree, node.getNumber());
                probNoEventsOnBranch[index] = Math.exp(-poissonRate * distance);
                if (tree.isRoot(node)) {
                    probNoEventsFromRoot[node.getNumber()] = 1.0;
                } else {
                    probNoEventsFromRoot[node.getNumber()] = probNoEventsFromRoot[tree.getParent(node).getNumber()]
                            * probNoEventsOnBranch[tree.getParent(node).getNumber() - externalNodeCount];
                }

                H[node.getNumber()] = probNoEventsFromRoot[node.getNumber()]
                        * -Math.expm1(-poissonRate * distance) / pi[node.getNumber()];
            } else {
                probNoEventsFromRoot[node.getNumber()] = probNoEventsFromRoot[tree.getParent(node).getNumber()]
                        * probNoEventsOnBranch[tree.getParent(node).getNumber() - externalNodeCount];
                H[node.getNumber()] = probNoEventsFromRoot[node.getNumber()] / pi[node.getNumber()];
            }
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
            final int nodeI = nodeOfState[i];
            for (int j = 0; j < stateCount; j++) {
                final int nodeJ = nodeOfState[j];
                final double piJ = pi[nodeJ];
                if (i == j) {
                    matrix[i * stateCount + j] = piJ * HSums[nodeI];
                } else {
                    matrix[i * stateCount + j] = piJ * HSums[MRCAs[getPairNodeIndex(nodeI, nodeJ)]];
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
        for (int i = 0; i < tree.getChildCount(node); i++) {
            postOrderTraversalNode(tree, tree.getChild(node, i), order);
        }
        order.add(node.getNumber());
    }

    private void getMRCAs() {
        for (int i = 0; i < postOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(postOrderTraversal.get(i));
            List<Integer> descendants = new ArrayList<>();
            if (tree.isExternal(node)) {
                descendants.add(node.getNumber());
            } else {
                for (int a = 0; a < tree.getChildCount(node); a++) {
                    List<Integer> left = descendantTips[tree.getChild(node, a).getNumber()];
                    for (int b = a + 1; b < tree.getChildCount(node); b++) {
                        List<Integer> right = descendantTips[tree.getChild(node, b).getNumber()];
                        for (int leftNode : left) {
                            for (int rightNode : right) {
                                MRCAs[getPairNodeIndex(leftNode, rightNode)] = node.getNumber();
                            }
                        }
                    }
                    descendants.addAll(left);
                }
            }
            descendantTips[node.getNumber()] = descendants;
        }
    }

    private int getPairNodeIndex(int leftNodeNumber, int rightNodeNumber) {
        final int larger = leftNodeNumber > rightNodeNumber ? leftNodeNumber : rightNodeNumber;
        final int smaller = leftNodeNumber < rightNodeNumber ? leftNodeNumber : rightNodeNumber;
        return larger - smaller - 1 + (2 * stateCount - smaller - 1) * smaller / 2;
    }

    public boolean isBinary() {
        return binary;
    }

    public int getBundleHeaderLength() {
        return 1 + 6 * (stateCount - 1);
    }

    public int getBundleBlockLength() {
        return (stateCount - 1) + stateCount;
    }

    public int getBundleLength(int categoryCount) {
        return getBundleHeaderLength() + categoryCount * getBundleBlockLength();
    }

    public void getBundleHeader(double[] out) {
        if (!binary) {
            throw new IllegalStateException("The hierarchical BEAGLE representation needs a binary character tree");
        }
        final int internalCount = stateCount - 1;
        out[0] = BUNDLE_VERSION * BUNDLE_VERSION_SHIFT + internalCount;
        for (int rank = 0; rank < internalCount; rank++) {
            final NodeRef node = tree.getNode(internalNodeOfBundleRank[rank]);
            final int parent = bundleNumber[node.getNumber()];
            final int left = bundleNumber[tree.getChild(node, 0).getNumber()];
            final int right = bundleNumber[tree.getChild(node, 1).getNumber()];
            final int post = 1 + 3 * rank;
            final int pre = 1 + 3 * internalCount + 3 * (internalCount - 1 - rank);
            out[post] = parent;
            out[post + 1] = left;
            out[post + 2] = right;
            out[pre] = parent;
            out[pre + 1] = left;
            out[pre + 2] = right;
        }
    }

    public void fillBundleBlock(double distance, double[] out, int offset) {
        getCachedValues(distance);
        final int internalCount = stateCount - 1;
        for (int rank = 0; rank < internalCount; rank++) {
            out[offset + rank] = H[internalNodeOfBundleRank[rank]];
        }
        for (int state = 0; state < stateCount; state++) {
            out[offset + internalCount + state] = probNoEventsFromRoot[nodeOfState[state]];
        }
    }

    public void fillBundle(double edgeLength, double[] categoryRates, double[] out) {
        getBundleHeader(out);
        for (int l = 0; l < categoryRates.length; l++) {
            fillBundleBlock(edgeLength * categoryRates[l], out, getBundleHeaderLength() + l * getBundleBlockLength());
        }
    }

    public double[] getBeagleModelTriples() {
        computeBeagleModel();
        return beagleModelTriples;
    }

    public double[] getBeagleModelLambdas() {
        computeBeagleModel();
        return beagleModelLambdas;
    }

    public double[] getBeagleModelUnused() {
        computeBeagleModel();
        return beagleModelUnused;
    }

    private void computeBeagleModel() {
        if (beagleModelKnown) {
            return;
        }
        final int internalCount = stateCount - 1;
        if (beagleModelTriples == null) {
            beagleModelTriples = new double[stateCount * stateCount];
            beagleModelLambdas = new double[stateCount];
            beagleModelUnused = new double[stateCount * stateCount];
        }
        final double[] header = new double[getBundleHeaderLength()];
        getBundleHeader(header);
        System.arraycopy(header, 1, beagleModelTriples, 0, 3 * internalCount);
        beagleModelLambdas[0] = BUNDLE_VERSION * BUNDLE_VERSION_SHIFT + internalCount;
        for (int rank = 0; rank < internalCount; rank++) {
            beagleModelLambdas[1 + rank] = hierarchicalRates.getRate(tree, internalNodeOfBundleRank[rank]);
        }
        beagleModelKnown = true;
    }

    private double[] getCumulativeRates(boolean divideByPi) {
        computePi();
        final double[] cumulative = new double[tree.getNodeCount()];
        for (int i = 0; i < preOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(preOrderTraversal.get(i));
            if (tree.isExternal(node)) {
                continue;
            }
            final double rate = hierarchicalRates.getRate(tree, node.getNumber());
            final double term = divideByPi ? rate / pi[node.getNumber()] : rate;
            cumulative[node.getNumber()] = tree.isRoot(node) ? term : cumulative[tree.getParent(node).getNumber()] + term;
        }
        return cumulative;
    }

    @Override
    public EigenDecomposition getEigenDecomposition() {
        if (!eigenDecompositionKnown) {
            eigenDecomposition = computeEigenDecomposition();
            eigenDecompositionKnown = true;
        }
        return eigenDecomposition;
    }

    private EigenDecomposition computeEigenDecomposition() {
        final double[] cumulativeRates = getCumulativeRates(false);
        final double[] evec = new double[stateCount * stateCount];
        final double[] ievc = new double[stateCount * stateCount];
        final double[] eval = new double[stateCount];

        for (int state = 0; state < stateCount; state++) {
            evec[state * stateCount] = 1.0;
            ievc[state] = pi[nodeOfState[state]];
        }

        int column = 1;
        final double[] v = new double[stateCount];
        for (int i = 0; i < postOrderTraversal.size(); i++) {
            NodeRef node = tree.getNode(postOrderTraversal.get(i));
            if (tree.isExternal(node)) {
                continue;
            }
            final int childCount = tree.getChildCount(node);
            double leadingPi = 0.0;
            for (int r = 0; r < childCount - 1; r++) {
                leadingPi += pi[tree.getChild(node, r).getNumber()];
                final int next = tree.getChild(node, r + 1).getNumber();
                Arrays.fill(v, 0.0);
                for (int c = 0; c <= r; c++) {
                    for (int tipNode : descendantTips[tree.getChild(node, c).getNumber()]) {
                        v[stateOfNode[tipNode]] = 1.0 / leadingPi;
                    }
                }
                for (int tipNode : descendantTips[next]) {
                    v[stateOfNode[tipNode]] = -1.0 / pi[next];
                }
                double norm = 0.0;
                for (int state = 0; state < stateCount; state++) {
                    norm += pi[nodeOfState[state]] * v[state] * v[state];
                }
                for (int state = 0; state < stateCount; state++) {
                    evec[state * stateCount + column] = v[state];
                    ievc[column * stateCount + state] = pi[nodeOfState[state]] * v[state] / norm;
                }
                eval[column] = -cumulativeRates[node.getNumber()];
                column++;
            }
        }
        return new EigenDecomposition(evec, ievc, eval);
    }

    @Override
    public FrequencyModel getFrequencyModel() {
        return frequency;
    }

    @Override
    public void getInfinitesimalMatrix(double[] matrix) {
        final double[] sigma = getCumulativeRates(true);
        for (int i = 0; i < stateCount; i++) {
            final int nodeI = nodeOfState[i];
            double offDiagonalSum = 0.0;
            for (int j = 0; j < stateCount; j++) {
                if (i != j) {
                    final int nodeJ = nodeOfState[j];
                    final double rate = pi[nodeJ] * sigma[MRCAs[getPairNodeIndex(nodeI, nodeJ)]];
                    matrix[i * stateCount + j] = rate;
                    offDiagonalSum += rate;
                }
            }
            matrix[i * stateCount + i] = -offDiagonalSum;
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
        eigenDecompositionKnown = false;
        beagleModelKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        eigenDecompositionKnown = false;
        beagleModelKnown = false;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        eigenDecompositionKnown = false;
        beagleModelKnown = false;
        fireModelChanged();
    }

    public interface HierarchicalRateProvider {
        double getRate(TreeModel tree, int node);
        TreeModel getTreeModel();
        void addToModelListener(AbstractModel listener);

        class Default implements HierarchicalRateProvider {

            private final TreeParameterModel hierarchicalPoissonRates;
            private final TreeModel tree;

            public Default(Parameter rateParameter, TreeModel treeModel) {
                if (rateParameter.getDimension() != treeModel.getInternalNodeCount()) {
                    throw new RuntimeException("Incorrect parameter dimension");
                }
                this.hierarchicalPoissonRates = new TreeParameterModel(treeModel,
                        new CompoundParameter("Default.Poisson.rates", new Parameter[]{
                                new Parameter.Default(treeModel.getExternalNodeCount(), 1.0),
                                rateParameter}), true);
                this.tree = treeModel;

            }

            @Override
            public double getRate(TreeModel tree, int node) {
                return hierarchicalPoissonRates.getNodeValue(tree, tree.getNode(node));
            }

            @Override
            public TreeModel getTreeModel() {
                return tree;
            }

            @Override
            public void addToModelListener(AbstractModel listener) {
                listener.addModel(hierarchicalPoissonRates);
            }
        }
    }
}
