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
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeChangedEvent;
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
 * @author Marc A. Suchard
 */
public class HierarchicalSubstitutionModel extends AbstractModel implements SubstitutionModel, TreeTraitProvider {

    public static final int BUNDLE_VERSION = 1;
    public static final String RATE_TRAIT = "rate";
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

    // The structure of the character tree, rebuilt after a change of its topology (Default provider only; the
    // parametric providers refuse such a change).
    private List<Integer> postOrderTraversal;
    private List<Integer> preOrderTraversal;

    private boolean structureKnown = false;
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

    private Aggregation aggregation = Aggregation.AUTO;

    // A multifurcating character tree: the binary skeleton (tree) with a collapse state. The upload is the
    // canonical resolution of the contracted tree; uploadTriples and uploadSources are the last upload's structure, so
    // a change of the skeleton that leaves the contracted tree as it is (a neutral move) changes nothing.
    private final CharacterTreeContraction contraction;
    private double[] uploadTriples;
    private int[] uploadSources;
    private double[] storedUploadTriples;
    private int[] storedUploadSources;

    // Logs the rate of each internal node of the character tree, as branch rates are logged on the time tree.
    private final TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper(RATE_TRAIT, new TreeTrait.D() {
        public String getTraitName() {
            return RATE_TRAIT;
        }

        public Intent getIntent() {
            return Intent.NODE;
        }

        public Double getTrait(Tree characterTree, NodeRef node) {
            if (characterTree != tree || characterTree.isExternal(node)
                    || (contraction != null && contraction.isCollapsed(node.getNumber()))) {
                return null;
            }
            return rate(node.getNumber());
        }
    });

    public HierarchicalSubstitutionModel(String name, HierarchicalRateProvider hierarchicalRateProvider,
                                         FrequencyModel frequency) {
        super(name);
        if (hierarchicalRateProvider.getFrequencyModel() != null && hierarchicalRateProvider.getFrequencyModel() != frequency) {
            throw new IllegalArgumentException("The rate provider and the substitution model must share one FrequencyModel");
        }
        this.hierarchicalRates = hierarchicalRateProvider;
        this.tree = hierarchicalRateProvider.getTreeModel();
        this.contraction = hierarchicalRateProvider.getContraction();
        if (contraction != null && contraction.getTree() != tree) {
            throw new IllegalArgumentException("the contraction must be of the provider's character tree");
        }
        this.frequency = frequency;
        this.dataType = frequency.getDataType();
        this.stateCount = tree.getExternalNodeCount();
        this.externalNodeCount = tree.getExternalNodeCount();

        if (dataType.getStateCount() != stateCount) {
            throw new IllegalArgumentException("Character tree has " + stateCount + " tips but the data type has "
                    + dataType.getStateCount() + " states");
        }

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
        this.bundleNumber = new int[tree.getNodeCount()];
        this.internalNodeOfBundleRank = new int[tree.getInternalNodeCount()];
        buildStructure();

        addModel(tree);
        tree.addModelRestoreListener(this);
        hierarchicalRates.addToModelListener(this);
        addModel(frequency);
        if (contraction != null) {
            if (contraction.getKeep() != null) {
                addVariable(contraction.getKeep());
            }
            recordUploadStructure();
        }
    }

    public CharacterTreeContraction getContraction() {
        return contraction;
    }

    /**
     * The number of internal nodes of the character tree as the model sees it: the contracted tree's when the tree
     * is multifurcating, otherwise the tree's (S - 1 for a binary tree).
     */
    public int getContractedInternalNodeCount() {
        return contraction != null ? contraction.getContractedInternalNodeCount() : tree.getInternalNodeCount();
    }

    /**
     * The normalised rate of a character-tree node (skeleton node number); 0 at a collapsed node.
     */
    public double getCladeRate(int nodeNumber) {
        return rate(nodeNumber);
    }

    private double rate(int nodeNumber) {
        if (contraction != null && contraction.isCollapsed(nodeNumber)) {
            return 0.0;
        }
        return hierarchicalRates.getRate(tree, nodeNumber);
    }

    private void recordUploadStructure() {
        uploadTriples = contraction.getTriples().clone();
        uploadSources = contraction.getLambdaSource().clone();
    }

    private boolean uploadStructureUnchanged() {
        return Arrays.equals(uploadTriples, contraction.getTriples())
                && Arrays.equals(uploadSources, contraction.getLambdaSource());
    }

    public HierarchicalSubstitutionModel(String name, HierarchicalRateProvider hierarchicalRateProvider) {
        this(name, hierarchicalRateProvider, frequencyOf(hierarchicalRateProvider));
    }

    private static FrequencyModel frequencyOf(HierarchicalRateProvider hierarchicalRateProvider) {
        FrequencyModel frequency = hierarchicalRateProvider.getFrequencyModel();
        if (frequency == null) {
            throw new IllegalArgumentException("rate provider " + hierarchicalRateProvider
                    + " has no FrequencyModel; use the three-argument constructor");
        }
        return frequency;
    }

    private void buildStructure() {
        postOrderTraversal = getPostOrderTraversalOrder(tree);
        List<Integer> reversedOrder = new ArrayList<>(postOrderTraversal);
        Collections.reverse(reversedOrder);
        preOrderTraversal = reversedOrder;

        getMRCAs();

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
        structureKnown = true;
    }

    private void ensureStructure() {
        if (!structureKnown) {
            buildStructure();
        }
    }

    /**
     * True for an event of the character tree that changes its topology: a tree edit is reported as a
     * TreeChangedEvent without a node parameter. A change of a node height or another node parameter (reported
     * as a parameter event and as the parameter itself) does not change anything the model reads.
     */
    public static boolean changesTopology(Object event) {
        if (event instanceof TreeChangedEvent) {
            final TreeChangedEvent treeEvent = (TreeChangedEvent) event;
            return !(treeEvent.isNodeParameterChanged() || treeEvent.isHeightChanged());
        }
        return false;
    }

    public Aggregation getAggregation() {
        return aggregation;
    }

    public void setAggregation(Aggregation aggregation) {
        if (aggregation == null) {
            throw new IllegalArgumentException("aggregation must not be null");
        }
        this.aggregation = aggregation;
    }

    public HierarchicalRateProvider getRateProvider() {
        return hierarchicalRates;
    }

    private void computePi() {
        ensureStructure();
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
                final double poissonRate = rate(node.getNumber());
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
            return AbstractHierarchicalRateProvider.stateOfCode(dataType, tree.getTaxonId(nodeNumber));
        }
        throw new RuntimeException("Should not be called for internal node on character tree.");
    }

    @Override
    public void getTransitionProbabilities(double distance, double[] matrix) {
        ensureStructure();
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
        ensureStructure();
        final int internalCount = stateCount - 1;
        out[0] = BUNDLE_VERSION * BUNDLE_VERSION_SHIFT + internalCount;
        if (contraction != null) {
            final double[] triples = contraction.getTriples();
            for (int rank = 0; rank < internalCount; rank++) {
                final int post = 1 + 3 * rank;
                final int pre = 1 + 3 * internalCount + 3 * (internalCount - 1 - rank);
                for (int c = 0; c < 3; c++) {
                    out[post + c] = triples[3 * rank + c];
                    out[pre + c] = triples[3 * rank + c];
                }
            }
            return;
        }
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
        ensureStructure();
        getCachedValues(distance);
        final int internalCount = stateCount - 1;
        if (contraction != null) {
            final int[] sources = contraction.getLambdaSource();
            for (int rank = 0; rank < internalCount; rank++) {
                out[offset + rank] = sources[rank] >= 0 ? H[sources[rank]] : 0.0;
            }
        } else {
            for (int rank = 0; rank < internalCount; rank++) {
                out[offset + rank] = H[internalNodeOfBundleRank[rank]];
            }
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

    private static final double GENERATOR_TOLERANCE = 8.0 * Math.ulp(1.0);

    /**
     * The generator has non-negative off-diagonal entries iff the clade exchangeability
     * sigma(m) = sum over the path from m to the root of lambda_k / pi_k is non-negative at every internal node
     * (docs/generator-constraints). Negative clade rates are allowed as long as this holds.
     */
    public boolean isValidGenerator() {
        return getInvalidGeneratorDescription() == null;
    }

    public String getInvalidGeneratorDescription() {
        computePi();
        final double[] sigma = new double[tree.getNodeCount()];
        final double[] scale = new double[tree.getNodeCount()];
        for (int nodeNumber : preOrderTraversal) {
            final NodeRef node = tree.getNode(nodeNumber);
            if (tree.isExternal(node)) {
                continue;
            }
            final double term = rate(nodeNumber) / pi[nodeNumber];
            final NodeRef parent = tree.getParent(node);
            sigma[nodeNumber] = (parent == null ? 0.0 : sigma[parent.getNumber()]) + term;
            scale[nodeNumber] = (parent == null ? 0.0 : scale[parent.getNumber()]) + Math.abs(term);
            if (!(sigma[nodeNumber] >= -GENERATOR_TOLERANCE * scale[nodeNumber])) {
                return "clade " + cladeName(nodeNumber) + " has exchangeability " + sigma[nodeNumber]
                        + " < 0, so the hierarchical rates do not give a valid generator";
            }
        }
        return null;
    }

    public boolean hasNonNegativeRates() {
        return getNegativeRateDescription() == null;
    }

    public String getNegativeRateDescription() {
        final double[] lambdas = getBeagleModelLambdas();
        for (int rank = 0; rank < stateCount - 1; rank++) {
            if (!(lambdas[1 + rank] >= 0.0)) {
                final int node = contraction != null ? contraction.getLambdaSource()[rank] : internalNodeOfBundleRank[rank];
                return "clade " + cladeName(node) + " has rate " + lambdas[1 + rank];
            }
        }
        return null;
    }

    private String cladeName(int nodeNumber) {
        final List<String> codes = new ArrayList<>();
        collectTipIds(tree.getNode(nodeNumber), codes);
        Collections.sort(codes);
        final StringBuilder sb = new StringBuilder();
        for (String code : codes) {
            sb.append(code);
        }
        return sb.toString();
    }

    private void collectTipIds(NodeRef node, List<String> codes) {
        if (tree.isExternal(node)) {
            codes.add(tree.getTaxonId(node.getNumber()));
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                collectTipIds(tree.getChild(node, i), codes);
            }
        }
    }

    private void computeBeagleModel() {
        if (beagleModelKnown) {
            return;
        }
        final int internalCount = stateCount - 1;
        if (beagleModelTriples == null) {
            // the plugin reads 3 (S - 1) triples and ignores the inverse eigenvectors; JNI copies whole arrays, so
            // S x S arrays would cost two S^2 copies per model upload (2 x 134 MB at S = 4096)
            beagleModelTriples = new double[Math.max(1, 3 * internalCount)];
            beagleModelLambdas = new double[stateCount];
            beagleModelUnused = new double[1];
        }
        final double[] header = new double[getBundleHeaderLength()];
        getBundleHeader(header);
        System.arraycopy(header, 1, beagleModelTriples, 0, 3 * internalCount);
        beagleModelLambdas[0] = BUNDLE_VERSION * BUNDLE_VERSION_SHIFT + internalCount;
        if (contraction != null) {
            final int[] sources = contraction.getLambdaSource();
            for (int rank = 0; rank < internalCount; rank++) {
                beagleModelLambdas[1 + rank] = sources[rank] >= 0 ? rate(sources[rank]) : 0.0;
            }
        } else {
            for (int rank = 0; rank < internalCount; rank++) {
                beagleModelLambdas[1 + rank] = hierarchicalRates.getRate(tree, internalNodeOfBundleRank[rank]);
            }
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
            final double rate = rate(node.getNumber());
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
        if (model == tree) {
            if (!changesTopology(object)) {
                return;
            }
            structureKnown = false;
            if (contraction != null) {
                contraction.topologyChanged();
                if (uploadStructureUnchanged()) {
                    return;
                }
                recordUploadStructure();
            }
        } else if (hierarchicalRates.isCollapsedRateEvent(model, index)) {
            return;
        }
        eigenDecompositionKnown = false;
        beagleModelKnown = false;
        fireModelChanged();
    }

    @Override
    protected void storeState() {
        if (contraction != null) {
            storedUploadTriples = uploadTriples.clone();
            storedUploadSources = uploadSources.clone();
        }
    }

    @Override
    protected void restoreState() {
        eigenDecompositionKnown = false;
        beagleModelKnown = false;
        if (contraction != null) {
            uploadTriples = storedUploadTriples;
            uploadSources = storedUploadSources;
        }
    }

    @Override
    public void modelRestored(Model model) {
        if (model == tree) {
            structureKnown = false;
            if (contraction != null) {
                contraction.topologyChanged();
            }
        }
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (contraction != null && variable == contraction.getKeep()) {
            recordUploadStructure();
        }
        eigenDecompositionKnown = false;
        beagleModelKnown = false;
        fireModelChanged();
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    /**
     * Whether the BEAGLE plugin may aggregate the unobserved states of each site pattern: AUTO lets the plugin
     * decide, ON and OFF force it. All hierarchical models of one likelihood must agree.
     */
    public enum Aggregation {
        AUTO, ON, OFF;

        public static Aggregation parse(String value) {
            for (Aggregation a : values()) {
                if (a.name().equalsIgnoreCase(value)) {
                    return a;
                }
            }
            throw new IllegalArgumentException("aggregate must be auto, on or off; found '" + value + "'");
        }
    }

    /**
     * The states a tip observes at one site pattern, as the BEAGLE plugin reads its partials: none when every entry
     * is equal (a gap, or a constant uncertain column), otherwise the states with a non-zero entry, or every state
     * when no entry is zero.
     */
    public static void markObservedStates(double[] partials, int offset, int stateCount, boolean[] observed) {
        final double first = partials[offset];
        boolean constant = true;
        boolean zero = false;
        for (int i = 0; i < stateCount; i++) {
            final double v = partials[offset + i];
            constant = constant && v == first;
            zero = zero || v == 0.0;
        }
        if (constant) {
            return;
        }
        for (int i = 0; i < stateCount; i++) {
            if (!zero || partials[offset + i] != 0.0) {
                observed[i] = true;
            }
        }
    }

    /**
     * What the plugin's aggregation sees for one alignment, given the first 3 (S - 1) entries of
     * {@link #getBeagleModelTriples} (the post-order triples): observed[k][x] says whether some tip observes state x at
     * pattern k (states path: codes below the state count; partials path: {@link #markObservedStates}). Each distinct
     * observed set is a class; its reduced state count S_c is 1 for an empty set and otherwise one more than the
     * number of character-tree nodes on the paths from the observed states to the root.
     */
    public static AggregationSummary summarizeAggregation(double[] postOrderTriples, int stateCount, int categoryCount,
                                                          boolean[][] observed) {
        final int nodeCount = 2 * stateCount - 1;
        final int[] parentOf = new int[nodeCount];
        Arrays.fill(parentOf, -1);
        for (int t = 0; t < 3 * (stateCount - 1); t += 3) {
            final int parent = (int) postOrderTriples[t];
            parentOf[(int) postOrderTriples[t + 1]] = parent;
            parentOf[(int) postOrderTriples[t + 2]] = parent;
        }
        final java.util.Map<String, Integer> reducedByClass = new java.util.HashMap<>();
        final boolean[] marked = new boolean[nodeCount];
        long sum = 0;
        int max = 0;
        for (boolean[] states : observed) {
            final StringBuilder key = new StringBuilder();
            for (int x = 0; x < stateCount; x++) {
                if (states[x]) {
                    key.append(x).append(',');
                }
            }
            Integer reduced = reducedByClass.get(key.toString());
            if (reduced == null) {
                Arrays.fill(marked, false);
                int count = 0;
                for (int x = 0; x < stateCount; x++) {
                    if (!states[x]) {
                        continue;
                    }
                    for (int m = parentOf[x]; m >= 0 && !marked[m]; m = parentOf[m]) {
                        marked[m] = true;
                        count++;
                    }
                }
                reduced = count + 1;
                reducedByClass.put(key.toString(), reduced);
            }
            sum += reduced;
            max = Math.max(max, reduced);
        }
        return new AggregationSummary(observed.length, reducedByClass.size(), stateCount, categoryCount, sum, max);
    }

    public static final class AggregationSummary {

        // the plugin's automatic rule (BeagleCPUHierarchicalImpl createInstance and buildClasses, default kernel)
        private static final int MIN_STATES = 9;
        private static final int FUSED_MAX_STATES = 128;
        private static final int FUSED_MIN_PATTERNS = 1000;

        public final int patternCount;
        public final int classCount;
        public final int stateCount;
        public final int categoryCount;
        public final long reducedStateSum;
        public final int maxReducedStates;

        AggregationSummary(int patternCount, int classCount, int stateCount, int categoryCount, long reducedStateSum,
                           int maxReducedStates) {
            this.patternCount = patternCount;
            this.classCount = classCount;
            this.stateCount = stateCount;
            this.categoryCount = categoryCount;
            this.reducedStateSum = reducedStateSum;
            this.maxReducedStates = maxReducedStates;
        }

        public double getMeanReducedStates() {
            return patternCount == 0 ? 0.0 : (double) reducedStateSum / patternCount;
        }

        public double getRatio() {
            return patternCount == 0 ? 1.0 : (double) reducedStateSum / ((double) patternCount * stateCount);
        }

        public double getThreshold() {
            return categoryCount <= 1 ? 0.08 : (categoryCount <= 3 ? 0.15 : 0.25);
        }

        public boolean automaticRuleAggregates() {
            final boolean fusedBand = stateCount <= FUSED_MAX_STATES && patternCount >= FUSED_MIN_PATTERNS;
            return stateCount >= MIN_STATES && !fusedBand && getRatio() < getThreshold();
        }

        public boolean aggregates(Aggregation aggregation) {
            switch (aggregation) {
                case ON:
                    return true;
                case OFF:
                    return false;
                default:
                    return automaticRuleAggregates();
            }
        }

        public String describe(Aggregation aggregation) {
            return String.format("%d patterns in %d classes of observed states, reduced state count mean %.2f and max %d "
                            + "of %d, r = %.4f (automatic threshold %.2f at %d categor%s); aggregate=%s: %s",
                    patternCount, classCount, getMeanReducedStates(), maxReducedStates, stateCount, getRatio(),
                    getThreshold(), categoryCount, categoryCount == 1 ? "y" : "ies",
                    aggregation.name().toLowerCase(), describeDecision(aggregation));
        }

        private String describeDecision(Aggregation aggregation) {
            if (aggregation == Aggregation.ON) {
                return "aggregated";
            }
            if (aggregation == Aggregation.OFF) {
                return "not aggregated";
            }
            if (stateCount < MIN_STATES) {
                return "not aggregated (too few states)";
            }
            if (stateCount <= FUSED_MAX_STATES && patternCount >= FUSED_MIN_PATTERNS) {
                return "not aggregated (the fused kernel applies)";
            }
            return getRatio() < getThreshold() ? "aggregated" : "not aggregated (r above the threshold)";
        }
    }

    public interface HierarchicalRateProvider {
        double getRate(TreeModel tree, int node);
        TreeModel getTreeModel();
        void addToModelListener(AbstractModel listener);

        default FrequencyModel getFrequencyModel() {
            return null;
        }

        /**
         * The collapse state of a multifurcating character tree, or null for a tree whose internal nodes are all
         * clades. A provider with a contraction returns rate 0 at its collapsed nodes.
         */
        default CharacterTreeContraction getContraction() {
            return null;
        }

        /**
         * True for an event of this provider that changes only the raw rate of a collapsed node, which changes
         * nothing the model reads.
         */
        default boolean isCollapsedRateEvent(Model model, int index) {
            return false;
        }

        class Default implements HierarchicalRateProvider {

            private final TreeParameterModel hierarchicalPoissonRates;
            private final TreeModel tree;
            private final Parameter rateParameter;
            private final FrequencyModel frequency;
            private final boolean normalised;
            private final int tips;
            private final int[] parentRate;
            private final int[][] childRates;
            private final int[][] cladeStates;
            private boolean structureKnown = false;
            private final double[] lastRates;
            private final double[] lastFrequencies;
            private double normalization = Double.NaN;
            private final CharacterTreeContraction contraction;

            public Default(Parameter rateParameter, TreeModel treeModel) {
                this(rateParameter, treeModel, null, false);
            }

            public Default(Parameter rateParameter, TreeModel treeModel, FrequencyModel frequency, boolean normalised) {
                this(rateParameter, treeModel, frequency, normalised, null);
            }

            /**
             * Rates on a multifurcating character tree held as a binary skeleton with a collapse state: one raw rate
             * per skeleton node, the effective rate of a collapsed node is 0. With a keep parameter (the topology is
             * sampled) the raw rates are bounded below by 0; with fixed flags the generator bounds apply. The
             * normalisation sums the kept clades in the canonical order of the contracted tree and each clade's
             * states in ascending order, so it does not depend on how the skeleton resolves a multifurcation.
             */
            public Default(Parameter rateParameter, CharacterTreeContraction contraction, FrequencyModel frequency,
                           boolean normalised) {
                this(rateParameter, contraction.getTree(), frequency, normalised, contraction);
            }

            private Default(Parameter rateParameter, TreeModel treeModel, FrequencyModel frequency, boolean normalised,
                            CharacterTreeContraction contraction) {
                if (contraction != null && contraction.getTree() != treeModel) {
                    throw new IllegalArgumentException("the contraction must be of the provider's character tree");
                }
                this.contraction = contraction;
                if (rateParameter.getDimension() != treeModel.getInternalNodeCount()) {
                    throw new IllegalArgumentException("ratesParameter has dimension " + rateParameter.getDimension()
                            + " but the character tree has " + treeModel.getInternalNodeCount() + " internal nodes");
                }
                if (normalised && frequency == null) {
                    throw new IllegalArgumentException("a normalised ratesParameter needs the FrequencyModel");
                }
                this.hierarchicalPoissonRates = new TreeParameterModel(treeModel,
                        new CompoundParameter("Default.Poisson.rates", new Parameter[]{
                                new Parameter.Default(treeModel.getExternalNodeCount(), 1.0),
                                rateParameter}), true);
                this.tree = treeModel;
                this.rateParameter = rateParameter;
                this.frequency = frequency;
                this.normalised = normalised;
                this.tips = treeModel.getExternalNodeCount();
                final int internal = treeModel.getInternalNodeCount();
                if (frequency != null) {
                    this.parentRate = new int[internal];
                    this.childRates = new int[internal][];
                    this.cladeStates = new int[internal][];
                    this.lastRates = new double[internal];
                    this.lastFrequencies = new double[frequency.getDataType().getStateCount()];
                    final ModelListener listener = new ModelListener() {
                        public void modelChangedEvent(Model model, Object object, int index) {
                            if (changesTopology(object)) {
                                structureKnown = false;
                            }
                        }

                        public void modelRestored(Model model) {
                            structureKnown = false;
                        }
                    };
                    treeModel.addModelListener(listener);
                    treeModel.addModelRestoreListener(listener);
                } else {
                    this.parentRate = null;
                    this.childRates = null;
                    this.cladeStates = null;
                    this.lastRates = null;
                    this.lastFrequencies = null;
                }
                rateParameter.addBounds(frequency == null || (contraction != null && contraction.getKeep() != null)
                        ? new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, rateParameter.getDimension())
                        : new GeneratorBounds());
            }

            @Override
            public CharacterTreeContraction getContraction() {
                return contraction;
            }

            @Override
            public boolean isCollapsedRateEvent(Model model, int index) {
                return contraction != null && model == hierarchicalPoissonRates && index >= tips
                        && index < tree.getNodeCount() && contraction.isCollapsed(index);
            }

            private double effectiveRate(int index) {
                if (contraction != null && contraction.isCollapsed(tips + index)) {
                    return 0.0;
                }
                return rateParameter.getParameterValue(index);
            }

            /**
             * Parents, internal children and clade states of the internal nodes, by rate index; rebuilt after a
             * change of the character tree's topology (the rates stay attached to the node numbers).
             */
            private void ensureStructure() {
                if (structureKnown) {
                    return;
                }
                for (int i = 0; i < parentRate.length; i++) {
                    final NodeRef node = tree.getNode(tips + i);
                    final NodeRef parent = tree.getParent(node);
                    parentRate[i] = parent == null ? -1 : parent.getNumber() - tips;
                    final List<Integer> children = new ArrayList<>();
                    for (int c = 0; c < tree.getChildCount(node); c++) {
                        final NodeRef child = tree.getChild(node, c);
                        if (!tree.isExternal(child)) {
                            children.add(child.getNumber() - tips);
                        }
                    }
                    childRates[i] = new int[children.size()];
                    for (int c = 0; c < children.size(); c++) {
                        childRates[i][c] = children.get(c);
                    }
                    final List<Integer> states = new ArrayList<>();
                    collectStates(tree, node, frequency.getDataType(), states);
                    cladeStates[i] = new int[states.size()];
                    for (int k = 0; k < states.size(); k++) {
                        cladeStates[i][k] = states.get(k);
                    }
                }
                structureKnown = true;
                normalization = Double.NaN;
            }

            /**
             * The exact lower limits of the clade rates (docs/generator-constraints): with the other rates held
             * fixed, lambda_m must keep sigma(n) >= 0 at m and at every internal node n below m, which gives
             * lambda_m >= -pi_m (sigma(parent(m)) + min(0, min_n sum_{k on the path from n up to a child of m}
             * lambda_k / pi_k)). All rates within their limits is exactly the valid-generator cone, which
             * contains negative rates (HKY with kappa < 1). Frequency moves are not covered by these limits.
             */
            private final class GeneratorBounds implements Bounds<Double> {

                private double cladeFrequency(int i) {
                    double sum = 0.0;
                    for (int state : cladeStates[i]) {
                        sum += frequency.getFrequency(state);
                    }
                    return sum;
                }

                private double term(int i) {
                    return effectiveRate(i) / cladeFrequency(i);
                }

                private double minimumBelow(int i, double above) {
                    double minimum = Double.POSITIVE_INFINITY;
                    for (int c : childRates[i]) {
                        final double sum = above + term(c);
                        minimum = Math.min(minimum, Math.min(sum, minimumBelow(c, sum)));
                    }
                    return minimum;
                }

                @Override
                public Double getUpperLimit(int dimension) {
                    return Double.POSITIVE_INFINITY;
                }

                @Override
                public Double getLowerLimit(int dimension) {
                    ensureStructure();
                    if (contraction != null && contraction.isCollapsed(tips + dimension)) {
                        return 0.0;
                    }
                    double sigmaParent = 0.0;
                    for (int k = parentRate[dimension]; k >= 0; k = parentRate[k]) {
                        sigmaParent += term(k);
                    }
                    final double below = Math.min(0.0, minimumBelow(dimension, 0.0));
                    return -cladeFrequency(dimension) * (sigmaParent + below);
                }

                @Override
                public int getBoundsDimension() {
                    return rateParameter.getDimension();
                }
            }

            private static void collectStates(TreeModel tree, NodeRef node, DataType dataType, List<Integer> states) {
                if (tree.isExternal(node)) {
                    states.add(AbstractHierarchicalRateProvider.stateOfCode(dataType, tree.getTaxonId(node.getNumber())));
                } else {
                    for (int i = 0; i < tree.getChildCount(node); i++) {
                        collectStates(tree, tree.getChild(node, i), dataType, states);
                    }
                }
            }

            public boolean isNormalised() {
                return normalised;
            }

            public double getNormalization() {
                if (!normalised) {
                    return 1.0;
                }
                ensureStructure();
                boolean changed = Double.isNaN(normalization);
                for (int i = 0; i < lastRates.length; i++) {
                    final double rate = effectiveRate(i);
                    if (rate != lastRates[i]) {
                        lastRates[i] = rate;
                        changed = true;
                    }
                }
                for (int i = 0; i < lastFrequencies.length; i++) {
                    final double f = frequency.getFrequency(i);
                    if (f != lastFrequencies[i]) {
                        lastFrequencies[i] = f;
                        changed = true;
                    }
                }
                if (changed && contraction != null) {
                    double mu = 0.0;
                    for (int node : contraction.getKeptInternalNodes()) {
                        double piClade = 0.0, sumSquares = 0.0;
                        for (int state : contraction.sortedStatesBelow(node)) {
                            piClade += lastFrequencies[state];
                            sumSquares += lastFrequencies[state] * lastFrequencies[state];
                        }
                        mu += lastRates[node - tips] * (piClade - sumSquares / piClade);
                    }
                    normalization = mu;
                } else if (changed) {
                    double mu = 0.0;
                    for (int i = 0; i < cladeStates.length; i++) {
                        double piClade = 0.0, sumSquares = 0.0;
                        for (int state : cladeStates[i]) {
                            piClade += lastFrequencies[state];
                            sumSquares += lastFrequencies[state] * lastFrequencies[state];
                        }
                        mu += lastRates[i] * (piClade - sumSquares / piClade);
                    }
                    normalization = mu;
                }
                return normalization;
            }

            @Override
            public double getRate(TreeModel tree, int node) {
                if (contraction != null && contraction.isCollapsed(node)) {
                    return 0.0;
                }
                final double rate = hierarchicalPoissonRates.getNodeValue(tree, tree.getNode(node));
                return normalised ? rate / getNormalization() : rate;
            }

            @Override
            public TreeModel getTreeModel() {
                return tree;
            }

            @Override
            public FrequencyModel getFrequencyModel() {
                return frequency;
            }

            @Override
            public void addToModelListener(AbstractModel listener) {
                listener.addModel(hierarchicalPoissonRates);
            }
        }
    }
}
