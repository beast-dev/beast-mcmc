/*
 * EpochLikelihoodTraversalTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.branchratemodel.TimeVaryingBranchRateModel.EpochTimeProvider;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.AugmentedNodeRegistry;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.NodeOperation;
import dr.evomodel.treedatalikelihood.EpochLikelihoodTraversal;
import dr.evomodel.treedatalikelihood.EpochSubstitutionModelDelegate;
import dr.evomodel.treedatalikelihood.TreeTraversal.TraversalType;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.*;

/**
 * Checks the operations of EpochLikelihoodTraversal by executing them in plain Java and comparing the resulting
 * likelihood with a direct calculation that multiplies the matrices of the epochs along each branch.
 *
 * @author Marc A Suchard
 */
public class EpochLikelihoodTraversalTest extends TestCase {

    private static final int STATES = 4;
    private static final double TOLERANCE = 1E-12;

    public EpochLikelihoodTraversalTest(String name) {
        super(name);
    }

    /**
     * ((A, B):1.0, C):3.0 with tips at 0 and transition times 0.5, 1.5, 2.0 and 4.0
     */
    public void testSmallTree() {
        for (TraversalType type : new TraversalType[]{TraversalType.POST_ORDER, TraversalType.REVERSE_LEVEL_ORDER}) {

            SimpleNode a = tip("A", 0.0);
            SimpleNode b = tip("B", 0.0);
            SimpleNode c = tip("C", 0.0);
            SimpleNode x = parent(1.0, a, b);
            SimpleNode root = parent(3.0, x, c);

            Scenario s = new Scenario(new Random(1), new SimpleTree(root), new double[]{0.5, 1.5, 2.0, 4.0}, type);
            s.evaluate();

            AugmentedNodeRegistry registry = s.registry;
            int nodeA = s.number("A");
            int nodeB = s.number("B");
            int nodeC = s.number("C");
            int nodeX = s.number(x);

            // the boundary above the root is never crossed
            assertEquals(1, registry.getChainLength(nodeA));
            assertEquals(1, registry.getChainLength(nodeB));
            assertEquals(2, registry.getChainLength(nodeX));
            assertEquals(3, registry.getChainLength(nodeC));

            assertEquals(0, registry.getFirstBoundary(nodeA));
            assertEquals(1, registry.getFirstBoundary(nodeX));
            assertEquals(0, registry.getFirstBoundary(nodeC));

            // the epoch of the segment above each original and augmented node
            assertEquals(0, registry.getMatrixEpoch(nodeA));
            assertEquals(1, registry.getMatrixEpoch(registry.getChainTop(nodeA)));
            assertEquals(1, registry.getMatrixEpoch(nodeX));
            assertEquals(3, registry.getMatrixEpoch(registry.getChainTop(nodeX)));
            assertEquals(0, registry.getMatrixEpoch(nodeC));
            assertEquals(3, registry.getMatrixEpoch(registry.getChainTop(nodeC)));

            // 2 internal nodes plus 1 + 1 + 2 + 3 augmented nodes; 2 + 2 + 3 + 4 segments
            assertEquals(9, s.lastNodeOperations.size());
            assertEquals(11, s.lastBranchOperations.size());

            // depth below the root of the augmented tree: X sits below 2 augmented nodes, C below 3, and the tips
            // below X (or C) and their own augmented node
            int[] chain = new int[4];

            assertEquals(3, levelOf(s.lastNodeOperations, nodeX));
            registry.copyChain(nodeX, chain);
            assertEquals(2, levelOf(s.lastNodeOperations, chain[0]));
            assertEquals(1, levelOf(s.lastNodeOperations, chain[1]));

            registry.copyChain(nodeC, chain);
            assertEquals(3, levelOf(s.lastNodeOperations, chain[0]));
            assertEquals(2, levelOf(s.lastNodeOperations, chain[1]));
            assertEquals(1, levelOf(s.lastNodeOperations, chain[2]));

            registry.copyChain(nodeA, chain);
            assertEquals(4, levelOf(s.lastNodeOperations, chain[0]));

            // the root reads the top augmented node of each child
            NodeOperation rootOp = null;
            for (NodeOperation op : s.lastNodeOperations) {
                if (op.getLevel() == 0) {
                    assertNull(rootOp);
                    rootOp = op;
                }
            }
            assertNotNull(rootOp);
            assertEquals(new HashSet<Integer>(Arrays.asList(registry.getChainTop(nodeX), registry.getChainTop(nodeC))),
                    new HashSet<Integer>(Arrays.asList(rootOp.getLeftChild(), rootOp.getRightChild())));

            // single-child operations
            int degree2 = 0;
            for (NodeOperation op : s.lastNodeOperations) {
                if (op.getRightChild() == NodeOperation.NO_CHILD) {
                    assertTrue(registry.isAugmented(op.getNodeNumber()));
                    ++degree2;
                }
            }
            assertEquals(7, degree2);

            assertEquals(s.reference(), s.likelihood, TOLERANCE);
        }
    }

    /**
     * A node exactly at a transition time is not crossed
     */
    public void testNodeAtTransitionTime() {
        Random random = new Random(2);

        SimpleNode a = tip("A", 0.0);
        SimpleNode b = tip("B", 0.0);
        SimpleNode c = tip("C", 0.0);
        SimpleNode x = parent(1.5, a, b);
        SimpleNode root = parent(3.0, x, c);

        Scenario s = new Scenario(random, new SimpleTree(root), new double[]{0.5, 1.5, 2.0}, TraversalType.POST_ORDER);
        s.evaluate();

        assertEquals(1, s.registry.getChainLength(s.number("A"))); // 0.5 only
        assertEquals(1, s.registry.getChainLength(s.number(x)));   // 2.0 only
        assertEquals(2, s.registry.getMatrixEpoch(s.number(x)));   // 1.5 starts epoch 2
        assertEquals(s.reference(), s.likelihood, TOLERANCE);
    }

    public void testRandomTreesWithAcceptAndReject() {
        for (int seed = 0; seed < 6; ++seed) {
            for (TraversalType type : new TraversalType[]{TraversalType.POST_ORDER, TraversalType.REVERSE_LEVEL_ORDER}) {

                Random random = new Random(seed);
                int tips = 4 + random.nextInt(20);
                int boundaries = 1 + random.nextInt(8);

                double[] times = new double[boundaries];
                double time = 0.0;
                for (int i = 0; i < boundaries; ++i) {
                    time += 0.05 + 0.4 * random.nextDouble();
                    times[i] = time;
                }

                Scenario s = new Scenario(random, randomTree(random, tips, seed % 2 == 1), times, type);
                s.evaluate();
                assertEquals(s.reference(), s.likelihood, TOLERANCE);

                for (int iteration = 0; iteration < 150; ++iteration) {

                    final double accepted = s.likelihood;
                    s.storeState();

                    if (random.nextInt(6) == 0) {
                        s.proposeTransitionTime();
                    } else {
                        s.proposeNodeHeight();
                    }

                    s.evaluate();
                    assertEquals("seed " + seed + " iteration " + iteration, s.reference(), s.likelihood, TOLERANCE);
                    s.checkRegistry();

                    if (random.nextBoolean()) {
                        s.restoreState();
                        s.evaluate();

                        // nothing to recompute after a rejection, and the old likelihood is back
                        assertEquals(0, s.lastNodeOperations.size());
                        assertEquals(0, s.lastBranchOperations.size());
                        assertEquals(accepted, s.likelihood, TOLERANCE);
                        assertEquals(s.reference(), s.likelihood, TOLERANCE);
                    }
                }
            }
        }
    }

    private static int levelOf(List<NodeOperation> operations, int nodeNumber) {
        for (NodeOperation op : operations) {
            if (op.getNodeNumber() == nodeNumber) {
                return op.getLevel();
            }
        }
        throw new IllegalArgumentException("No operation for node " + nodeNumber);
    }

    private static SimpleNode tip(String name, double height) {
        SimpleNode node = new SimpleNode();
        node.setTaxon(new Taxon(name));
        node.setHeight(height);
        return node;
    }

    private static SimpleNode parent(double height, SimpleNode left, SimpleNode right) {
        SimpleNode node = new SimpleNode();
        node.setHeight(height);
        node.addChild(left);
        node.addChild(right);
        return node;
    }

    private static SimpleTree randomTree(Random random, int tipCount, boolean heterochronous) {
        List<SimpleNode> active = new ArrayList<SimpleNode>();
        for (int i = 0; i < tipCount; ++i) {
            active.add(tip("t" + i, heterochronous ? 0.5 * random.nextDouble() : 0.0));
        }

        while (active.size() > 1) {
            SimpleNode left = active.remove(random.nextInt(active.size()));
            SimpleNode right = active.remove(random.nextInt(active.size()));
            active.add(parent(Math.max(left.getHeight(), right.getHeight()) + 0.05 + 0.5 * random.nextDouble(),
                    left, right));
        }

        return new SimpleTree(active.get(0));
    }

    /**
     * A tree, an epoch model, and a plain Java stand-in for BEAGLE that executes the operations.
     */
    private static class Scenario {

        final Random random;
        final DefaultTreeModel tree;
        final Parameter epochTimes;
        final List<SubstitutionModel> models = new ArrayList<SubstitutionModel>();
        final BranchRateModel rates = new StrictClockBranchRates(new Parameter.Default(1.7));
        final EpochSubstitutionModelDelegate delegate;
        final AugmentedNodeRegistry registry;
        final EpochLikelihoodTraversal traversal;
        final TraversalType type;

        final double[][] tipPartials;

        // "device" buffers; arrays are replaced and never modified so that a store is a shallow copy
        double[][] partials;
        double[][] matrices;

        List<NodeOperation> lastNodeOperations;
        List<BranchOperation> lastBranchOperations;
        double likelihood;

        Scenario(Random random, SimpleTree simpleTree, double[] transitionTimes, TraversalType type) {
            this.random = random;
            this.type = type;
            this.tree = new DefaultTreeModel(simpleTree);
            this.epochTimes = new Parameter.Default(transitionTimes);

            for (int i = 0; i <= transitionTimes.length; ++i) {
                double[] frequencies = new double[STATES];
                double sum = 0.0;
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] = 0.5 + random.nextDouble();
                    sum += frequencies[j];
                }
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] /= sum;
                }
                models.add(new HKY(0.5 + 4.0 * random.nextDouble(),
                        new FrequencyModel(Nucleotides.INSTANCE, frequencies)));
            }

            EpochBranchModel branchModel = new EpochBranchModel(tree, models, epochTimes);
            EpochTimeProvider provider = new EpochTimeProvider.ParameterWrapper(epochTimes);

            delegate = new EpochSubstitutionModelDelegate(tree, branchModel, provider);
            registry = delegate.getAugmentedNodeRegistry();
            traversal = new EpochLikelihoodTraversal(tree, rates, type, delegate);

            tipPartials = new double[tree.getExternalNodeCount()][STATES];
            for (double[] tipPartial : tipPartials) {
                tipPartial[random.nextInt(STATES)] = 1.0;
            }

            partials = new double[registry.getTotalNodeCount()][];
            matrices = new double[registry.getTotalNodeCount()][];
            for (int i = 0; i < tipPartials.length; ++i) {
                partials[i] = tipPartials[i];
            }
        }

        int number(String taxonId) {
            for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                if (tree.getNodeTaxon(tree.getExternalNode(i)).getId().equals(taxonId)) {
                    return tree.getExternalNode(i).getNumber();
                }
            }
            throw new IllegalArgumentException(taxonId);
        }

        int number(SimpleNode simpleNode) {
            // internal nodes are matched by height, which is unique in these trees
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                if (!tree.isExternal(tree.getNode(i)) && tree.getNodeHeight(tree.getNode(i)) == simpleNode.getHeight()) {
                    return i;
                }
            }
            throw new IllegalArgumentException("Unknown node");
        }

        /**
         * Like TreeDataLikelihood: collect the operations, execute them, and mark every node as up to date
         */
        void evaluate() {
            traversal.dispatchTreeTraversalCollectBranchAndNodeOperations();

            lastBranchOperations = new ArrayList<BranchOperation>(traversal.getBranchOperations());
            lastNodeOperations = new ArrayList<NodeOperation>(traversal.getNodeOperations());

            execute();
            traversal.setAllNodesUpdated();

            double[] rootPartials = partials[tree.getRoot().getNumber()];
            double[] frequencies = delegate.getRootStateFrequencies();
            likelihood = 0.0;
            for (int i = 0; i < STATES; ++i) {
                likelihood += frequencies[i] * rootPartials[i];
            }
        }

        private void execute() {

            for (BranchOperation op : lastBranchOperations) {
                double[] matrix = new double[STATES * STATES];
                models.get(registry.getMatrixEpoch(op.getBranchNumber())).getTransitionProbabilities(
                        op.getBranchLength(), matrix);
                matrices[op.getBranchNumber()] = matrix;
            }

            Map<Integer, Integer> level = new HashMap<Integer, Integer>();
            for (NodeOperation op : lastNodeOperations) {
                level.put(op.getNodeNumber(), op.getLevel());
            }

            int previousLevel = Integer.MAX_VALUE;
            for (NodeOperation op : lastNodeOperations) {

                // operations at one level depend only on the level below
                assertTrue(op.getLevel() >= 0);
                for (int child : new int[]{op.getLeftChild(), op.getRightChild()}) {
                    if (child != NodeOperation.NO_CHILD && level.containsKey(child)) {
                        assertEquals(op.getLevel() + 1, (int) level.get(child));
                    }
                }

                if (type == TraversalType.REVERSE_LEVEL_ORDER) {
                    assertTrue(op.getLevel() <= previousLevel);
                    previousLevel = op.getLevel();
                }

                double[] result = multiply(matrices[op.getLeftChild()], partials[op.getLeftChild()]);
                if (op.getRightChild() != NodeOperation.NO_CHILD) {
                    double[] other = multiply(matrices[op.getRightChild()], partials[op.getRightChild()]);
                    for (int i = 0; i < STATES; ++i) {
                        result[i] *= other[i];
                    }
                }
                partials[op.getNodeNumber()] = result;
            }
        }

        double reference() {
            double[] frequencies = delegate.getRootStateFrequencies();
            double[] rootPartials = referencePartials(tree.getRoot());

            double sum = 0.0;
            for (int i = 0; i < STATES; ++i) {
                sum += frequencies[i] * rootPartials[i];
            }
            return sum;
        }

        /**
         * Convolves the epochs along each branch directly from the heights
         */
        private double[] referencePartials(NodeRef node) {
            if (tree.isExternal(node)) {
                return tipPartials[node.getNumber()];
            }

            final double height = tree.getNodeHeight(node);
            final double[] transitionTimes = epochTimes.getParameterValues();

            double[] product = null;
            for (int c = 0; c < 2; ++c) {
                NodeRef child = tree.getChild(node, c);
                double[] vector = referencePartials(child);
                double rate = rates.getBranchRate(tree, child);

                double low = tree.getNodeHeight(child);
                for (int epoch = 0; epoch <= transitionTimes.length; ++epoch) {
                    double start = Math.max(low, epoch == 0 ? Double.NEGATIVE_INFINITY : transitionTimes[epoch - 1]);
                    double end = Math.min(height, epoch == transitionTimes.length ?
                            Double.POSITIVE_INFINITY : transitionTimes[epoch]);
                    if (end > start) {
                        double[] matrix = new double[STATES * STATES];
                        models.get(epoch).getTransitionProbabilities(rate * (end - start), matrix);
                        vector = multiply(matrix, vector);
                    }
                }

                if (product == null) {
                    product = vector;
                } else {
                    for (int i = 0; i < STATES; ++i) {
                        product[i] *= vector[i];
                    }
                }
            }
            return product;
        }

        private static double[] multiply(double[] matrix, double[] vector) {
            double[] result = new double[STATES];
            for (int i = 0; i < STATES; ++i) {
                for (int j = 0; j < STATES; ++j) {
                    result[i] += matrix[i * STATES + j] * vector[j];
                }
            }
            return result;
        }

        /**
         * Chains must partition the free list: no augmented node is lost or used twice
         */
        void checkRegistry() {
            Set<Integer> used = new HashSet<Integer>();
            int[] chain = new int[registry.getBoundaryCount()];
            for (int n = 0; n < tree.getNodeCount(); ++n) {
                int count = registry.copyChain(n, chain);
                for (int j = 0; j < count; ++j) {
                    assertTrue(registry.isAugmented(chain[j]));
                    assertTrue(used.add(chain[j]));
                    assertEquals(registry.getFirstBoundary(n) + j + 1, registry.getMatrixEpoch(chain[j]));
                }
            }
            assertEquals(registry.getCapacity(), used.size() + registry.getFreeCount());
        }

        // Moves, stored so that they can be undone

        private double[][] storedPartials;
        private double[][] storedMatrices;
        private NodeRef movedNode;
        private double movedNodeHeight;
        private int movedBoundary;
        private double movedBoundaryTime;

        void storeState() {
            delegate.storeState();

            storedPartials = partials.clone();
            storedMatrices = matrices.clone();
            movedNode = null;
            movedBoundary = -1;
        }

        void restoreState() {
            delegate.restoreState();

            partials = storedPartials;
            matrices = storedMatrices;
            if (movedNode != null) {
                tree.setNodeHeight(movedNode, movedNodeHeight);
            }
            if (movedBoundary >= 0) {
                epochTimes.setParameterValue(movedBoundary, movedBoundaryTime);
            }
        }

        void proposeNodeHeight() {
            NodeRef node;
            do {
                node = tree.getNode(random.nextInt(tree.getNodeCount()));
            } while (tree.isRoot(node));

            double low = 0.0;
            if (!tree.isExternal(node)) {
                low = Math.max(tree.getNodeHeight(tree.getChild(node, 0)), tree.getNodeHeight(tree.getChild(node, 1)));
            }
            double high = tree.getNodeHeight(tree.getParent(node));

            movedNode = node;
            movedNodeHeight = tree.getNodeHeight(node);

            tree.setNodeHeight(node, low + (high - low) * random.nextDouble());
            traversal.updateNodeAndChildren(node); // the branch above node and the branches above its children
        }

        void proposeTransitionTime() {
            int i = random.nextInt(epochTimes.getDimension());
            double low = i == 0 ? 0.0 : epochTimes.getParameterValue(i - 1);
            double high = i == epochTimes.getDimension() - 1 ? epochTimes.getParameterValue(i) + 1.0 :
                    epochTimes.getParameterValue(i + 1);

            movedBoundary = i;
            movedBoundaryTime = epochTimes.getParameterValue(i);

            epochTimes.setParameterValue(i, low + (high - low) * (0.01 + 0.98 * random.nextDouble()));
        }
    }
}
