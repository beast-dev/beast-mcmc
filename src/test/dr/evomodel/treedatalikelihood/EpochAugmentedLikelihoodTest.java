/*
 * EpochAugmentedLikelihoodTest.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate.LikelihoodException;
import dr.evomodel.treedatalikelihood.EpochEvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.NodeOperation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Compares the likelihood of an epoch model calculated with degree-2 nodes at the epoch transition times against the
 * likelihood calculated by convolving the matrices of the epochs along each branch, through BEAGLE and along a
 * sequence of proposals that are accepted and rejected.
 *
 * @author Marc A Suchard
 */
public class EpochAugmentedLikelihoodTest extends TraceCorrelationAssert {

    private static final double TOLERANCE = 1E-9;

    public EpochAugmentedLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testAugmentedNodesAreUsed() {
        Fixture augmented = new Fixture(0, true, PartialsRescalingScheme.NONE, 1);
        Fixture convolved = new Fixture(0, false, PartialsRescalingScheme.NONE, 1);

        assertTrue(augmented.delegate.getEvolutionaryProcessDelegate() instanceof EpochEvolutionaryProcessDelegate);
        assertFalse(convolved.delegate.getEvolutionaryProcessDelegate() instanceof EpochEvolutionaryProcessDelegate);

        assertEquals(convolved.likelihood.getLogLikelihood(), augmented.likelihood.getLogLikelihood(),
                TOLERANCE * Math.abs(convolved.likelihood.getLogLikelihood()));
    }

    /**
     * The second child of a degree-2 node is a buffer of ones, which any transition matrix would leave unchanged, so
     * check that the matrix is the identity itself
     */
    public void testIdentityMatrix() {
        for (int categories : new int[]{1, 4}) {
            Fixture augmented = new Fixture(0, true, PartialsRescalingScheme.NONE, categories);
            augmented.likelihood.getLogLikelihood();

            final int identityIndex = ((EpochEvolutionaryProcessDelegate)
                    augmented.delegate.getEvolutionaryProcessDelegate()).getIdentityMatrixIndex();

            final int stateCount = 4;
            final double[] matrix = new double[categories * stateCount * stateCount];
            augmented.delegate.getBeagleInstance().getTransitionMatrix(identityIndex, matrix);

            for (int i = 0; i < matrix.length; ++i) {
                final int row = (i / stateCount) % stateCount;
                final int column = i % stateCount;
                assertEquals("element " + i + " of " + categories + " categories", row == column ? 1.0 : 0.0,
                        matrix[i], 1E-14);
            }
        }
    }

    /**
     * With the epochs convolved, an epoch branch model passes every change of the tree on to the likelihood, which
     * then updates all nodes, because the substitution models on a branch depend on the node heights. With degree-2
     * nodes at the transition times, only the nodes that change, and their ancestors, are updated.
     */
    public void testOnlyChangedNodesAreUpdated() {
        for (boolean augment : new boolean[]{true, false}) {
            Fixture fixture = new Fixture(1, augment, PartialsRescalingScheme.NONE, 1);
            fixture.likelihood.getLogLikelihood();
            final int all = fixture.delegate.lastNodeOperationCount;

            NodeRef cherry = null;
            for (int i = 0; i < fixture.tree.getInternalNodeCount(); ++i) {
                NodeRef node = fixture.tree.getInternalNode(i);
                if (fixture.tree.isExternal(fixture.tree.getChild(node, 0)) &&
                        fixture.tree.isExternal(fixture.tree.getChild(node, 1))) {
                    cherry = node;
                    break;
                }
            }
            assertNotNull(cherry);

            final double low = Math.max(fixture.tree.getNodeHeight(fixture.tree.getChild(cherry, 0)),
                    fixture.tree.getNodeHeight(fixture.tree.getChild(cherry, 1)));
            final double high = fixture.tree.getNodeHeight(fixture.tree.getParent(cherry));

            fixture.store();
            fixture.tree.setNodeHeight(cherry, 0.5 * (low + high));
            fixture.likelihood.getLogLikelihood();

            if (augment) {
                assertTrue(fixture.delegate.lastNodeOperationCount + " of " + all,
                        fixture.delegate.lastNodeOperationCount < all);
            } else {
                assertEquals(all, fixture.delegate.lastNodeOperationCount);
            }
        }
    }

    public void testWithProposals() {
        PartialsRescalingScheme[] schemes = {PartialsRescalingScheme.NONE, PartialsRescalingScheme.DYNAMIC,
                PartialsRescalingScheme.ALWAYS};

        for (int seed = 0; seed < 3; ++seed) {
            for (PartialsRescalingScheme scheme : schemes) {
                for (int categories : new int[]{1, 4}) {

                    Fixture augmented = new Fixture(seed, true, scheme, categories);
                    Fixture convolved = new Fixture(seed, false, scheme, categories);

                    final String label = "seed " + seed + ", " + scheme.getText() + ", " + categories + " categories";
                    assertSame(label, convolved, augmented);

                    for (int iteration = 0; iteration < 60; ++iteration) {

                        final double accepted = convolved.likelihood.getLogLikelihood();

                        augmented.store();
                        convolved.store();

                        augmented.propose();
                        convolved.propose();

                        assertSame(label + ", iteration " + iteration, convolved, augmented);

                        if (iteration % 3 == 0) {
                            augmented.likelihood.acceptModelState();
                            convolved.likelihood.acceptModelState();
                        } else {
                            augmented.likelihood.restoreModelState();
                            convolved.likelihood.restoreModelState();

                            assertSame(label + ", restored", convolved, augmented);
                            assertEquals(label + ", restored", accepted, augmented.likelihood.getLogLikelihood(),
                                    TOLERANCE * Math.abs(accepted));
                        }
                    }
                }
            }
        }
    }

    private void assertSame(String label, Fixture expected, Fixture actual) {
        final double logL = expected.likelihood.getLogLikelihood();
        assertEquals(label, logL, actual.likelihood.getLogLikelihood(), TOLERANCE * Math.abs(logL));
    }

    /**
     * Counts the node operations of the last evaluation
     */
    private static class CountingDelegate extends BeagleDataLikelihoodDelegate {

        int lastNodeOperationCount;

        CountingDelegate(Tree tree, PatternList patternList, BranchModel branchModel,
                         SiteRateModel siteRateModel, PartialsRescalingScheme scheme, PreOrderSettings settings) {
            super(tree, patternList, branchModel, siteRateModel, false, false, scheme, false, settings);
        }

        @Override
        public double calculateLikelihood(List<BranchOperation> branchOperations, List<NodeOperation> nodeOperations,
                                          int rootNodeNumber) throws LikelihoodException {
            lastNodeOperationCount = nodeOperations.size();
            return super.calculateLikelihood(branchOperations, nodeOperations, rootNodeNumber);
        }
    }

    /**
     * A tree over the taxa of the alignment, an epoch model with a substitution model for each epoch, and the
     * likelihood of them with the epochs handled by augmented nodes or by convolution. Two fixtures with the same
     * seed are identical and make the same proposals.
     */
    private class Fixture {

        final DefaultTreeModel tree;
        final Parameter epochTimes;
        final CountingDelegate delegate;
        final TreeDataLikelihood likelihood;

        private final Random proposals;

        Fixture(long seed, boolean augment, PartialsRescalingScheme scheme, int categories) {
            Random random = new Random(seed);
            proposals = new Random(seed + 1000);

            tree = new DefaultTreeModel(randomTree(random));

            final int epochCount = 3 + random.nextInt(5);
            final double[] times = new double[epochCount - 1];
            double time = 0.0;
            for (int i = 0; i < times.length; ++i) {
                time += 0.05 + 0.3 * random.nextDouble();
                times[i] = time;
            }
            epochTimes = new Parameter.Default(times);

            List<SubstitutionModel> models = new ArrayList<SubstitutionModel>();
            for (int i = 0; i < epochCount; ++i) {
                double[] frequencies = new double[4];
                double sum = 0.0;
                for (int j = 0; j < 4; ++j) {
                    frequencies[j] = 0.5 + random.nextDouble();
                    sum += frequencies[j];
                }
                for (int j = 0; j < 4; ++j) {
                    frequencies[j] /= sum;
                }
                models.add(new HKY(0.5 + 4.0 * random.nextDouble(),
                        new FrequencyModel(Nucleotides.INSTANCE, frequencies)));
            }

            EpochBranchModel branchModel = new EpochBranchModel(tree, models, epochTimes);
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);

            SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

            delegate = new CountingDelegate(tree, patterns, branchModel, siteRateModel, scheme,
                    new PreOrderSettings(false, false, false, false, false, false, augment));

            likelihood = new TreeDataLikelihood(delegate, tree, new DefaultBranchRateModel());
        }

        void store() {
            likelihood.storeModelState();
        }

        /**
         * Moves a node height or an epoch transition time
         */
        void propose() {
            if (proposals.nextInt(6) == 0) {
                int i = proposals.nextInt(epochTimes.getDimension());
                double low = i == 0 ? 0.0 : epochTimes.getParameterValue(i - 1);
                double high = i == epochTimes.getDimension() - 1 ? epochTimes.getParameterValue(i) + 0.5 :
                        epochTimes.getParameterValue(i + 1);
                epochTimes.setParameterValue(i, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
            } else {
                NodeRef node;
                do {
                    node = tree.getInternalNode(proposals.nextInt(tree.getInternalNodeCount()));
                } while (tree.isRoot(node));

                double low = Math.max(tree.getNodeHeight(tree.getChild(node, 0)),
                        tree.getNodeHeight(tree.getChild(node, 1)));
                double high = tree.getNodeHeight(tree.getParent(node));

                tree.setNodeHeight(node, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
            }
        }

        private SimpleTree randomTree(Random random) {
            List<SimpleNode> active = new ArrayList<SimpleNode>();
            for (int i = 0; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.0);
                active.add(tip);
            }

            while (active.size() > 1) {
                SimpleNode left = active.remove(random.nextInt(active.size()));
                SimpleNode right = active.remove(random.nextInt(active.size()));

                SimpleNode parent = new SimpleNode();
                parent.setHeight(Math.max(left.getHeight(), right.getHeight()) + 0.02 + 0.2 * random.nextDouble());
                parent.addChild(left);
                parent.addChild(right);
                active.add(parent);
            }

            return new SimpleTree(active.get(0));
        }
    }
}
