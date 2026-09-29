/*
 * EpochPreOrderTest.java
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

import beagle.Beagle;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.AugmentedNodeRegistry;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.EpochEvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Checks the pre-order partials and the branch-rate gradient of an epoch model with degree-2 nodes at the epoch
 * transition times.
 * <p>
 * At every node of the augmented tree, augmented nodes included, the pre-order and post-order partials give the
 * likelihood of each pattern. The gradient with respect to the rate of a branch, which is the sum over the segments
 * of the branch, is compared with central differences of the log likelihood.
 *
 * @author Marc A Suchard
 */
public class EpochPreOrderTest extends TraceCorrelationAssert {

    private static final int STATES = 4;

    public EpochPreOrderTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testBottomPartials() {
        for (boolean spectral : new boolean[]{false, true}) {
            for (int seed = 0; seed < 3; ++seed) {
                for (int categories : new int[]{1, 4}) {
                    Fixture f = new Fixture(seed, PartialsRescalingScheme.NONE, categories, spectral);
                    f.gradient.getGradientLogDensity(); // computes the pre-order partials

                    checkPartials(f, false);
                }
            }
        }
    }

    /**
     * Top pre-order partials are only computed exactly by the spectral BEAGLE representation. On the standard CPU
     * implementation they differ from this identity by about 1E-5 relative at every node, also without degree-2 nodes.
     */
    public void testTopPartials() {
        for (int seed = 0; seed < 3; ++seed) {
            for (int categories : new int[]{1, 4}) {
                Fixture f = new Fixture(seed, PartialsRescalingScheme.NONE, categories, true);
                f.computeTopPartials();

                checkPartials(f, true);
            }
        }
    }

    public void testBranchRateGradient() {
        PartialsRescalingScheme[] schemes = {PartialsRescalingScheme.NONE, PartialsRescalingScheme.DYNAMIC,
                PartialsRescalingScheme.ALWAYS};

        for (boolean spectral : new boolean[]{false, true}) {
            for (int seed = 0; seed < 3; ++seed) {
                for (PartialsRescalingScheme scheme : schemes) {
                    for (int categories : new int[]{1, 4}) {

                        checkGradient(new Fixture(seed, scheme, categories, spectral),
                                "seed " + seed + ", " + scheme.getText() + ", " + categories + " categories" +
                                        (spectral ? ", spectral" : ""));
                    }
                }
            }
        }
    }

    private void checkGradient(Fixture f, String label) {

        double[] gradient = f.gradient.getGradientLogDensity();
        assertEquals(f.rates.getDimension(), gradient.length);

        final double h = 1E-5;
        for (int i = 0; i < gradient.length; ++i) {
            final double rate = f.rates.getParameterValue(i);

            f.rates.setParameterValue(i, rate + h);
            final double up = f.likelihood.getLogLikelihood();
            f.rates.setParameterValue(i, rate - h);
            final double down = f.likelihood.getLogLikelihood();
            f.rates.setParameterValue(i, rate);

            final double numerical = (up - down) / (2 * h);
            assertEquals(label + ", branch " + i, numerical, gradient[i],
                    1E-5 * Math.max(1.0, Math.abs(numerical)));
        }
    }

    /**
     * At each node, the pre-order and post-order partials give the likelihood of each pattern in each category.
     * A bottom pre-order partial lies at the node, below the branch above it. A top pre-order partial lies above
     * the branch, so the post-order partial is propagated along the branch first.
     */
    private void checkPartials(Fixture f, boolean top) {
        final Beagle beagle = f.delegate.getBeagleInstance();
        final Tree tree = f.tree;
        final int size = f.categories * f.patternCount * STATES;
        final int preOrderOffset = f.delegate.getPartialBufferCount();

        final double[] pre = new double[size];
        final double[] post = new double[size];
        final double[] matrix = new double[f.categories * STATES * STATES];

        // the likelihood at the root
        final int root = tree.getRoot().getNumber();
        beagle.getPartials(preOrderOffset + root, Beagle.NONE, pre);
        beagle.getPartials(f.delegate.getPartialBufferIndex(root), Beagle.NONE, post);
        final double[] expected = new double[f.categories * f.patternCount];
        for (int i = 0; i < expected.length; ++i) {
            for (int s = 0; s < STATES; ++s) {
                expected[i] += pre[i * STATES + s] * post[i * STATES + s];
            }
            assertTrue(expected[i] > 0.0);
        }

        final AugmentedNodeRegistry registry = f.registry;
        final int[] chain = new int[registry.getBoundaryCount()];

        int checked = 0;
        int augmented = 0;
        for (int n = 0; n < tree.getNodeCount(); ++n) {
            if (tree.isRoot(tree.getNode(n))) {
                continue;
            }

            final int count = registry.copyChain(n, chain);
            for (int j = -1; j < count; ++j) {
                final int id = (j < 0) ? n : chain[j];

                beagle.getPartials(preOrderOffset + id, Beagle.NONE, pre);
                if (id < tree.getExternalNodeCount()) {
                    setTipPartials(f, id, post); // tips have states, not partials
                } else {
                    beagle.getPartials(f.delegate.getPartialBufferIndex(id), Beagle.NONE, post);
                }
                if (top) {
                    beagle.getTransitionMatrix(f.delegate.getEvolutionaryProcessDelegate().getMatrixIndex(id),
                            matrix);
                }

                for (int c = 0; c < f.categories; ++c) {
                    for (int p = 0; p < f.patternCount; ++p) {
                        final int offset = (c * f.patternCount + p) * STATES;

                        double sum = 0.0;
                        for (int s = 0; s < STATES; ++s) {
                            double value = post[offset + s];
                            if (top) {
                                value = 0.0;
                                for (int t = 0; t < STATES; ++t) {
                                    value += matrix[(c * STATES + s) * STATES + t] * post[offset + t];
                                }
                            }
                            sum += pre[offset + s] * value;
                        }

                        final double reference = expected[c * f.patternCount + p];
                        assertEquals("node " + id + (registry.isAugmented(id) ? " (augmented)" : ""), reference,
                                sum, 1E-9 * reference);
                    }
                }

                ++checked;
                if (j >= 0) {
                    ++augmented;
                }
            }
        }

        assertTrue(augmented > 0);
        assertEquals(tree.getNodeCount() - 1 + augmented, checked);
    }

    private static void setTipPartials(Fixture f, int tip, double[] post) {
        final int taxon = f.patterns.getTaxonIndex(f.tree.getTaxonId(tip));
        for (int c = 0; c < f.categories; ++c) {
            for (int p = 0; p < f.patternCount; ++p) {
                final int state = f.patterns.getPatternState(taxon, p);
                for (int s = 0; s < STATES; ++s) {
                    post[(c * f.patternCount + p) * STATES + s] = (state == s || state >= STATES) ? 1.0 : 0.0;
                }
            }
        }
    }

    private class Fixture {

        final DefaultTreeModel tree;
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;
        final Parameter rates;
        final BranchRateGradientForDiscreteTrait gradient;
        final AugmentedNodeRegistry registry;

        final int categories;
        final int patternCount;
        final SitePatterns patterns;

        Fixture(long seed, PartialsRescalingScheme scheme, int categories, boolean spectral) {
            Random random = new Random(seed);
            this.categories = categories;

            tree = new DefaultTreeModel(randomTree(random));

            final int epochCount = 3 + random.nextInt(5);
            final double[] times = new double[epochCount - 1];
            double time = 0.0;
            for (int i = 0; i < times.length; ++i) {
                time += 0.05 + 0.3 * random.nextDouble();
                times[i] = time;
            }
            final Parameter epochTimes = new Parameter.Default(times);

            List<SubstitutionModel> models = new ArrayList<SubstitutionModel>();
            for (int i = 0; i < epochCount; ++i) {
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
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);

            patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
            patternCount = patterns.getPatternCount();

            rates = new Parameter.Default(tree.getNodeCount() - 1, 1.0);
            for (int i = 0; i < rates.getDimension(); ++i) {
                rates.setParameterValue(i, 0.5 + random.nextDouble());
            }
            ArbitraryBranchRates branchRates = new ArbitraryBranchRates(tree, rates,
                    ArbitraryBranchRates.make(false, false, false), false);

            // pre-order partials and the derivative with respect to the branch rates, with augmented epoch nodes
            delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel, siteRateModel,
                    false, false, scheme, false,
                    new PreOrderSettings(true, true, false, false, spectral, false, true));

            likelihood = new TreeDataLikelihood(delegate, tree, branchRates);
            registry = ((EpochEvolutionaryProcessDelegate) delegate.getEvolutionaryProcessDelegate())
                    .getAugmentedNodeRegistry();

            gradient = new BranchRateGradientForDiscreteTrait("test", likelihood, delegate, rates, false);
        }

        /**
         * Computes top pre-order partials, which the branch-rate gradient does not use
         */
        void computeTopPartials() {
            AbstractBeagleGradientDelegate topDelegate = new AbstractBeagleGradientDelegate("top", tree, delegate) {

                @Override
                protected DiscretePartialsType getPreOrderType() {
                    return DiscretePartialsType.TOP;
                }

                @Override
                protected int getGradientLength() {
                    return 0;
                }

                @Override
                protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
                }

                @Override
                protected void constructTraits(Helper treeTraitHelper) {
                }
            };

            new ProcessSimulation(likelihood, topDelegate).cacheSimulatedTraits(null);
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
