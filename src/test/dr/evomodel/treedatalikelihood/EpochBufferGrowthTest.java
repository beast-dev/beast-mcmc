/*
 * EpochBufferGrowthTest.java
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
import beagle.BeagleInfo;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.AugmentedNodeRegistry;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.EpochEvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.evomodel.treedatalikelihood.discrete.LogCtmcRateGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Twin likelihoods of one epoch model with degree-2 nodes at the epoch transition times: one keeps BEAGLE buffers
 * for every augmented node the tree can have (beagle.augmented.preallocate), the other starts with none and grows
 * them on demand (beagleEnsureBufferCounts, BEAGLE 4.2.0 or later). Over scripted node-height, tip-date and
 * epoch-time proposals, accepted and rejected, a rejected proposal that grew the buffers among them, the twins give
 * identical (==) log likelihoods, log-rate and branch-rate gradients, and pre-order partials of every node, augmented
 * nodes included. Standard and spectral BEAGLE implementations, with and without pre-order, without rescaling and
 * with dynamic rescaling (whose first evaluation is repeated, after a rescaling exception).
 * <p>
 * The growing twin grows at least twice, and at most once per evaluation. With a BEAGLE older than 4.2.0, or with
 * -Dbeagle.augmented.preallocate=true (testFallback), both twins preallocate and the same comparisons hold.
 *
 * @author Marc A Suchard
 */
public class EpochBufferGrowthTest extends TraceCorrelationAssert {

    private static final int STATES = 4;
    private static final int RATES = STATES * (STATES - 1);
    private static final String PREALLOCATE = "beagle.augmented.preallocate";

    public EpochBufferGrowthTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testGrowingTwinEqualsPreallocated() {
        runAll(false);
    }

    /**
     * With growth turned off, as with a BEAGLE older than 4.2.0, neither twin grows and they still agree
     */
    public void testFallback() {
        final String previous = System.getProperty(PREALLOCATE);
        System.setProperty(PREALLOCATE, "true");
        try {
            runAll(true);
        } finally {
            restoreProperty(previous);
        }
    }

    private void runAll(boolean fallback) {
        for (boolean spectral : new boolean[]{false, true}) {
            for (boolean preOrder : new boolean[]{false, true}) {
                for (PartialsRescalingScheme scheme : new PartialsRescalingScheme[]{PartialsRescalingScheme.NONE,
                        PartialsRescalingScheme.DYNAMIC}) {
                    for (int seed = 0; seed < 2; ++seed) {
                        run(seed, spectral, preOrder, scheme, fallback);
                    }
                }
            }
        }
    }

    private void run(long seed, boolean spectral, boolean preOrder, PartialsRescalingScheme scheme,
                     boolean fallback) {
        final String label = "seed " + seed + (spectral ? ", spectral" : ", standard") +
                (preOrder ? ", pre-order" : "") + ", " + scheme.getText();

        Fixture preallocated = new Fixture(seed, spectral, preOrder, scheme, true);
        Fixture growing = new Fixture(seed, spectral, preOrder, scheme, false);

        final int capacity = preallocated.registry.getCapacity();
        assertEquals(label, capacity, preallocated.delegate.getAllocatedAugmentedCount());
        // the growing twin starts with no augmented buffers, unless growth is unavailable or turned off
        final boolean grows = growing.delegate.getAllocatedAugmentedCount() == 0;
        if (fallback) {
            assertFalse(label + ": grows although turned off", grows);
        } else {
            assertEquals(label + ": growth with BEAGLE " + BeagleInfo.getVersion(), isGrowthAvailable(), grows);
        }
        if (!grows) {
            assertEquals(label, capacity, growing.delegate.getAllocatedAugmentedCount());
        }

        // the first evaluation (repeated after a rescaling exception with dynamic rescaling) grows once
        compare(label + ", first evaluation", preallocated, growing, grows);

        // move the first transition time to just above the tips, so that almost every lineage crosses it, which
        // needs more augmented nodes: reject it, then make it again and accept it
        for (boolean accept : new boolean[]{false, true}) {
            preallocated.store();
            growing.store();
            preallocated.moveFirstTransitionTime();
            growing.moveFirstTransitionTime();
            compare(label + ", near the tips" + (accept ? ", accepted" : ", rejected"), preallocated, growing, grows);
            if (accept) {
                preallocated.accept();
                growing.accept();
            } else {
                preallocated.restore();
                growing.restore();
                compare(label + ", restored", preallocated, growing, grows);
            }
        }

        for (int iteration = 0; iteration < 30; ++iteration) {
            preallocated.store();
            growing.store();
            preallocated.propose();
            growing.propose();
            compare(label + ", iteration " + iteration, preallocated, growing, grows);
            if (iteration % 3 == 0) {
                preallocated.accept();
                growing.accept();
            } else {
                preallocated.restore();
                growing.restore();
                compare(label + ", iteration " + iteration + ", restored", preallocated, growing, grows);
            }
        }

        assertEquals(label, 0, preallocated.delegate.getGrowthEventCount());
        if (grows) {
            assertTrue(label + ": grew " + growing.delegate.getGrowthEventCount() + " times",
                    growing.delegate.getGrowthEventCount() >= 2);
            assertTrue(label, growing.delegate.getAllocatedAugmentedCount() <= capacity);
        } else {
            assertEquals(label, 0, growing.delegate.getGrowthEventCount());
        }
    }

    private static boolean isGrowthAvailable() { // beagleEnsureBufferCounts from BEAGLE 4.2.0
        final int[] version = BeagleInfo.getVersionNumbers();
        return version.length >= 2 && (version[0] > 4 || (version[0] == 4 && version[1] >= 2));
    }

    /**
     * Evaluates both twins and compares everything they compute; the growing twin grows at most once
     */
    private void compare(String label, Fixture expected, Fixture actual, boolean grows) {
        final int growthBefore = actual.delegate.getGrowthEventCount();
        final double expectedLogL = expected.likelihood.getLogLikelihood();
        final double actualLogL = actual.likelihood.getLogLikelihood();
        final int growth = actual.delegate.getGrowthEventCount() - growthBefore;

        assertTrue(label + ": grew " + growth + " times in one evaluation", growth <= (grows ? 1 : 0));
        assertTrue(label + ": log likelihood " + actualLogL + " instead of " + expectedLogL,
                Double.compare(expectedLogL, actualLogL) == 0);
        assertTrue(label + ": augmented buffers " + actual.delegate.getAllocatedAugmentedCount() + " for " +
                        actual.registry.getHighWaterMark() + " nodes used",
                actual.delegate.getAllocatedAugmentedCount() >= actual.registry.getHighWaterMark());

        if (expected.branchRateGradient != null) {
            assertTrue(label + ": branch-rate gradients differ",
                    Arrays.equals(expected.branchRateGradient.getGradientLogDensity(),
                            actual.branchRateGradient.getGradientLogDensity()));
            for (int m = 0; m < expected.logRateGradients.size(); ++m) {
                assertTrue(label + ": log-rate gradients of model " + m + " differ",
                        Arrays.equals(expected.logRateGradients.get(m).getGradientLogDensity(),
                                actual.logRateGradients.get(m).getGradientLogDensity()));
            }
            assertEquals(label + ": growth while computing gradients", growthBefore + growth,
                    actual.delegate.getGrowthEventCount());
            comparePreOrderPartials(label, expected, actual);
        }
    }

    // the pre-order partials of every original and augmented node, as the last gradient left them
    private void comparePreOrderPartials(String label, Fixture expected, Fixture actual) {
        final int size = expected.categories * expected.patternCount * STATES;
        final double[] e = new double[size];
        final double[] a = new double[size];
        final int[] chain = new int[expected.registry.getBoundaryCount()];
        final int[] actualChain = new int[chain.length];
        for (int n = 0; n < expected.tree.getNodeCount(); ++n) {
            final int count = expected.registry.copyChain(n, chain);
            assertEquals(label + ": chains differ", count, actual.registry.copyChain(n, actualChain));
            for (int j = -1; j < count; ++j) {
                final int id = (j < 0) ? n : chain[j];
                assertEquals(label + ": chains differ", id, (j < 0) ? n : actualChain[j]);
                expected.delegate.getBeagleInstance().getPartials(expected.delegate.getPreOrderPartialIndex(id),
                        Beagle.NONE, e);
                actual.delegate.getBeagleInstance().getPartials(actual.delegate.getPreOrderPartialIndex(id),
                        Beagle.NONE, a);
                assertTrue(label + ": pre-order partials of " + id + " differ", Arrays.equals(e, a));
            }
        }
    }

    private static void restoreProperty(String previous) {
        if (previous == null) {
            System.clearProperty(PREALLOCATE);
        } else {
            System.setProperty(PREALLOCATE, previous);
        }
    }

    /**
     * A tree with tips at different times over the taxa of the alignment, three epoch transition times well above
     * most nodes, alternating log-rate substitution models, branch rates, and the likelihood with its gradients.
     * Two fixtures with the same seed are identical and make the same proposals.
     */
    private class Fixture {

        final DefaultTreeModel tree;
        final Parameter epochTimes;
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;
        final AugmentedNodeRegistry registry;
        final BranchRateGradientForDiscreteTrait branchRateGradient;
        final List<LogCtmcRateGradient> logRateGradients = new ArrayList<LogCtmcRateGradient>();

        final int categories = 2;
        final int patternCount;

        private final Random proposals;

        Fixture(long seed, boolean spectral, boolean preOrder, PartialsRescalingScheme scheme, boolean preallocate) {
            Random random = new Random(seed);
            proposals = new Random(seed + 1000);

            tree = new DefaultTreeModel(randomTree(random));
            final double root = tree.getNodeHeight(tree.getRoot());
            epochTimes = new Parameter.Default(new double[]{0.55 * root, 0.7 * root, 0.85 * root});

            List<LogRateSubstitutionModel> models = new ArrayList<LogRateSubstitutionModel>();
            for (String name : new String[]{"A", "B"}) {
                double[] frequencies = new double[STATES];
                double sum = 0.0;
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] = 0.5 + random.nextDouble();
                    sum += frequencies[j];
                }
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] /= sum;
                }
                double[] values = new double[RATES];
                for (int k = 0; k < RATES; ++k) {
                    values[k] = random.nextGaussian();
                }
                Parameter logRates = new Parameter.Default("logRates." + name, values);
                LogRateSubstitutionModel model = new LogRateSubstitutionModel("model." + name, Nucleotides.INSTANCE,
                        new FrequencyModel(Nucleotides.INSTANCE, frequencies),
                        new LogAdditiveCtmcRateProvider.DataAugmented.Basic(logRates.getId(), logRates));
                model.setNormalization(true);
                model.setScaleRatesByFrequencies(false);
                models.add(model);
            }
            List<SubstitutionModel> epochModels = new ArrayList<SubstitutionModel>();
            for (int e = 0; e <= epochTimes.getDimension(); ++e) {
                epochModels.add(models.get(e % 2));
            }

            EpochBranchModel branchModel = new EpochBranchModel(tree, epochModels, epochTimes);
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);
            SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
            patternCount = patterns.getPatternCount();

            Parameter rates = new Parameter.Default(tree.getNodeCount() - 1, 1.0);
            for (int i = 0; i < rates.getDimension(); ++i) {
                rates.setParameterValue(i, 0.5 + random.nextDouble());
            }
            ArbitraryBranchRates branchRates = new ArbitraryBranchRates(tree, rates,
                    ArbitraryBranchRates.make(false, false, false), false);

            // the property is read when the delegate is created
            final String previous = System.getProperty(PREALLOCATE);
            if (preallocate) {
                System.setProperty(PREALLOCATE, "true");
            }
            try {
                delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel, siteRateModel,
                        false, false, scheme, false,
                        new PreOrderSettings(preOrder, preOrder, false, false, spectral, false, true));
            } finally {
                restoreProperty(previous);
            }
            likelihood = new TreeDataLikelihood(delegate, tree, branchRates);
            registry = ((EpochEvolutionaryProcessDelegate) delegate.getEvolutionaryProcessDelegate())
                    .getAugmentedNodeRegistry();

            if (preOrder) {
                branchRateGradient = new BranchRateGradientForDiscreteTrait("branchRates", likelihood, delegate,
                        rates, false);
                for (LogRateSubstitutionModel model : models) {
                    logRateGradients.add(new LogCtmcRateGradient("logRates." + model.getId(), likelihood, delegate,
                            model, spectral ? ApproximationMode.EXACT_SPECTRAL : ApproximationMode.FIRST_ORDER,
                            false));
                }
            } else {
                branchRateGradient = null;
            }
        }

        void store() {
            likelihood.storeModelState();
        }

        void accept() {
            likelihood.acceptModelState();
        }

        void restore() {
            likelihood.restoreModelState();
        }

        void moveFirstTransitionTime() {
            double tips = 0.0;
            for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                tips = Math.max(tips, tree.getNodeHeight(tree.getExternalNode(i)));
            }
            epochTimes.setParameterValue(0, tips + 0.01 * (epochTimes.getParameterValue(1) - tips));
        }

        /**
         * Moves an internal node height, a tip date or an epoch transition time
         */
        void propose() {
            final int kind = proposals.nextInt(6);
            if (kind == 0) {
                final int i = proposals.nextInt(epochTimes.getDimension());
                final double low = i == 0 ? 0.0 : epochTimes.getParameterValue(i - 1);
                final double high = i == epochTimes.getDimension() - 1 ? epochTimes.getParameterValue(i) + 0.1 :
                        epochTimes.getParameterValue(i + 1);
                epochTimes.setParameterValue(i, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
            } else {
                NodeRef node;
                if (kind == 1) {
                    node = tree.getExternalNode(proposals.nextInt(tree.getExternalNodeCount()));
                } else {
                    do {
                        node = tree.getInternalNode(proposals.nextInt(tree.getInternalNodeCount()));
                    } while (tree.isRoot(node));
                }
                final double low = tree.isExternal(node) ? 0.0 : Math.max(tree.getNodeHeight(tree.getChild(node, 0)),
                        tree.getNodeHeight(tree.getChild(node, 1)));
                final double high = tree.getNodeHeight(tree.getParent(node));
                tree.setNodeHeight(node, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
            }
        }

        private SimpleTree randomTree(Random random) {
            List<SimpleNode> active = new ArrayList<SimpleNode>();
            for (int i = 0; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.05 * random.nextDouble());
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
