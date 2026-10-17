/*
 * MATBranchHistoryModelTest.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package test.dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.operators.AbstractTreeOperator;
import dr.evomodel.operators.FixedHeightSubtreePruneRegraftOperator;
import dr.evomodel.operators.SubtreeJumpOperator;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.MAT.*;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.math.MathUtils;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tests the branch history models of the MAT likelihood (single hit, integrated hits, explicit hits
 * and timed hits) and the resampling of explicit histories.
 */
public class MATBranchHistoryModelTest extends TraceCorrelationAssert {

    private static final String[] MODELS = BranchHistoryModel.NAMES;

    private SitePatterns patterns;
    private TreeModel treeModel;
    private HKY hky;
    private Parameter clockRate;
    private StrictClockBranchRates branchRates;

    public MATBranchHistoryModelTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        MathUtils.setSeed(1234);

        String[][] taxaSequence = new String[2][];
        taxaSequence[0] = PRIMATES_TAXON_SEQUENCE[0].clone();
        taxaSequence[1] = PRIMATES_TAXON_SEQUENCE[1].clone();
        taxaSequence[1][0] = mask(taxaSequence[1][0], 0, 40, '-');
        taxaSequence[1][1] = mask(taxaSequence[1][1], 20, 60, 'N');
        createAlignment(taxaSequence, Nucleotides.INSTANCE);

        treeModel = createPrimateTreeModel();
        patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.CONSTANT_ONLY);

        FrequencyModel frequencyModel = new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(new double[]{0.3, 0.2, 0.15, 0.35}));
        hky = new HKY(new Parameter.Default(4.0), frequencyModel);
        // a fast clock so that branches are long enough for multiple hits to matter
        clockRate = new Parameter.Default(3.0);
        branchRates = new StrictClockBranchRates(clockRate);
    }

    private static String mask(String sequence, int from, int to, char c) {
        StringBuilder sb = new StringBuilder(sequence);
        for (int i = from; i < to; i++) {
            sb.setCharAt(i, c);
        }
        return sb.toString();
    }

    private MATDataLikelihoodDelegate createDelegate(String model) {
        return new MATDataLikelihoodDelegate(treeModel, patterns, new HomogeneousBranchModel(hky), null, branchRates, null,
                BranchHistoryModel.create(model, 4));
    }

    // ------------------------------------------------------------------------------------------
    // Path machinery
    // ------------------------------------------------------------------------------------------

    private static class TestModel extends ExplicitHitsHistoryModel {
        TestModel(HKY hky) {
            super(4);
            double[] q = new double[16];
            hky.getInfinitesimalMatrix(q);
            double[] exit = new double[4];
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    if (a != b) {
                        exit[a] += q[a * 4 + b];
                    }
                }
            }
            setRates(hky, q, exit);
        }

        double untimedPathLogDensity(int[] states, int k, double length, double nu) {
            return pathLogDensity(states, new double[k + 1], k, length, nu);
        }

        int sample(int a, int b, double length, double nu, int[][] states) {
            return samplePath(a, b, length, nu, states, new double[1][]);
        }

        double[] transitionProbabilities(double distance) {
            return getTransitionProbabilities(distance).clone();
        }
    }

    public void testOrderedTimesIntegral() throws Exception {
        java.lang.reflect.Method method = AbstractBranchHistoryModel.class.getDeclaredMethod(
                "logOrderedTimesIntegral", double[].class, int.class, double.class);
        method.setAccessible(true);
        double length = 0.7;
        double[][] rateSets = {{1.0, 2.5, 0.4}, {1.3, 1.3, 1.3}, {0.2, 3.0, 1.1, 2.0}};
        for (double[] c : rateSets) {
            int k = c.length - 1;
            double expected = numericalOrderedIntegral(c, length, 400);
            double actual = Math.exp((Double) method.invoke(null, c, k, length));
            // the midpoint rule above is accurate to about 1e-5
            assertEquals(expected, actual, 1E-4 * Math.abs(expected));
        }
    }

    /**
     * Midpoint-rule integral over 0 < t_1 < ... < t_k < L (k = 2 or 3) of exp(-sum c_i (t_{i+1} - t_i)).
     */
    private static double numericalOrderedIntegral(double[] c, double length, int steps) {
        double h = length / steps;
        double sum = 0.0;
        if (c.length == 3) {
            for (int i = 0; i < steps; i++) {
                double t1 = (i + 0.5) * h;
                for (int j = i; j < steps; j++) {
                    double t2 = (j + 0.5) * h;
                    double weight = i == j ? 0.5 : 1.0;
                    sum += weight * Math.exp(-c[0] * t1 - c[1] * (t2 - t1) - c[2] * (length - t2));
                }
            }
            return sum * h * h;
        }
        int coarse = steps / 4;
        h = length / coarse;
        for (int i = 0; i < coarse; i++) {
            double t1 = (i + 0.5) * h;
            for (int j = i; j < coarse; j++) {
                double t2 = (j + 0.5) * h;
                for (int l = j; l < coarse; l++) {
                    double t3 = (l + 0.5) * h;
                    double weight = (i == j ? 0.5 : 1.0) * (j == l ? 0.5 : 1.0) * (i == j && j == l ? 2.0 / 3.0 : 1.0);
                    sum += weight * Math.exp(-c[0] * t1 - c[1] * (t2 - t1) - c[2] * (t3 - t2) - c[3] * (length - t3));
                }
            }
        }
        return sum * h * h * h;
    }

    /**
     * The untimed path densities, summed over all paths (up to 8 hits), give the transition probabilities.
     */
    public void testPathDensitiesSumToTransitionProbabilities() {
        TestModel model = new TestModel(hky);
        double length = 0.4;
        double nu = 1.5;
        double[] p = model.transitionProbabilities(nu * length);
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                double total = 0.0;
                for (int k = 0; k <= 8; k++) {
                    total += sumPaths(model, new int[k + 1], 0, k, a, b, length, nu);
                }
                assertEquals("P(" + a + "," + b + ")", p[a * 4 + b], total, 1E-7);
            }
        }
    }

    private double sumPaths(TestModel model, int[] states, int i, int k, int a, int b, double length, double nu) {
        if (i == 0) {
            states[0] = a;
            return k == 0 ? (a == b ? Math.exp(model.untimedPathLogDensity(states, 0, length, nu)) : 0.0) :
                    sumPaths(model, states, 1, k, a, b, length, nu);
        }
        double total = 0.0;
        for (int x = 0; x < 4; x++) {
            if (x == states[i - 1] || (i == k && x != b)) {
                continue;
            }
            states[i] = x;
            total += i == k ? Math.exp(model.untimedPathLogDensity(states, k, length, nu)) :
                    sumPaths(model, states, i + 1, k, a, b, length, nu);
        }
        return total;
    }

    /**
     * The uniformization sampler draws paths with their exact conditional probabilities.
     */
    public void testPathSampler() {
        TestModel model = new TestModel(hky);
        double length = 0.5;
        double nu = 2.0;
        int[][] fromTo = {{0, 0}, {0, 2}, {1, 3}};
        int samples = 100000;
        for (int[] ends : fromTo) {
            int a = ends[0], b = ends[1];
            double pab = model.transitionProbabilities(nu * length)[a * 4 + b];
            Map<String, Integer> counts = new HashMap<String, Integer>();
            Map<String, Double> expected = new HashMap<String, Double>();
            int[][] states = new int[1][];
            for (int i = 0; i < samples; i++) {
                int k = model.sample(a, b, length, nu, states);
                String key = java.util.Arrays.toString(java.util.Arrays.copyOf(states[0], k + 1));
                counts.merge(key, 1, Integer::sum);
                if (!expected.containsKey(key)) {
                    expected.put(key, Math.exp(model.untimedPathLogDensity(java.util.Arrays.copyOf(states[0], k + 1), k, length, nu)) / pab);
                }
            }
            for (Map.Entry<String, Double> e : expected.entrySet()) {
                double probability = e.getValue();
                if (probability * samples < 50) {
                    continue;
                }
                double observed = counts.get(e.getKey()) / (double) samples;
                double se = Math.sqrt(probability * (1 - probability) / samples);
                assertEquals(a + "->" + b + " path " + e.getKey(), probability, observed, 5 * se);
            }
        }
    }

    // ------------------------------------------------------------------------------------------
    // The models in the likelihood
    // ------------------------------------------------------------------------------------------

    /**
     * For every model, the per-site likelihoods sum to the (aggregated, cached) likelihood, and the
     * likelihood with the MAT summed out (likelihood minus the log full conditional of the MAT) does not
     * change when the whole MAT is redrawn from its full conditional.
     */
    public void testMarginalIsInvariantToResampling() {
        for (String model : MODELS) {
            MATDataLikelihoodDelegate delegate = createDelegate(model);
            TreeDataLikelihood likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
            MutationHistoryResampler resampler = new MutationHistoryResampler(delegate);

            double marginal = likelihood.getLogLikelihood() - resampler.getLogFullConditional();
            int multipleHitSites = 0;
            for (int i = 0; i < 20; i++) {
                likelihood.storeModelState();
                resampler.resampleAll();
                likelihood.makeDirty();
                double logL = likelihood.getLogLikelihood();
                assertEquals(model + ": site log likelihoods", logL, sum(delegate.getSiteLogLikelihoods()), 1E-6);
                assertEquals(model + ": marginal", marginal, logL - resampler.getLogFullConditional(), 1E-6);
                multipleHitSites += countMultipleHitSites(delegate.getMutationStore());
                likelihood.acceptModelState();
            }
            if (delegate.getHistoryModel().hasExplicitHistories()) {
                assertTrue(model + " should sample some multiple hits", multipleHitSites > 0);
            } else {
                assertEquals(0, multipleHitSites);
            }
        }
    }

    /**
     * The summed-out likelihoods of the models agree: integrated hits on the variable sites is the
     * standard likelihood, and the explicit models differ from it only in the constant sites (which
     * are assumed to have no hits, rather than summed over hidden hits).
     */
    public void testMarginalsAgreeAcrossModels() {
        double[] marginals = new double[MODELS.length];
        for (int m = 0; m < MODELS.length; m++) {
            MATDataLikelihoodDelegate delegate = createDelegate(MODELS[m]);
            TreeDataLikelihood likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
            marginals[m] = likelihood.getLogLikelihood() - new MutationHistoryResampler(delegate).getLogFullConditional();
        }

        // contribution of the constant sites (fixed states) under integrated hits and under no hits
        double[] pi = hky.getFrequencyModel().getFrequencies();
        double[] q = new double[16];
        hky.getInfinitesimalMatrix(q);
        double[] p = new double[16];
        double integratedConstant = 0.0;
        double noHitConstant = 0.0;
        for (int i = 0; i < 4; i++) {
            double w = patterns.getPatternWeight(i);
            integratedConstant += w * Math.log(pi[i]);
            noHitConstant += w * Math.log(pi[i]);
            for (int n = 0; n < treeModel.getNodeCount(); n++) {
                NodeRef node = treeModel.getNode(n);
                if (treeModel.isRoot(node)) {
                    continue;
                }
                double length = treeModel.getBranchLength(node) * branchRates.getBranchRate(treeModel, node);
                hky.getTransitionProbabilities(length, p);
                integratedConstant += w * Math.log(p[i * 4 + i]);
                noHitConstant += w * q[i * 4 + i] * length;
            }
        }

        int integrated = indexOf(IntegratedHitsHistoryModel.NAME);
        for (String model : new String[]{ExplicitHitsHistoryModel.NAME, TimedHitsHistoryModel.NAME}) {
            assertEquals(model, marginals[integrated] - integratedConstant + noHitConstant, marginals[indexOf(model)], 1E-6);
        }

        // integrated hits on the variable sites == the standard (Felsenstein) likelihood of those sites
        BeagleDataLikelihoodDelegate beagle = new BeagleDataLikelihoodDelegate(treeModel, patterns,
                new HomogeneousBranchModel(hky), new GammaSiteRateModel("noGamma"), false, false,
                PartialsRescalingScheme.DEFAULT, false, PreOrderSettings.getDefault());
        TreeDataLikelihood standard = new TreeDataLikelihood(beagle, treeModel, branchRates);
        standard.getLogLikelihood();
        double[] siteLogL = beagle.getSiteLogLikelihoods();
        double variable = 0.0;
        for (int i = 4; i < patterns.getPatternCount(); i++) {
            variable += siteLogL[i];
        }
        assertEquals(variable, marginals[integrated] - integratedConstant, 1E-6);
    }

    /**
     * SPR moves with history resampling keep the MAT consistent for every model, and rejected moves
     * restore it exactly.
     */
    public void testResamplingMoves() {
        for (String model : MODELS) {
            treeModel = createPrimateTreeModel();
            MATDataLikelihoodDelegate delegate = createDelegate(model);
            TreeDataLikelihood likelihood = new TreeDataLikelihood(delegate, treeModel, branchRates);
            MutationHistoryResampler resampler = new MutationHistoryResampler(delegate);
            List<AbstractTreeOperator> operators = new ArrayList<AbstractTreeOperator>();
            operators.add(new SubtreeLeapOperator(treeModel, 1.0, 0.05, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                    AdaptationMode.ADAPTATION_OFF, 0.234));
            operators.add(new SubtreeJumpOperator(treeModel, 1.0, 1.0, 0.234, true, AdaptationMode.ADAPTATION_OFF));
            operators.add(new FixedHeightSubtreePruneRegraftOperator(treeModel, 1.0));
            for (AbstractTreeOperator operator : operators) {
                operator.setSubtreeMoveHook(resampler);
            }

            likelihood.getLogLikelihood();
            for (int iteration = 0; iteration < 100; iteration++) {
                double logLBefore = likelihood.getLogLikelihood();
                likelihood.storeModelState();
                double logHR = operators.get(iteration % operators.size()).doOperation();
                assertFalse(model, Double.isNaN(logHR));
                double logL = likelihood.getLogLikelihood();
                assertEquals(model + ": site log likelihoods after a move", logL, sum(delegate.getSiteLogLikelihoods()), 1E-6);
                if (MathUtils.nextBoolean()) {
                    likelihood.restoreModelState();
                    likelihood.makeDirty();
                    assertEquals(model + ": restore", logLBefore, likelihood.getLogLikelihood(), 1E-8);
                } else {
                    likelihood.acceptModelState();
                }
            }
            assertEquals(model, 0, delegate.getRebuildCount());
        }
    }

    private static int indexOf(String model) {
        for (int i = 0; i < MODELS.length; i++) {
            if (MODELS[i].equals(model)) {
                return i;
            }
        }
        throw new IllegalArgumentException(model);
    }

    private static double sum(double[] x) {
        double total = 0.0;
        for (double v : x) {
            total += v;
        }
        return total;
    }

    private static int countMultipleHitSites(MutationStore store) {
        int count = 0;
        for (int n = 0; n < store.getNodeCount(); n++) {
            BranchMutations mutations = store.getBranchMutations(n);
            int k = 0;
            while (k < mutations.getSubstitutionCount()) {
                int end = mutations.getSiteEnd(k);
                if (end - k > 1) {
                    count++;
                }
                k = end;
            }
        }
        return count;
    }
}
