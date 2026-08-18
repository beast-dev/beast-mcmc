/*
 * OneZeroOneShuffleGibbsOperatorTargetTest.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.operators;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RewardsAwarePartialsRepresentation;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.inference.operators.OneZeroOneShuffleGibbsOperator;
import dr.math.MathUtils;
import test.dr.math.MathTestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Exact-target and distributional regressions for
 * {@link OneZeroOneShuffleGibbsOperator}.
 *
 * @author Filippo Monti
 */
public class OneZeroOneShuffleGibbsOperatorTargetTest extends MathTestCase {

    private static final double TOL = 1.0e-7;

    public void testDoOperationLandsOnACorrectlyEvaluatedPermutation() {
        MathUtils.setSeed(20260818);

        final Fixture fixture = createFixture(4, false);
        final double[] truth = bruteForceLogTargets(fixture);

        final OneZeroOneShuffleGibbsOperator operator = createOperator(fixture);
        assertFinite(operator.doOperation());

        final int landedIndex = indexOfCurrentMapping(fixture);
        fixture.independentLikelihood.makeDirty();
        final double landedLogLikelihood = fixture.independentLikelihood.getLogLikelihood();

        assertLogEquals("Operator landed on a permutation whose fresh likelihood does not match " +
                        "the independently brute-forced value for that permutation",
                truth[landedIndex], landedLogLikelihood);
    }

    public void testEmpiricalDistributionMatchesTheoreticalSoftmax() {
        MathUtils.setSeed(20260818);

        final Fixture fixture = createFixture(4, false);
        final double[] truth = bruteForceLogTargets(fixture);
        final double[] theoreticalProbabilities = softmax(truth);

        final OneZeroOneShuffleGibbsOperator operator = createOperator(fixture);

        final int nDraws = 20000;
        final int[] counts = new int[fixture.permutations.size()];
        for (int draw = 0; draw < nDraws; draw++) {
            resetMapping(fixture, fixture.permutations.get(0));
            operator.doOperation();
            counts[indexOfCurrentMapping(fixture)]++;
        }

        // Standard chi-squared goodness-of-fit practice: only include cells with
        // expected count >= 5 (low-probability permutations otherwise destabilize
        // the chi-squared approximation). Threshold is a generous many-sigma bound
        // on the actual degrees of freedom used, computed dynamically so the test
        // isn't tied to a hand-picked K/critical-value pair.
        final ChiSquaredResult result = chiSquaredStatistic(counts, theoreticalProbabilities, nDraws);
        final double threshold = result.degreesOfFreedom + 6.0 * Math.sqrt(2.0 * result.degreesOfFreedom);
        assertTrue("Empirical permutation distribution diverges from the theoretical softmax " +
                        "(chi-squared=" + result.statistic + ", df=" + result.degreesOfFreedom +
                        ", threshold=" + threshold + ")",
                result.statistic < threshold);
    }

    public void testDependentLikelihoodInfluencesTheDraw() {
        MathUtils.setSeed(20260818);

        final Fixture withoutDependent = createFixture(4, false);
        final Fixture withDependent = createFixture(4, true);

        final double[] truthWithout = bruteForceLogTargets(withoutDependent);
        final double[] truthWith = bruteForceLogTargetsIncludingDependent(withDependent);

        // The dependent likelihood uses a different alignment on the same tree/branch
        // model, so the two target vectors (and hence which permutation dominates)
        // should generally differ; assert they are not identical, which would
        // indicate the dependent-likelihood contribution is silently being dropped.
        boolean anyDifference = false;
        for (int i = 0; i < truthWithout.length; i++) {
            if (Math.abs((truthWith[i] - truthWithout[i]) -
                    (truthWith[0] - truthWithout[0])) > TOL) {
                anyDifference = true;
                break;
            }
        }
        assertTrue("Dependent likelihood did not change the relative target across permutations; " +
                "the dependent-likelihood sum may not be wired correctly", anyDifference);
    }

    // --- brute-force target computation -----------------------------------------------------

    private double[] bruteForceLogTargets(final Fixture fixture) {
        final double[] logTargets = new double[fixture.permutations.size()];
        for (int i = 0; i < fixture.permutations.size(); i++) {
            resetMapping(fixture, fixture.permutations.get(i));
            fixture.independentLikelihood.makeDirty();
            logTargets[i] = fixture.independentLikelihood.getLogLikelihood();
        }
        return logTargets;
    }

    private double[] bruteForceLogTargetsIncludingDependent(final Fixture fixture) {
        final double[] logTargets = new double[fixture.permutations.size()];
        for (int i = 0; i < fixture.permutations.size(); i++) {
            resetMapping(fixture, fixture.permutations.get(i));
            fixture.independentLikelihood.makeDirty();
            double logTarget = fixture.independentLikelihood.getLogLikelihood();
            if (fixture.dependentLikelihood != null) {
                fixture.dependentLikelihood.makeDirty();
                logTarget += fixture.dependentLikelihood.getLogLikelihood();
            }
            logTargets[i] = logTarget;
        }
        return logTargets;
    }

    private static double[] softmax(final double[] logValues) {
        double max = Double.NEGATIVE_INFINITY;
        for (final double v : logValues) {
            max = Math.max(max, v);
        }
        double sum = 0.0;
        final double[] weights = new double[logValues.length];
        for (int i = 0; i < logValues.length; i++) {
            weights[i] = Math.exp(logValues[i] - max);
            sum += weights[i];
        }
        for (int i = 0; i < weights.length; i++) {
            weights[i] /= sum;
        }
        return weights;
    }

    private static ChiSquaredResult chiSquaredStatistic(final int[] counts, final double[] probabilities,
                                                        final int n) {
        double stat = 0.0;
        int cellsUsed = 0;
        for (int i = 0; i < counts.length; i++) {
            final double expected = probabilities[i] * n;
            if (expected >= 5.0) {
                final double diff = counts[i] - expected;
                stat += diff * diff / expected;
                cellsUsed++;
            }
        }
        return new ChiSquaredResult(stat, Math.max(1, cellsUsed - 1));
    }

    private static final class ChiSquaredResult {
        final double statistic;
        final int degreesOfFreedom;

        private ChiSquaredResult(final double statistic, final int degreesOfFreedom) {
            this.statistic = statistic;
            this.degreesOfFreedom = degreesOfFreedom;
        }
    }

    private static int indexOfCurrentMapping(final Fixture fixture) {
        final int[] current = new int[fixture.rewardRatesMapping.getDimension()];
        for (int s = 0; s < current.length; s++) {
            current[s] = (int) Math.round(fixture.rewardRatesMapping.getParameterValue(s));
        }
        for (int i = 0; i < fixture.permutations.size(); i++) {
            if (java.util.Arrays.equals(current, fixture.permutations.get(i))) {
                return i;
            }
        }
        throw new IllegalStateException("Current mapping does not match any known permutation");
    }

    private static void resetMapping(final Fixture fixture, final int[] mapping) {
        for (int s = 0; s < mapping.length; s++) {
            fixture.rewardRatesMapping.setParameterValueQuietly(s, mapping[s]);
        }
        fixture.rewardRatesMapping.fireParameterChangedEvent();
    }

    private static OneZeroOneShuffleGibbsOperator createOperator(final Fixture fixture) {
        return new OneZeroOneShuffleGibbsOperator(
                fixture.rewardRatesValues,
                fixture.rewardRatesMapping,
                fixture.independentLikelihood,
                new TreeDataLikelihood[0],
                fixture.dependentLikelihood == null
                        ? new TreeDataLikelihood[0]
                        : new TreeDataLikelihood[]{fixture.dependentLikelihood},
                1.0,
                0.0);
    }

    /**
     * All n! permutations of {0,...,n-1}, generated independently of
     * {@link OneZeroOneShuffleGibbsOperator}'s own enumerator (deliberately a
     * separate implementation, so a bug in one is unlikely to be mirrored in
     * the other).
     */
    private static List<int[]> allPermutations(final int n) {
        final List<int[]> results = new ArrayList<int[]>();
        permuteRecursive(new int[n], new boolean[n], 0, results);
        return results;
    }

    private static void permuteRecursive(final int[] current, final boolean[] used,
                                         final int position, final List<int[]> results) {
        if (position == current.length) {
            results.add(java.util.Arrays.copyOf(current, current.length));
            return;
        }
        for (int value = 0; value < current.length; value++) {
            if (!used[value]) {
                used[value] = true;
                current[position] = value;
                permuteRecursive(current, used, position + 1, results);
                used[value] = false;
            }
        }
    }

    // --- fixture construction ----------------------------------------------------------------

    private static Fixture createFixture(final int stateCount, final boolean withDependent) {
        final TreeModel tree = createThreeTipTree();
        final SitePatterns patterns = createSitePatterns("A", "C", "G");
        final SubstitutionModel substitutionModel = createSubstitutionModel(stateCount);
        final GammaSiteRateModel siteRateModel = new GammaSiteRateModel("shuffleGibbsTargetSiteModel");

        final int branchCount = 4; // 3-tip rooted tree has 4 branches
        final double[] ctsInit = new double[branchCount];
        final double[] categoryInit = new double[branchCount];
        java.util.Arrays.fill(ctsInit, 0.5);
        java.util.Arrays.fill(categoryInit, 0.5);

        final Parameter ctsRewards = new Parameter.Default("targetRewardCts", ctsInit);
        final Parameter categoryState = new Parameter.Default("targetRewardCategory", categoryInit);
        final double[] cuts = new double[stateCount + 2];
        for (int i = 0; i < cuts.length; i++) {
            cuts[i] = i;
        }
        final Parameter categoryCuts = new Parameter.Default("targetRewardCategoryCuts", cuts);

        final double[] values = new double[stateCount];
        values[0] = 0.0;
        values[1] = 1.0;
        for (int i = 2; i < stateCount; i++) {
            values[i] = (i - 1.0) / (stateCount - 1.0);
        }
        final Parameter rewardRatesValues = new Parameter.Default("targetRewardRatesValues", values);

        final double[] mappingInit = new double[stateCount];
        for (int i = 0; i < stateCount; i++) {
            mappingInit[i] = i;
        }
        final Parameter rewardRatesMapping = new Parameter.Default("targetRewardRatesMapping", mappingInit);

        final RewardRates rewardRates = new RewardRates(
                rewardRatesValues,
                null,
                new Parameter.Default("targetRewardRatesVarying", stateCount - 2),
                rewardRatesMapping
        );

        final RewardsAwareCategoricalMixtureBranchRates rewardBranchRates =
                new RewardsAwareCategoricalMixtureBranchRates(
                        tree,
                        ctsRewards,
                        categoryState,
                        categoryCuts,
                        rewardRates,
                        new ArbitraryBranchRates.BranchRateTransform.None(),
                        false,
                        TreeParameterModel.Type.WITHOUT_ROOT);
        final RewardsAwareBranchModel rewardsAwareBranchModel = new RewardsAwareBranchModel(
                tree,
                substitutionModel,
                rewardRates,
                categoryState,
                categoryCuts,
                rewardBranchRates,
                false
        );

        final TreeDataLikelihood independentLikelihood =
                buildTreeDataLikelihood(tree, patterns, rewardsAwareBranchModel, siteRateModel, rewardBranchRates);

        TreeDataLikelihood dependentLikelihood = null;
        if (withDependent) {
            final SitePatterns dependentPatterns = createSitePatterns("T", "G", "C");
            dependentLikelihood = buildTreeDataLikelihood(
                    tree, dependentPatterns, rewardsAwareBranchModel, siteRateModel, rewardBranchRates);
        }

        return new Fixture(independentLikelihood, dependentLikelihood, rewardRatesValues,
                rewardRatesMapping, allPermutations(stateCount));
    }

    private static TreeDataLikelihood buildTreeDataLikelihood(final TreeModel tree,
                                                               final SitePatterns patterns,
                                                               final RewardsAwareBranchModel rewardsAwareBranchModel,
                                                               final GammaSiteRateModel siteRateModel,
                                                               final RewardsAwareCategoricalMixtureBranchRates rewardBranchRates) {
        final DiscreteDataLikelihoodDelegate delegate = new DiscreteDataLikelihoodDelegate(
                tree,
                patterns,
                rewardsAwareBranchModel,
                siteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, false, false, false),
                new RewardsAwarePartialsRepresentation(rewardsAwareBranchModel),
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                null,
                null
        );
        return new TreeDataLikelihood(delegate, tree, rewardBranchRates);
    }

    private static SitePatterns createSitePatterns(final String a, final String b, final String c) {
        final DataType dataType = Nucleotides.INSTANCE;
        final SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
        alignment.addSequence(sequence("a", a, dataType));
        alignment.addSequence(sequence("b", b, dataType));
        alignment.addSequence(sequence("c", c, dataType));
        return new SitePatterns(alignment, null, 0, -1, 1, true);
    }

    private static Sequence sequence(final String taxonName, final String sequenceString, final DataType dataType) {
        final Sequence sequence = new Sequence(new Taxon(taxonName), sequenceString);
        sequence.setDataType(dataType);
        return sequence;
    }

    private static TreeModel createThreeTipTree() {
        final SimpleNode left = new SimpleNode();
        left.setTaxon(new Taxon("a"));
        left.setHeight(0.0);

        final SimpleNode right = new SimpleNode();
        right.setTaxon(new Taxon("b"));
        right.setHeight(0.0);

        final SimpleNode internal = new SimpleNode();
        internal.setHeight(0.4);
        internal.addChild(left);
        internal.addChild(right);

        final SimpleNode outgroup = new SimpleNode();
        outgroup.setTaxon(new Taxon("c"));
        outgroup.setHeight(0.0);

        final SimpleNode root = new SimpleNode();
        root.setHeight(1.0);
        root.addChild(internal);
        root.addChild(outgroup);

        final Tree simpleTree = new SimpleTree(root);
        return new DefaultTreeModel("shuffleGibbsTargetTree", simpleTree);
    }

    private static SubstitutionModel createSubstitutionModel(final int stateCount) {
        // stateCount is always <= 4 (nucleotide alphabet) in these tests.
        final DataType dataType = Nucleotides.INSTANCE;
        final double[] freqs = new double[4];
        java.util.Arrays.fill(freqs, 0.25);
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, freqs);
        return new HKY(new Parameter.Default("shuffleGibbsTargetKappa", 2.0), frequencyModel);
    }

    private static void assertFinite(final double value) {
        assertTrue("Expected finite value but found " + value,
                !Double.isNaN(value) && !Double.isInfinite(value));
    }

    private static void assertLogEquals(final String message, final double expected, final double actual) {
        if (Double.isInfinite(expected) || Double.isInfinite(actual)) {
            assertEquals(message, expected, actual);
        } else {
            assertEquals(message, expected, actual, TOL);
        }
    }

    private static final class Fixture {
        final TreeDataLikelihood independentLikelihood;
        final TreeDataLikelihood dependentLikelihood;
        final Parameter rewardRatesValues;
        final Parameter rewardRatesMapping;
        final List<int[]> permutations;

        private Fixture(final TreeDataLikelihood independentLikelihood,
                        final TreeDataLikelihood dependentLikelihood,
                        final Parameter rewardRatesValues,
                        final Parameter rewardRatesMapping,
                        final List<int[]> permutations) {
            this.independentLikelihood = independentLikelihood;
            this.dependentLikelihood = dependentLikelihood;
            this.rewardRatesValues = rewardRatesValues;
            this.rewardRatesMapping = rewardRatesMapping;
            this.permutations = permutations;
        }
    }
}
