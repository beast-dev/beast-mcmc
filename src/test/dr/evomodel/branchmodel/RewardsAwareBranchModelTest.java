/*
 * RewardsAwareBranchModelTest.java
 *
 * Copyright (c) 2002-2024 the BEAST Development Team
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

package test.dr.evomodel.branchmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStates;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.markovjumps.SericolaSeriesMarkovRewardFastModel;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.RewardsMixtureIndicatorAndAtomIndicesOperatorParser;
import dr.xml.XMLSyntaxRule;
import test.dr.math.MathTestCase;

/**
 * Focused smoke tests for the reward-aware transition-matrix path.
 *
 * @author Filippo Monti
 */
public class RewardsAwareBranchModelTest extends MathTestCase {

    private static final double TOL = 1.0e-10;

    public void testGibbsOperatorParserDeclaresTreeDataLikelihood() {
        XMLSyntaxRule[] rules = new RewardsMixtureIndicatorAndAtomIndicesOperatorParser().getSyntaxRules();

        assertTrue(acceptsElementClass(rules, TreeDataLikelihood.class));
        assertTrue(acceptsElementName(rules, "indicatorZ"));
        assertTrue(acceptsElementName(rules, "atomIndex"));
        assertTrue(acceptsElementName(rules, "dependentCtmcLikelihoods"));
    }

    public void testContinuousTransitionMatrixIsFiniteAndNonNegative() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );

        double[] matrix = fixture.branchModel.getTransitionMatrix(fixture.tree.getNode(0));

        assertEquals(4, matrix.length);
        assertFiniteNonNegative(matrix);
        assertTrue(accumulate(matrix) > 0.0);
    }

    public void testAtomicTransitionMatrixHasSingleNoJumpMass() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{1.0, 0.0},
                new double[]{0.0, 1.0}
        );

        final int atomState = 0;
        final double branchLength = fixture.tree.getBranchLength(fixture.tree.getNode(0));
        final double[] infinitesimalMatrix = infinitesimalMatrix(fixture.substitutionModel);
        final double expectedLogScale =
                infinitesimalMatrix[atomState * fixture.branchModel.getStateCount() + atomState] * branchLength;
        final double expected = Math.exp(expectedLogScale);
        final double zeroUniformizationEventMass =
                Math.exp(-fixture.branchModel.getUniformizationRate() * branchLength);

        assertTrue(Math.abs(expected - zeroUniformizationEventMass) > 1.0e-8);
        assertEquals(expectedLogScale,
                fixture.branchModel.getAtomicBranchLogScaleForState(fixture.tree.getNode(0).getNumber(), atomState),
                TOL);

        double[] matrix = fixture.branchModel.getTransitionMatrix(fixture.tree.getNode(0));

        assertEquals(expected, matrix[0], TOL);
        assertEquals(0.0, matrix[1], TOL);
        assertEquals(0.0, matrix[2], TOL);
        assertEquals(0.0, matrix[3], TOL);
    }

    public void testSericolaDerivativeMatchesFiniteDifference() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double rewardProportion = 0.45;
        double time = 0.9;
        double h = 1.0e-6;

        double[] plus = new double[4];
        double[] minus = new double[4];
        double[] differential = new double[4];

        sericola.computePdfInto(rewardProportion + h, time, plus);
        sericola.computePdfInto(rewardProportion - h, time, minus);
        sericola.computePdfDerivativeWrtRewardProportionInto(rewardProportion, time, differential, false);

        for (int i = 0; i < differential.length; i++) {
            double finiteDifference = (plus[i] - minus[i]) / (2.0 * h);
            double tolerance = Math.max(1.0e-6, Math.abs(finiteDifference) * 1.0e-5);
            assertEquals("entry " + i, finiteDifference, differential[i], tolerance);
        }
    }

    public void testSericolaInfinitesimalMatrixAdjointMatchesFiniteDifferenceThroughRates() {
        Fixture fixture = createFixtureWithSubstitutionRates(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double rewardProportion = 0.45;
        double time = 0.9;
        double[] upstream = new double[]{0.70, -0.25, 0.40, 1.10};
        double[] qAdjoint = new double[4];

        sericola.computePdfGradientWrtInfinitesimalMatrixInto(
                rewardProportion,
                time,
                upstream,
                qAdjoint);

        for (int p = 0; p < fixture.substitutionRates.getDimension(); p++) {
            double oldValue = fixture.substitutionRates.getParameterValue(p);
            double step = 1.0e-6 * Math.max(1.0, Math.abs(oldValue));

            fixture.substitutionRates.setParameterValue(p, oldValue + step);
            double objectivePlus = sericolaObjective(sericola, rewardProportion, time, upstream);
            double[] qPlus = infinitesimalMatrix(fixture.substitutionModel);

            fixture.substitutionRates.setParameterValue(p, oldValue - step);
            double objectiveMinus = sericolaObjective(sericola, rewardProportion, time, upstream);
            double[] qMinus = infinitesimalMatrix(fixture.substitutionModel);

            fixture.substitutionRates.setParameterValue(p, oldValue);

            double finiteDifference = (objectivePlus - objectiveMinus) / (2.0 * step);
            double chainRule = 0.0;
            for (int i = 0; i < qAdjoint.length; i++) {
                chainRule += qAdjoint[i] * (qPlus[i] - qMinus[i]) / (2.0 * step);
            }

            double tolerance = Math.max(1.0e-6, Math.abs(finiteDifference) * 2.0e-5);
            assertEquals("rate parameter " + p, finiteDifference, chainRule, tolerance);
        }
    }

    public void testSericolaRewardRateAdjointMatchesFiniteDifferenceThroughValues() {
        Fixture fixture = createFixture(
                createTwoStateSubstitutionFixture(),
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0},
                new double[]{0.2, 0.8}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double rewardProportion = 0.45;
        double time = 0.9;
        double[] upstream = new double[]{0.70, -0.25, 0.40, 1.10};
        double[] rewardRateAdjoint = new double[fixture.rewardRateValues.getDimension()];

        sericola.computePdfGradientWrtRewardRatesInto(
                rewardProportion,
                time,
                upstream,
                rewardRateAdjoint);

        for (int p = 0; p < fixture.rewardRateValues.getDimension(); p++) {
            double oldValue = fixture.rewardRateValues.getParameterValue(p);
            double step = 1.0e-6;

            fixture.rewardRateValues.setParameterValue(p, oldValue + step);
            double objectivePlus = sericolaObjective(sericola, rewardProportion, time, upstream);

            fixture.rewardRateValues.setParameterValue(p, oldValue - step);
            double objectiveMinus = sericolaObjective(sericola, rewardProportion, time, upstream);

            fixture.rewardRateValues.setParameterValue(p, oldValue);

            double finiteDifference = (objectivePlus - objectiveMinus) / (2.0 * step);
            double tolerance = Math.max(2.0e-6, Math.abs(finiteDifference) * 5.0e-5);
            assertEquals("reward-rate value " + p, finiteDifference, rewardRateAdjoint[p], tolerance);
        }
    }

    public void testVectorizedPdfWithSharedTimeMatchesScalarCalls() {
        Fixture fixture = createFixture(
                new double[]{0.25, 0.75},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double[] rewardProportions = new double[]{0.2, 0.45, 0.8};
        double time = 0.9;
        double[][] vectorized = newMatrixRows(rewardProportions.length, 4);

        sericola.computePdfInto(rewardProportions, time, vectorized);

        assertVectorizedPdfMatchesScalarCalls(sericola, rewardProportions, new double[]{time}, vectorized);
    }

    public void testVectorizedPdfWithPerEntryTimesMatchesScalarCalls() {
        Fixture fixture = createFixture(
                new double[]{0.25, 0.75},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double[] rewardProportions = new double[]{0.2, 0.45, 0.8};
        double[] times = new double[]{0.4, 0.9, 1.3};
        double[][] vectorized = newMatrixRows(rewardProportions.length, 4);

        sericola.computePdfInto(rewardProportions, times, vectorized);

        assertVectorizedPdfMatchesScalarCalls(sericola, rewardProportions, times, vectorized);
    }

    public void testVectorizedDerivativeWithSharedTimeMatchesScalarCalls() {
        Fixture fixture = createFixture(
                new double[]{0.25, 0.75},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double[] rewardProportions = new double[]{0.2, 0.45, 0.8};
        double time = 0.9;
        double[][] vectorized = newMatrixRows(rewardProportions.length, 4);

        sericola.computePdfDerivativeWrtRewardProportionInto(rewardProportions, new double[]{time}, vectorized, false);

        assertVectorizedDerivativeMatchesScalarCalls(sericola, rewardProportions, new double[]{time}, vectorized);
    }

    public void testVectorizedDerivativeWithPerEntryTimesMatchesScalarCalls() {
        Fixture fixture = createFixture(
                new double[]{0.25, 0.75},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double[] rewardProportions = new double[]{0.2, 0.45, 0.8};
        double[] times = new double[]{0.4, 0.9, 1.3};
        double[][] vectorized = newMatrixRows(rewardProportions.length, 4);

        sericola.computePdfDerivativeWrtRewardProportionInto(rewardProportions, times, vectorized, false);

        assertVectorizedDerivativeMatchesScalarCalls(sericola, rewardProportions, times, vectorized);
    }

    public void testDerivativeContractionMatchesMaterializedDerivative() {
        Fixture fixture = createFixture(
                new double[]{0.25, 0.75},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double rewardProportion = 0.45;
        double time = 0.9;
        double[] pre = new double[]{0.70, 0.30};
        double[] post = new double[]{0.15, 0.85};
        double[] materialized = new double[4];

        sericola.computePdfDerivativeWrtRewardProportionInto(rewardProportion, time, materialized, false);

        double expected = bilinearFormStable(pre, materialized, post, 2);
        double observed = sericola.contractPdfDerivativeWrtRewardProportion(
                rewardProportion,
                time,
                pre,
                post,
                false);

        assertEquals(expected, observed, Math.max(1.0e-10, Math.abs(expected) * 1.0e-10));
    }

    public void testVectorizedDerivativeContractionMatchesMaterializedDerivative() {
        Fixture fixture = createFixture(
                new double[]{0.25, 0.75},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        double[] rewardProportions = new double[]{0.2, 0.45, 0.8};
        double[] times = new double[]{0.4, 0.9, 1.3};
        double[][] pre = new double[][]{
                {0.80, 0.20},
                {0.55, 0.45},
                {0.10, 0.90}
        };
        double[][] post = new double[][]{
                {0.30, 0.70},
                {0.65, 0.35},
                {0.40, 0.60}
        };
        double[] contracted = new double[rewardProportions.length];

        sericola.contractPdfDerivativeWrtRewardProportionInto(
                rewardProportions,
                times,
                pre,
                post,
                contracted,
                false);

        double[] materialized = new double[4];
        for (int t = 0; t < rewardProportions.length; t++) {
            sericola.computePdfDerivativeWrtRewardProportionInto(
                    rewardProportions[t],
                    times[t],
                    materialized,
                    false);
            double expected = bilinearFormStable(pre[t], materialized, post[t], 2);
            assertEquals("row " + t, expected, contracted[t],
                    Math.max(1.0e-10, Math.abs(expected) * 1.0e-10));
        }
    }

    public void testVectorizedDerivativeRejectsBoundaryRewardProportion() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        try {
            sericola.computePdfDerivativeWrtRewardProportionInto(
                    new double[]{0.25, 0.0, 0.75},
                    new double[]{0.9, 0.9, 0.9},
                    newMatrixRows(3, 4),
                    false);
            fail("Expected UnsupportedOperationException for a boundary rewardProportion in a vectorized derivative call");
        } catch (UnsupportedOperationException e) {
            assertTrue(e.getMessage().contains("Boundary Values for rewardProportion touched"));
            assertTrue(e.getMessage().contains("rewardProportion=0.0"));
        }
    }

    public void testVectorizedDerivativeContractionRejectsBoundaryRewardProportion() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        try {
            sericola.contractPdfDerivativeWrtRewardProportionInto(
                    new double[]{0.25, 0.0, 0.75},
                    new double[]{0.9, 0.9, 0.9},
                    new double[][]{{0.7, 0.3}, {0.7, 0.3}, {0.7, 0.3}},
                    new double[][]{{0.2, 0.8}, {0.2, 0.8}, {0.2, 0.8}},
                    new double[3],
                    false);
            fail("Expected UnsupportedOperationException for a boundary rewardProportion in a vectorized derivative-contraction call");
        } catch (UnsupportedOperationException e) {
            assertTrue(e.getMessage().contains("Boundary Values for rewardProportion touched"));
            assertTrue(e.getMessage().contains("rewardProportion=0.0"));
        }
    }

    public void testRewardProportionOutsideSupportThrowsDiagnosticException() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        assertIllegalRewardProportion(sericola, -0.1, "< min(alpha)");
        assertIllegalRewardProportion(sericola, 1.1, "> max(alpha)");
    }

    public void testSericolaDerivativeRejectsBoundaryRewardProportion() {
        Fixture fixture = createFixture(
                new double[]{0.35, 0.65},
                new double[]{0.0, 0.0},
                new double[]{0.0, 1.0}
        );
        SericolaSeriesMarkovRewardFastModel sericola = fixture.branchModel.getSericolaModel();

        assertBoundaryDerivativeRejected(sericola, 0.0);
        assertBoundaryDerivativeRejected(sericola, 1.0);
    }

    private static Fixture createFixture(double[] rewardProportionValues, double[] indicatorValues, double[] atomValues) {
        return createFixture(createTwoStateSubstitutionFixture(), rewardProportionValues, indicatorValues, atomValues);
    }

    private static Fixture createFixtureWithSubstitutionRates(
            double[] rewardProportionValues,
            double[] indicatorValues,
            double[] atomValues) {

        return createFixture(createTwoStateSubstitutionFixture(), rewardProportionValues, indicatorValues, atomValues);
    }

    private static Fixture createFixture(
            SubstitutionFixture substitutionFixture,
            double[] rewardProportionValues,
            double[] indicatorValues,
            double[] atomValues) {

        return createFixture(substitutionFixture, rewardProportionValues, indicatorValues, atomValues,
                new double[]{0.0, 1.0});
    }

    private static Fixture createFixture(
            SubstitutionFixture substitutionFixture,
            double[] rewardProportionValues,
            double[] indicatorValues,
            double[] atomValues,
            double[] rewardRateValueArray) {

        TreeModel tree = createTwoTipTree();
        SubstitutionModel substitutionModel = substitutionFixture.substitutionModel;

        Parameter rewardProportion = new Parameter.Default("rewardProportion", rewardProportionValues);
        ArbitraryBranchRates branchRates = new ArbitraryBranchRates(
                tree,
                rewardProportion,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false
        );

        Parameter rewardRateValues = new Parameter.Default("rewardRates", rewardRateValueArray);
        Parameter rewardRateInternalValues = new Parameter.Default("rewardRatesInternal", new double[0]);
        Parameter rewardRateMapping = new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0});
        RewardRates rewardRates = new RewardRates(rewardRateValues, null, rewardRateInternalValues, rewardRateMapping);
        Parameter indicator = new Parameter.Default("indicator", indicatorValues);
        Parameter atomIndices = new Parameter.Default("atomIndices", atomValues);

        RewardsAwareBranchModel branchModel = new RewardsAwareBranchModel(
                tree,
                substitutionModel,
                rewardRates,
                indicator,
                branchRates,
                atomIndices,
                false
        );

        return new Fixture(tree, branchModel, substitutionModel, substitutionFixture.rates, rewardRateValues);
    }

    private static double sericolaObjective(
            SericolaSeriesMarkovRewardFastModel sericola,
            double rewardProportion,
            double time,
            double[] upstream) {

        double[] density = new double[upstream.length];
        sericola.computePdfInto(rewardProportion, time, density);
        double objective = 0.0;
        for (int i = 0; i < density.length; i++) {
            objective += upstream[i] * density[i];
        }
        return objective;
    }

    private static double[] infinitesimalMatrix(SubstitutionModel substitutionModel) {
        int n = substitutionModel.getDataType().getStateCount();
        double[] matrix = new double[n * n];
        substitutionModel.getInfinitesimalMatrix(matrix);
        return matrix;
    }

    private static void assertVectorizedPdfMatchesScalarCalls(
            SericolaSeriesMarkovRewardFastModel sericola,
            double[] rewardProportions,
            double[] times,
            double[][] vectorized) {

        final boolean singleTime = times.length == 1;
        double[] scalar = new double[4];

        for (int t = 0; t < rewardProportions.length; t++) {
            double time = singleTime ? times[0] : times[t];
            sericola.computePdfInto(rewardProportions[t], time, scalar);
            assertMatrixEntryEquals("row " + t, scalar, vectorized[t]);
        }
    }

    private static void assertVectorizedDerivativeMatchesScalarCalls(
            SericolaSeriesMarkovRewardFastModel sericola,
            double[] rewardProportions,
            double[] times,
            double[][] vectorized) {

        final boolean singleTime = times.length == 1;
        double[] scalar = new double[4];

        for (int t = 0; t < rewardProportions.length; t++) {
            double time = singleTime ? times[0] : times[t];
            sericola.computePdfDerivativeWrtRewardProportionInto(rewardProportions[t], time, scalar, false);
            assertMatrixEntryEquals("row " + t, scalar, vectorized[t]);
        }
    }

    private static double[][] newMatrixRows(int rows, int columns) {
        double[][] matrix = new double[rows][];
        for (int i = 0; i < rows; i++) {
            matrix[i] = new double[columns];
        }
        return matrix;
    }

    private static void assertMatrixEntryEquals(String label, double[] expected, double[] observed) {
        assertEquals(expected.length, observed.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(label + " entry " + i, expected[i], observed[i], TOL);
        }
    }

    private static double bilinearFormStable(double[] pre, double[] D, double[] post, int n) {
        double acc = 0.0;
        double cAcc = 0.0;

        for (int i = 0; i < n; i++) {
            double pre_i = pre[i];
            if (pre_i == 0.0) continue;

            int rowBase = i * n;
            double rowDot = 0.0;
            double cRow = 0.0;

            for (int j = 0; j < n; j++) {
                double y = D[rowBase + j] * post[j] - cRow;
                double t = rowDot + y;
                cRow = (t - rowDot) - y;
                rowDot = t;
            }

            double y = pre_i * rowDot - cAcc;
            double t = acc + y;
            cAcc = (t - acc) - y;
            acc = t;
        }
        return acc;
    }

    private static void assertIllegalRewardProportion(
            SericolaSeriesMarkovRewardFastModel sericola,
            double rewardProportion,
            String expectedMessageFragment) {

        try {
            sericola.computePdfInto(rewardProportion, 0.9, new double[4]);
            fail("Expected IllegalArgumentException for rewardProportion=" + rewardProportion);
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage().contains("rewardProportion=" + rewardProportion));
            assertTrue(e.getMessage().contains(expectedMessageFragment));
            assertTrue(e.getMessage().contains("snapping to boundary is not supported"));
        }
    }

    private static void assertBoundaryDerivativeRejected(
            SericolaSeriesMarkovRewardFastModel sericola,
            double rewardProportion) {

        try {
            sericola.computePdfDerivativeWrtRewardProportionInto(rewardProportion, 0.9, new double[4], false);
            fail("Expected UnsupportedOperationException for boundary rewardProportion=" + rewardProportion);
        } catch (UnsupportedOperationException e) {
            assertTrue(e.getMessage().contains("Boundary Values for rewardProportion touched"));
            assertTrue(e.getMessage().contains("one-sided interior derivative is implemented but disabled"));
            assertTrue(e.getMessage().contains("rewardProportion=" + rewardProportion));
        }
    }

    private static TreeModel createTwoTipTree() {
        SimpleNode left = new SimpleNode();
        left.setTaxon(new Taxon("a"));
        left.setHeight(0.0);

        SimpleNode right = new SimpleNode();
        right.setTaxon(new Taxon("b"));
        right.setHeight(0.0);

        SimpleNode root = new SimpleNode();
        root.setHeight(1.0);
        root.addChild(left);
        root.addChild(right);

        Tree tree = new SimpleTree(root);
        return new DefaultTreeModel("rewardAwareTestTree", tree);
    }

    private static SubstitutionFixture createTwoStateSubstitutionFixture() {
        DataType dataType = TwoStates.INSTANCE;
        FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.5, 0.5});
        Parameter rates = new Parameter.Default("twoStateRates", new double[]{2.0, 3.0});
        ComplexSubstitutionModel substitutionModel =
                new ComplexSubstitutionModel("twoState", dataType, frequencyModel, rates);
        substitutionModel.setNormalization(false);
        return new SubstitutionFixture(substitutionModel, rates);
    }

    private static void assertFiniteNonNegative(double[] values) {
        for (int i = 0; i < values.length; i++) {
            assertTrue("entry " + i + " is not finite: " + values[i],
                    !Double.isNaN(values[i]) && !Double.isInfinite(values[i]));
            assertTrue("entry " + i + " is negative: " + values[i], values[i] >= -TOL);
        }
    }

    private static boolean acceptsElementClass(XMLSyntaxRule[] rules, Class c) {
        for (XMLSyntaxRule rule : rules) {
            if (rule.isLegalElementClass(c)) {
                return true;
            }
        }
        return false;
    }

    private static boolean acceptsElementName(XMLSyntaxRule[] rules, String name) {
        for (XMLSyntaxRule rule : rules) {
            if (rule.isLegalElementName(name)) {
                return true;
            }
        }
        return false;
    }

    private static final class Fixture {
        final TreeModel tree;
        final RewardsAwareBranchModel branchModel;
        final SubstitutionModel substitutionModel;
        final Parameter substitutionRates;
        final Parameter rewardRateValues;

        Fixture(
                TreeModel tree,
                RewardsAwareBranchModel branchModel,
                SubstitutionModel substitutionModel,
                Parameter substitutionRates,
                Parameter rewardRateValues) {
            this.tree = tree;
            this.branchModel = branchModel;
            this.substitutionModel = substitutionModel;
            this.substitutionRates = substitutionRates;
            this.rewardRateValues = rewardRateValues;
        }
    }

    private static final class SubstitutionFixture {
        final SubstitutionModel substitutionModel;
        final Parameter rates;

        SubstitutionFixture(SubstitutionModel substitutionModel, Parameter rates) {
            this.substitutionModel = substitutionModel;
            this.rates = rates;
        }
    }
}
