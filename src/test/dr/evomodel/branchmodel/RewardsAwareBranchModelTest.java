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
                new double[]{1.0, 0.0}
        );

        double[] matrix = fixture.branchModel.getTransitionMatrix(fixture.tree.getNode(0));
        double expected = Math.exp(-fixture.branchModel.getUniformizationRate() *
                fixture.tree.getBranchLength(fixture.tree.getNode(0)));

        assertEquals(0.0, matrix[0], TOL);
        assertEquals(0.0, matrix[1], TOL);
        assertEquals(0.0, matrix[2], TOL);
        assertEquals(expected, matrix[3], TOL);
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
        TreeModel tree = createTwoTipTree();
        SubstitutionModel substitutionModel = createTwoStateSubstitutionModel();

        Parameter rewardProportion = new Parameter.Default("rewardProportion", rewardProportionValues);
        ArbitraryBranchRates branchRates = new ArbitraryBranchRates(
                tree,
                rewardProportion,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false
        );

        Parameter rewardRateValues = new Parameter.Default("rewardRates", new double[]{0.0, 1.0});
        Parameter rewardRateInternalValues = new Parameter.Default("rewardRatesInternal", new double[0]);
        Parameter rewardRateMapping = new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0});
        Parameter indicator = new Parameter.Default("indicator", indicatorValues);
        Parameter atomIndices = new Parameter.Default("atomIndices", atomValues);

        RewardsAwareBranchModel branchModel = new RewardsAwareBranchModel(
                tree,
                substitutionModel,
                rewardRateValues,
                rewardRateInternalValues,
                rewardRateMapping,
                indicator,
                branchRates,
                atomIndices,
                false
        );

        return new Fixture(tree, branchModel);
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

    private static SubstitutionModel createTwoStateSubstitutionModel() {
        DataType dataType = TwoStates.INSTANCE;
        FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.5, 0.5});
        Parameter rates = new Parameter.Default("twoStateRates", new double[]{2.0, 3.0});
        ComplexSubstitutionModel substitutionModel =
                new ComplexSubstitutionModel("twoState", dataType, frequencyModel, rates);
        substitutionModel.setNormalization(false);
        return substitutionModel;
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

        Fixture(TreeModel tree, RewardsAwareBranchModel branchModel) {
            this.tree = tree;
            this.branchModel = branchModel;
        }
    }
}
