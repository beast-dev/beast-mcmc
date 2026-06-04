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

        double rho = 0.45;
        double time = 0.9;
        double h = 1.0e-6;

        double[] plus = new double[4];
        double[] minus = new double[4];
        double[] differential = new double[4];

        sericola.computePdfInto(rho + h, time, plus);
        sericola.computePdfInto(rho - h, time, minus);
        sericola.computePdfDerivativeWrtRhoInto(rho, time, differential, false);

        for (int i = 0; i < differential.length; i++) {
            double finiteDifference = (plus[i] - minus[i]) / (2.0 * h);
            double tolerance = Math.max(1.0e-6, Math.abs(finiteDifference) * 1.0e-5);
            assertEquals("entry " + i, finiteDifference, differential[i], tolerance);
        }
    }

    private static Fixture createFixture(double[] rhoValues, double[] indicatorValues, double[] atomValues) {
        TreeModel tree = createTwoTipTree();
        SubstitutionModel substitutionModel = createTwoStateSubstitutionModel();

        Parameter rho = new Parameter.Default("rho", rhoValues);
        ArbitraryBranchRates branchRates = new ArbitraryBranchRates(
                tree,
                rho,
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
