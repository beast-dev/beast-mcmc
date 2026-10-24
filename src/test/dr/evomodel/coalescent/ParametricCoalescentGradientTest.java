/*
 * ParametricCoalescentGradientTest.java
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

package test.dr.evomodel.coalescent;

import dr.evolution.coalescent.CoalescentGradient;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.coalescent.demographicmodel.ConstantPopulationModel;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.demographicmodel.ExpansionModel;
import dr.evomodel.coalescent.demographicmodel.ExponentialGrowthModel;
import dr.evomodel.coalescent.demographicmodel.LogisticGrowthModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Gradients of the constant, exponential, logistic and expansion coalescents with respect to node heights and the
 * demographic parameters, against central finite differences of the coalescent log-likelihood (plan 3, G1).
 */
public class ParametricCoalescentGradientTest extends TestCase {

    // heterochronous, so sampling intervals sit between coalescent ones
    private static final String NEWICK = "((A:0.12,B:0.32):0.31,(C:0.25,(D:0.07,E:0.17):0.18):0.18);";
    // two cherries at the same height
    private static final String TIED = "((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);";

    private static final double H = 1E-6;
    private static final double TOLERANCE = 1E-6;

    private enum Intervals {BIG_FAST, MAPPED, FAST}

    public ParametricCoalescentGradientTest(String name) {
        super(name);
    }

    public void testConstant() throws Exception {
        Parameter N = new Parameter.Default(0.7);
        checkAll(new ConstantPopulationModel(N, Units.Type.YEARS), "constant", N);
    }

    public void testExponential() throws Exception {
        for (double r : new double[]{0.5, -0.3, 2.5, 1E-7, 0.0}) {
            Parameter N0 = new Parameter.Default(0.4);
            Parameter growthRate = new Parameter.Default(r);
            checkAll(new ExponentialGrowthModel(N0, growthRate, Units.Type.YEARS, true),
                    "exponential r=" + r, N0, growthRate);
        }
    }

    public void testExponentialDoublingTime() throws Exception {
        Parameter N0 = new Parameter.Default(0.4);
        Parameter doublingTime = new Parameter.Default(0.6);
        checkAll(new ExponentialGrowthModel(N0, doublingTime, Units.Type.YEARS, false),
                "exponential doubling time", N0, doublingTime);
    }

    public void testLogistic() throws Exception {
        for (double r : new double[]{2.0, 6.0}) {
            Parameter N0 = new Parameter.Default(0.5);
            Parameter growthRate = new Parameter.Default(r);
            Parameter t50 = new Parameter.Default(0.6);
            checkAll(new LogisticGrowthModel(N0, growthRate, t50, 0.5, Units.Type.YEARS, true),
                    "logistic r=" + r, N0, growthRate, t50);
        }
    }

    public void testLogisticDoublingTime() throws Exception {
        Parameter N0 = new Parameter.Default(0.5);
        Parameter doublingTime = new Parameter.Default(0.4);
        Parameter t50 = new Parameter.Default(0.6);
        checkAll(new LogisticGrowthModel(N0, doublingTime, t50, 0.5, Units.Type.YEARS, false),
                "logistic doubling time", N0, doublingTime, t50);
    }

    public void testExpansion() throws Exception {
        for (double p : new double[]{0.05, 0.4, 0.9}) {
            Parameter N0 = new Parameter.Default(0.6);
            Parameter proportion = new Parameter.Default(p);
            Parameter growthRate = new Parameter.Default(3.0);
            checkAll(new ExpansionModel(N0, proportion, growthRate, Units.Type.YEARS, true),
                    "expansion p=" + p, N0, proportion, growthRate);
        }
    }

    public void testExpansionDoublingTime() throws Exception {
        Parameter N0 = new Parameter.Default(0.6);
        Parameter proportion = new Parameter.Default(0.3);
        Parameter doublingTime = new Parameter.Default(0.25);
        checkAll(new ExpansionModel(N0, proportion, doublingTime, Units.Type.YEARS, false),
                "expansion doubling time", N0, proportion, doublingTime);
    }

    // tied nodes share the derivative of moving them together
    public void testTiedHeights() throws Exception {
        for (Intervals type : Intervals.values()) {
            DefaultTreeModel tree = tree(TIED);
            DemographicModel model = new ExponentialGrowthModel(new Parameter.Default(0.4),
                    new Parameter.Default(1.5), Units.Type.YEARS, true);
            CoalescentLikelihood likelihood = new CoalescentLikelihood(intervals(tree, type), model);
            double[] analytic = new CoalescentGradient(likelihood, tree, null,
                    CoalescentGradient.Wrt.NODE_HEIGHTS, 1E-4).getGradientLogDensity();

            NodeRef left = tree.getParent(tree.getExternalNode(0));
            NodeRef right = tree.getParent(tree.getExternalNode(2));
            assertNotSame(left, right);
            final double height = tree.getNodeHeight(left);
            tree.setNodeHeight(left, height + H);
            tree.setNodeHeight(right, height + H);
            final double plus = likelihood.getLogLikelihood();
            tree.setNodeHeight(left, height - H);
            tree.setNodeHeight(right, height - H);
            final double minus = likelihood.getLogLikelihood();
            tree.setNodeHeight(left, height);
            tree.setNodeHeight(right, height);

            final int external = tree.getExternalNodeCount();
            final double leftValue = analytic[left.getNumber() - external];
            final double rightValue = analytic[right.getNumber() - external];
            assertEquals(type + " split equally", leftValue, rightValue, 0.0);
            assertEquals(type.toString(), (plus - minus) / (2 * H), leftValue + rightValue, TOLERANCE);
        }
    }

    public void testRefusals() throws Exception {
        DefaultTreeModel tree = tree(NEWICK);
        Parameter N0 = new Parameter.Default(0.4);
        DemographicModel model = new ExponentialGrowthModel(N0, new Parameter.Default(0.5), Units.Type.YEARS, true);
        CoalescentLikelihood likelihood = new CoalescentLikelihood(new BigFastTreeIntervals(tree), model);
        try {
            new CoalescentGradient(likelihood, tree, new Parameter.Default(1.0), CoalescentGradient.Wrt.PARAMETER, 1E-4);
            fail("an unrelated parameter must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        CoalescentGradient gradient = new CoalescentGradient(likelihood, tree, N0, CoalescentGradient.Wrt.PARAMETER, 1E-4);
        try {
            gradient.getDiagonalHessianLogDensity();
            fail("the Hessian wrt a demographic parameter is not implemented");
        } catch (UnsupportedOperationException expected) {
            // expected
        }
    }

    public void testReportPasses() throws Exception {
        DefaultTreeModel tree = tree(NEWICK);
        Parameter growthRate = new Parameter.Default(-0.3);
        DemographicModel model = new ExponentialGrowthModel(new Parameter.Default(0.4), growthRate,
                Units.Type.YEARS, true);
        CoalescentLikelihood likelihood = new CoalescentLikelihood(new BigFastTreeIntervals(tree), model);
        // throws if the numeric check fails; the negative growth rate must not be clamped at zero
        new CoalescentGradient(likelihood, tree, growthRate, CoalescentGradient.Wrt.PARAMETER, 1E-4).getReport();
        new CoalescentGradient(likelihood, tree, null, CoalescentGradient.Wrt.NODE_HEIGHTS, 1E-4).getReport();
    }

    private static void checkAll(DemographicModel model, String label, Parameter... parameters) throws Exception {
        for (Intervals type : Intervals.values()) {
            DefaultTreeModel tree = tree(NEWICK);
            CoalescentLikelihood likelihood = new CoalescentLikelihood(intervals(tree, type), model);
            checkNodeHeights(tree, likelihood, label + " " + type);
            checkNodeHeightDiagonalHessian(tree, likelihood, label + " " + type);
            for (Parameter parameter : parameters) {
                checkParameter(tree, likelihood, parameter, label + " " + type);
            }
        }
    }

    private static void checkNodeHeights(DefaultTreeModel tree, CoalescentLikelihood likelihood, String label) {
        double[] analytic = new CoalescentGradient(likelihood, tree, null,
                CoalescentGradient.Wrt.NODE_HEIGHTS, 1E-4).getGradientLogDensity();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            final double height = tree.getNodeHeight(node);
            tree.setNodeHeight(node, height + H);
            final double plus = likelihood.getLogLikelihood();
            tree.setNodeHeight(node, height - H);
            final double minus = likelihood.getLogLikelihood();
            tree.setNodeHeight(node, height);
            final double numeric = (plus - minus) / (2 * H);
            assertEquals(label + " node " + node.getNumber(), numeric, analytic[i],
                    TOLERANCE * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static void checkNodeHeightDiagonalHessian(DefaultTreeModel tree, CoalescentLikelihood likelihood,
                                                       String label) {
        CoalescentGradient gradient = new CoalescentGradient(likelihood, tree, null,
                CoalescentGradient.Wrt.NODE_HEIGHTS, 1E-4);
        double[] analytic = gradient.getDiagonalHessianLogDensity();
        final double h = 1E-5;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            final double height = tree.getNodeHeight(node);
            tree.setNodeHeight(node, height + h);
            final double plus = gradient.getGradientLogDensity()[i];
            tree.setNodeHeight(node, height - h);
            final double minus = gradient.getGradientLogDensity()[i];
            tree.setNodeHeight(node, height);
            final double numeric = (plus - minus) / (2 * h);
            assertEquals(label + " Hessian node " + node.getNumber(), numeric, analytic[i],
                    1E-5 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static void checkParameter(DefaultTreeModel tree, CoalescentLikelihood likelihood, Parameter parameter,
                                       String label) {
        final double analytic = new CoalescentGradient(likelihood, tree, parameter,
                CoalescentGradient.Wrt.PARAMETER, 1E-4).getGradientLogDensity()[0];
        final double value = parameter.getParameterValue(0);
        parameter.setParameterValue(0, value + H);
        final double plus = likelihood.getLogLikelihood();
        parameter.setParameterValue(0, value - H);
        final double minus = likelihood.getLogLikelihood();
        parameter.setParameterValue(0, value);
        final double numeric = (plus - minus) / (2 * H);
        assertEquals(label + " parameter value " + value, numeric, analytic,
                TOLERANCE * Math.max(1.0, Math.abs(numeric)));
    }

    private static IntervalList intervals(Tree tree, Intervals type) {
        switch (type) {
            case BIG_FAST:
                return new BigFastTreeIntervals((DefaultTreeModel) tree);
            case MAPPED:
                return new TreeIntervals(tree, true);
            default:
                return new TreeIntervals(tree, false); // BEAUti's default <treeIntervals>
        }
    }

    private static DefaultTreeModel tree(String newick) throws Exception {
        Tree tree = new NewickImporter(newick).importTree(null);
        return new DefaultTreeModel("treeModel", tree);
    }
}
