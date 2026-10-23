/*
 * CoalescentGrowthGradientTest.java
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
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.demographicmodel.ExponentialGrowthModel;
import dr.evomodel.coalescent.demographicmodel.LogisticGrowthModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Node-height gradients of the coalescent under exponential and logistic growth with a growth rate large enough
 * for the -log N(t) term to matter (plan 3, B3). The likelihood subtracts log N at each coalescent time, so the
 * gradient must subtract d log N / dt.
 */
public class CoalescentGrowthGradientTest extends TestCase {

    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";

    public CoalescentGrowthGradientTest(String name) {
        super(name);
    }

    public void testExponentialGrowth() throws Exception {
        for (double r : new double[]{0.5, -0.8, 3.0}) {
            DefaultTreeModel tree = tree();
            DemographicModel model = new ExponentialGrowthModel(new Parameter.Default(0.3),
                    new Parameter.Default(r), Units.Type.YEARS, true);
            checkNodeHeights(tree, model, "exponential r=" + r);
            checkLogDemographicGradient(model.getDemographicFunction(), "exponential r=" + r);
        }
    }

    public void testLogisticGrowth() throws Exception {
        DefaultTreeModel tree = tree();
        DemographicModel model = new LogisticGrowthModel(new Parameter.Default(0.3), new Parameter.Default(2.0),
                new Parameter.Default(0.2), 0.0, Units.Type.YEARS, true);
        checkNodeHeights(tree, model, "logistic");
        checkLogDemographicGradient(model.getDemographicFunction(), "logistic");
    }

    private static void checkNodeHeights(DefaultTreeModel tree, DemographicModel model, String label) {
        CoalescentLikelihood likelihood = new CoalescentLikelihood(new BigFastTreeIntervals(tree), model);
        CoalescentGradient gradient = new CoalescentGradient(likelihood, tree, null,
                CoalescentGradient.Wrt.NODE_HEIGHTS, 1E-4);
        double[] analytic = gradient.getGradientLogDensity();

        // analytic[i] is internal node externalCount + i
        final double h = 1E-6;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            final double height = tree.getNodeHeight(node);
            tree.setNodeHeight(node, height + h);
            final double plus = likelihood.getLogLikelihood();
            tree.setNodeHeight(node, height - h);
            final double minus = likelihood.getLogLikelihood();
            tree.setNodeHeight(node, height);
            final double numeric = (plus - minus) / (2 * h);
            assertEquals(label + " node " + node.getNumber(), numeric, analytic[i],
                    1E-6 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    // getLogDemographicGradient is d log N / dt, the convention the gradient now relies on
    private static void checkLogDemographicGradient(DemographicFunction demographic, String label) {
        final double h = 1E-6;
        for (double t : new double[]{0.05, 0.2, 0.4}) {
            final double numeric = (Math.log(demographic.getDemographic(t + h))
                    - Math.log(demographic.getDemographic(t - h))) / (2 * h);
            assertEquals(label + " t=" + t, numeric, demographic.getLogDemographicGradient(t), 1E-7);
        }
    }

    private static DefaultTreeModel tree() throws Exception {
        Tree newick = new NewickImporter(NEWICK).importTree(null);
        return new DefaultTreeModel("treeModel", newick);
    }
}
