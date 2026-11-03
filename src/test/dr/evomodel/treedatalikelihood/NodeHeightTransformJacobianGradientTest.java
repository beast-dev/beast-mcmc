/*
 * NodeHeightTransformJacobianGradientTest.java
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

import dr.evolution.io.NewickImporter;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightTransform;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * The full ratio transform's (withRoot) gradient with the Jacobian: updateGradientLogDensity against the finite
 * difference of l(h(y)) - logJacobian(h(y)), the root coordinate included (its Jacobian term had the wrong sign).
 *
 * @author Xiang Ji
 */
public class NodeHeightTransformJacobianGradientTest extends TestCase {

    private static NodeHeightTransform transform(String newick) throws Exception {
        TreeModel tree = new DefaultTreeModel("treeModel", new NewickImporter(newick).importTree(null));
        Parameter ratios = new Parameter.Default(tree.getInternalNodeCount() - 1);
        ratios.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ratios.getDimension()));
        return new NodeHeightTransform(ratios, tree, new StrictClockBranchRates(new Parameter.Default(1.0)), true);
    }

    // l(h) = sum c_k h_k, or 0
    private static double density(double[] h, boolean zero) {
        double sum = 0.0;
        for (int k = 0; !zero && k < h.length; ++k) {
            sum += (k % 2 == 0 ? 1.0 : -1.0) * (0.35 + 0.17 * k) * h[k];
        }
        return sum;
    }

    private static double[] densityGradient(int n, boolean zero) {
        double[] gradient = new double[n];
        for (int k = 0; !zero && k < n; ++k) {
            gradient[k] = (k % 2 == 0 ? 1.0 : -1.0) * (0.35 + 0.17 * k);
        }
        return gradient;
    }

    private static double[] analytic(NodeHeightTransform transform, boolean zero) {
        double[] h = transform.getNodeHeights().getParameterValues();
        return transform.updateGradientLogDensity(densityGradient(h.length, zero), h, 0, h.length);
    }

    private static void checkAgainstFiniteDifferences(String newick, boolean zero) throws Exception {
        NodeHeightTransform transform = transform(newick);
        double[] h0 = transform.getNodeHeights().getParameterValues();
        final int n = h0.length;
        double[] y0 = transform.transform(h0, 0, n);
        transform.inverse(y0, 0, n);
        double[] analytic = analytic(transform, zero);
        for (int j = 0; j < n; ++j) {
            final double step = 1e-6 * Math.max(1.0, Math.abs(y0[j]));
            double[] y = y0.clone();
            y[j] = y0[j] + step;
            double[] h = transform.inverse(y, 0, n);
            final double plus = density(h, zero) - transform.logJacobian(h, 0, n);
            y[j] = y0[j] - step;
            h = transform.inverse(y, 0, n);
            final double minus = density(h, zero) - transform.logJacobian(h, 0, n);
            transform.inverse(y0, 0, n);
            assertEquals(newick + " coordinate " + j, (plus - minus) / (2 * step), analytic[j], 1e-6);
        }
    }

    // contemporaneous tips, l = 0: -logJacobian = sum log(h_i) over the non-root internal nodes, h_i = H * prod(ratios),
    // so its derivative in the root height is (internal nodes - 1) / H
    public void testRootCoordinateClosedForm() throws Exception {
        assertEquals(1.0 / 2.0, analytic(transform("((A:1,B:1):1,C:2);"), true)[0], 1e-12);
        assertEquals(5.0 / 3.5, analytic(transform(
                "((((A:1,B:1):0.5,C:1.5):1.2,(D:0.7,E:0.7):2.0):0.8,(F:2.1,G:2.1):1.4);"), true)[0], 1e-12);
    }

    public void testDatedTipsAgainstFiniteDifferences() throws Exception {
        final String dated = "((((A:1.2,B:0.8):0.8,C:0.9):1.0,(D:1.25,E:0.6):1.5):1.0,((F:0.8,G:2.4):0.9,H:1.0):0.7);";
        checkAgainstFiniteDifferences(dated, true);
        checkAgainstFiniteDifferences(dated, false);
        checkAgainstFiniteDifferences("(((((A:1.0,B:0.4):0.8,C:0.5):0.6,(D:0.7,F:0.4):1.5):0.7,E:1.1):0.5,G:0.7);", false);
    }
}
