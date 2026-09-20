/*
 * IGCoalescentLikelihoodTest.java
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

import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.io.NewickImporter;
import dr.evomodel.coalescent.IGCoalescentLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.tree.DefaultTreeModel;
import junit.framework.TestCase;

/**
 * Checks that IGCoalescentLikelihood integrates the population size out of the *joint*
 * density of the genealogy, rather than multiplying together per-interval marginals (which
 * would ignore the fact that every interval shares the same, now-marginalised, population size).
 */
public class IGCoalescentLikelihoodTest extends TestCase {

    private static final double TOLERANCE = 1E-6;

    /**
     * A 3-taxon tree gives two coalescent intervals, with lineage counts 3 and 2 and
     * durations 0.5 and 0.5. The joint marginal density can be reduced by hand to a
     * single integral over an inverse power of theta, which is verified independently
     * below by numerical quadrature rather than by re-deriving the closed form under test.
     */
    public void testJointMarginalMatchesNumericalIntegration() throws Exception {

        NewickImporter importer = new NewickImporter("(A:1.0,(B:0.5,C:0.5):0.5);");
        DefaultTreeModel tree = new DefaultTreeModel(importer.importTree(null));
        IntervalList intervals = new TreeIntervals(tree);

        final double alpha = 2.0;
        final double beta = 3.0;

        final double logLikelihood = IGCoalescentLikelihood.calculateLogLikelihood(intervals, alpha, beta);

        final double expected = Math.log(numericallyIntegratedJointDensity(intervals, alpha, beta));

        assertEquals(expected, logLikelihood, TOLERANCE);
    }

    /**
     * Confirms that multiplying the single-interval marginal densities together (the bug
     * this class replaces) gives a different answer from the correct joint calculation,
     * on the same tree used above.
     */
    public void testJointDiffersFromProductOfPerIntervalMarginals() throws Exception {

        NewickImporter importer = new NewickImporter("(A:1.0,(B:0.5,C:0.5):0.5);");
        DefaultTreeModel tree = new DefaultTreeModel(importer.importTree(null));
        IntervalList intervals = new TreeIntervals(tree);

        final double alpha = 2.0;
        final double beta = 3.0;

        final double joint = IGCoalescentLikelihood.calculateLogLikelihood(intervals, alpha, beta);

        double productOfMarginals = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            final double duration = intervals.getInterval(i);
            final int lineageCount = intervals.getLineageCount(i);
            final double kChoose2 = lineageCount * (lineageCount - 1.0) / 2.0;
            productOfMarginals += IGCoalescentLikelihood.marginalLogDensity(duration, kChoose2, alpha, beta);
        }

        assertTrue(Math.abs(joint - productOfMarginals) > 1E-3);
    }

    /**
     * Brute-force reference: integrates the exact joint density of the intervals given theta,
     * against the inverse-Gamma(alpha, beta) prior on theta, by substituting u = 1/theta so the
     * domain becomes finite-ish and Simpson's rule converges quickly.
     */
    private static double numericallyIntegratedJointDensity(IntervalList intervals, double alpha, double beta) {

        final int steps = 200000;
        final double uMax = 200.0; // theta = 1/u ranges down to 0.005, far into the tail
        final double h = uMax / steps;

        double sum = 0.0;
        for (int i = 0; i <= steps; i++) {
            final double u = i * h;
            final double weight = (i == 0 || i == steps) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
            sum += weight * integrandInU(u, intervals, alpha, beta);
        }
        return sum * h / 3.0;
    }

    private static double integrandInU(double u, IntervalList intervals, double alpha, double beta) {
        if (u == 0.0) {
            return 0.0;
        }
        final double theta = 1.0 / u;

        double jointDensityGivenTheta = 1.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            final double duration = intervals.getInterval(i);
            final int lineageCount = intervals.getLineageCount(i);
            final double kChoose2 = lineageCount * (lineageCount - 1.0) / 2.0;

            // every interval contributes the "no coalescence yet" survival factor;
            // only an interval that actually ends in a coalescent event also contributes
            // the rate at which that event occurs.
            jointDensityGivenTheta *= Math.exp(-kChoose2 * duration / theta);
            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                jointDensityGivenTheta *= kChoose2 / theta;
            }
        }

        final double igPdf = Math.pow(beta, alpha) / gamma(alpha)
                * Math.pow(theta, -alpha - 1.0) * Math.exp(-beta / theta);

        // dtheta = -du/u^2
        return jointDensityGivenTheta * igPdf / (u * u);
    }

    private static double gamma(double x) {
        return Math.exp(dr.math.GammaFunction.lnGamma(x));
    }
}
