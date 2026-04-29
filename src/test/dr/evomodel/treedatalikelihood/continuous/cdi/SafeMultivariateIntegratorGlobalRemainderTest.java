/*
 * SafeMultivariateIntegratorGlobalRemainderTest.java
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

package dr.evomodel.treedatalikelihood.continuous.cdi;

import junit.framework.TestCase;

public class SafeMultivariateIntegratorGlobalRemainderTest extends TestCase {

    public void testGlobalRemainderAdjointAllowsZeroLengthObservedBranch() {
        final int dim = 1;
        final SafeMultivariateIntegrator integrator = new SafeMultivariateIntegrator(
                PrecisionType.FULL,
                1,
                dim,
                dim,
                2,
                1);

        integrator.setDiffusionPrecision(0, new double[]{1.0}, 0.0);
        integrator.updateBrownianDiffusionMatrices(
                0,
                new int[]{0},
                new double[]{0.0},
                null,
                1);

        integrator.setPostOrderPartial(0, observedPartial(1.0));
        integrator.setPostOrderPartial(1, observedPartial(1.0));

        final double[] adjoint = new double[dim * dim];
        integrator.computeGlobalRemainderAdjointWrtBranchVariance(0, 0, 1, adjoint);

        assertEquals(0.0, adjoint[0], 0.0);
    }

    private static double[] observedPartial(final double mean) {
        final int dim = 1;
        final double[] partial = new double[PrecisionType.FULL.getPartialsDimension(dim)];
        partial[PrecisionType.FULL.getMeanOffset(dim)] = mean;
        partial[PrecisionType.FULL.getPrecisionOffset(dim)] = Double.POSITIVE_INFINITY;
        partial[PrecisionType.FULL.getVarianceOffset(dim)] = 0.0;
        PrecisionType.FULL.fillEffDimInPartials(partial, 0, dim, dim);
        PrecisionType.FULL.fillNoDeterminantInPartials(partial, 0, dim);
        PrecisionType.FULL.fillRemainderInPartials(partial, 0, 0.0, dim);
        return partial;
    }
}
