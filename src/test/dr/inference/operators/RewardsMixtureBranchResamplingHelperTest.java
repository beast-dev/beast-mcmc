/*
 * RewardsMixtureBranchResamplingHelperTest.java
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

import dr.inference.operators.RewardsMixtureBranchResamplingHelper;
import test.dr.math.MathTestCase;

/**
 * Regression tests for branch-local reward-mixture weights.
 *
 * @author Filippo Monti
 */
public class RewardsMixtureBranchResamplingHelperTest extends MathTestCase {

    public void testImpossiblePreorderScaleAnnihilatesAtomicWeight() {
        final double logWeight = RewardsMixtureBranchResamplingHelper.logAtomicWeight(
                0.25,
                1.0,
                -0.5,
                Double.NEGATIVE_INFINITY,
                0.0
        );

        assertEquals(Double.NEGATIVE_INFINITY, logWeight);
    }

    public void testImpossiblePostorderScaleAnnihilatesContinuousWeight() {
        final double[] pre = new double[]{0.25, 0.25, 0.25, 0.25};
        final double[] post = new double[]{1.0, 1.0, 1.0, 1.0};
        final double[] identity = new double[]{
                1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 1.0
        };

        final double logWeight = RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                pre,
                identity,
                post,
                4,
                0.0,
                Double.NEGATIVE_INFINITY
        );

        assertEquals(Double.NEGATIVE_INFINITY, logWeight);
    }

    public void testAllZeroContinuousMatrixAnnihilatesContinuousWeight() {
        final double[] pre = new double[]{0.4, 0.6};
        final double[] post = new double[]{0.7, 0.3};
        final double[] matrix = new double[]{
                0.0, 0.0,
                0.0, 0.0
        };

        final double logWeight = RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                pre,
                matrix,
                post,
                2,
                -1.5,
                0.25
        );

        assertEquals(Double.NEGATIVE_INFINITY, logWeight);
    }

    public void testFiniteScalesAreAddedToAtomicWeight() {
        final double logWeight = RewardsMixtureBranchResamplingHelper.logAtomicWeight(
                0.25,
                0.50,
                -0.5,
                -2.0,
                3.0
        );

        assertEquals(Math.log(0.25) + Math.log(0.50) - 0.5 - 2.0 + 3.0, logWeight, 1.0e-12);
    }

    public void testHugeContinuousMatrixIsScaledBeforeBilinearWeight() {
        final double[] pre = new double[]{0.4, 0.6};
        final double[] post = new double[]{0.7, 0.3};
        final double scale = 2.0e150;
        final double[] matrix = new double[]{
                scale, 0.5 * scale,
                0.25 * scale, 0.75 * scale
        };
        final double preScale = -12.25;
        final double postScale = 3.5;

        final double scaledInner =
                pre[0] * (1.0 * post[0] + 0.5 * post[1]) +
                        pre[1] * (0.25 * post[0] + 0.75 * post[1]);
        final double expected = Math.log(scaledInner) + Math.log(scale) + preScale + postScale;

        final double logWeight = RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                pre,
                matrix,
                post,
                2,
                preScale,
                postScale
        );

        assertEquals(expected, logWeight, 1.0e-12);
    }

    public void testTinyContinuousMatrixIsScaledBeforeBilinearWeight() {
        final double[] pre = new double[]{0.25, 0.75};
        final double[] post = new double[]{0.8, 0.2};
        final double scale = 4.0e-150;
        final double[] matrix = new double[]{
                scale, 0.5 * scale,
                0.25 * scale, 0.125 * scale
        };

        final double scaledInner =
                pre[0] * (1.0 * post[0] + 0.5 * post[1]) +
                        pre[1] * (0.25 * post[0] + 0.125 * post[1]);
        final double expected = Math.log(scaledInner) + Math.log(scale);

        final double logWeight = RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                pre,
                matrix,
                post,
                2
        );

        assertEquals(expected, logWeight, 1.0e-12);
    }

    public void testNonFiniteContinuousMatrixEntryIsRejected() {
        final double[] pre = new double[]{0.4, 0.6};
        final double[] post = new double[]{0.7, 0.3};
        final double[] matrix = new double[]{
                1.0, Double.NaN,
                0.5, 0.25
        };

        try {
            RewardsMixtureBranchResamplingHelper.logContinuousWeight(
                    pre,
                    matrix,
                    post,
                    2
            );
            fail("Expected non-finite matrix entry to be rejected");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("Non-finite continuous transition matrix entry"));
        }
    }
}
