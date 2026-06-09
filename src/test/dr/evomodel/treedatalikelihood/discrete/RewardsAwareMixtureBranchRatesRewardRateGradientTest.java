/*
 * RewardsAwareMixtureBranchRatesRewardRateGradientTest.java
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

package test.dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.treedatalikelihood.discrete.RewardsAwareMixtureBranchRatesRewardRateGradient;
import dr.inference.model.Parameter;
import test.dr.math.MathTestCase;

/**
 * Regression tests for the atom-state to reward-rate-value mapping used by the
 * dependent reward-rate gradient.
 *
 * @author Filippo Monti
 */
public class RewardsAwareMixtureBranchRatesRewardRateGradientTest extends MathTestCase {

    public void testAtomStateIsMappedThroughRewardRatesMapping() {
        Parameter atomIndices = new Parameter.Default(new double[]{0.0, 1.0, 2.0, 3.0});
        Parameter rewardRatesMapping = new Parameter.Default(new double[]{0.0, 2.0, 3.0, 1.0});
        Parameter rewardRateValues = new Parameter.Default(new double[]{0.0, 1.0, 0.25, 0.75});

        assertEquals(0, map(atomIndices, rewardRatesMapping, rewardRateValues, 0));
        assertEquals(2, map(atomIndices, rewardRatesMapping, rewardRateValues, 1));
        assertEquals(3, map(atomIndices, rewardRatesMapping, rewardRateValues, 2));
        assertEquals(1, map(atomIndices, rewardRatesMapping, rewardRateValues, 3));
    }

    public void testMappedRewardRateIndexGivesVaryingValueOffset() {
        Parameter atomIndices = new Parameter.Default(new double[]{1.0, 2.0, 3.0});
        Parameter rewardRatesMapping = new Parameter.Default(new double[]{0.0, 2.0, 3.0, 1.0});
        Parameter rewardRateValues = new Parameter.Default(new double[]{0.0, 1.0, 0.25, 0.75});

        assertEquals(0, map(atomIndices, rewardRatesMapping, rewardRateValues, 0) - 2);
        assertEquals(1, map(atomIndices, rewardRatesMapping, rewardRateValues, 1) - 2);
        assertEquals(-1, map(atomIndices, rewardRatesMapping, rewardRateValues, 2) - 2);
    }

    public void testRejectsNonIntegerAtomStateIndex() {
        Parameter atomIndices = new Parameter.Default(new double[]{1.5});
        Parameter rewardRatesMapping = new Parameter.Default(new double[]{0.0, 2.0});
        Parameter rewardRateValues = new Parameter.Default(new double[]{0.0, 1.0, 0.25});

        assertInvalid(atomIndices, rewardRatesMapping, rewardRateValues, 0);
    }

    public void testRejectsMappedIndexOutsideRewardRateValues() {
        Parameter atomIndices = new Parameter.Default(new double[]{1.0});
        Parameter rewardRatesMapping = new Parameter.Default(new double[]{0.0, 3.0});
        Parameter rewardRateValues = new Parameter.Default(new double[]{0.0, 1.0, 0.25});

        assertInvalid(atomIndices, rewardRatesMapping, rewardRateValues, 0);
    }

    private static int map(
            Parameter atomIndices,
            Parameter rewardRatesMapping,
            Parameter rewardRateValues,
            int parameterIndex) {

        return RewardsAwareMixtureBranchRatesRewardRateGradient.mapAtomStateToRewardRateIndex(
                atomIndices,
                rewardRatesMapping,
                rewardRateValues,
                parameterIndex);
    }

    private static void assertInvalid(
            Parameter atomIndices,
            Parameter rewardRatesMapping,
            Parameter rewardRateValues,
            int parameterIndex) {

        try {
            map(atomIndices, rewardRatesMapping, rewardRateValues, parameterIndex);
            fail("Expected invalid reward-rate mapping to be rejected");
        } catch (IllegalStateException expected) {
            // expected
        }
    }
}
