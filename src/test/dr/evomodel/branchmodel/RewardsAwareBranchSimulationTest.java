/*
 * RewardsAwareBranchSimulationTest.java
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

import dr.evomodel.branchmodel.RewardsAwareBranchSimulation;
import test.dr.math.MathTestCase;

import java.util.Random;

/**
 * Smoke tests for reward-conditioned branch path simulation.
 *
 * @author Filippo Monti
 */
public class RewardsAwareBranchSimulationTest extends MathTestCase {

    private static final double TOL = 1.0e-8;

    public void testTwoKnotRewardDensityIsUniform() {
        assertEquals(1.0,
                RewardsAwareBranchSimulation.dividedDifferenceDensity(new double[]{0.0, 1.0}, 0.5),
                1.0e-12);
    }

    public void testTwoStateConditionalPathHonorsEndpointsTimeAndReward() {
        double[][] q = {
                {-2.0, 2.0},
                {3.0, -3.0}
        };
        double[] rewards = {0.0, 1.0};
        double branchLength = 1.0;
        double rho = 0.35;

        RewardsAwareBranchSimulation.Options options = new RewardsAwareBranchSimulation.Options();
        options.nParticles = 500;
        options.maxParticleMultiplier = 2;
        options.essThreshold = 0.01;
        options.component = RewardsAwareBranchSimulation.Component.CONTINUOUS;

        RewardsAwareBranchSimulation sampler =
                new RewardsAwareBranchSimulation(q, rewards, 1, branchLength);
        RewardsAwareBranchSimulation.PathSample sample =
                sampler.sample(0, rho, options, new Random(17));

        assertEquals(0, sample.getStates()[0]);
        assertEquals(1, sample.getStates()[sample.getStates().length - 1]);
        assertEquals(branchLength, sample.getTotalTime(), TOL);
        assertEquals(rho, sample.getRewardProportion(rewards), TOL);
        assertTrue(sample.getEffectiveSampleSizeFraction() > 0.0);
    }

    public void testAutoBoundaryRewardUsesAtomPath() {
        double[][] q = {
                {-1.5, 1.0, 0.5},
                {0.7, -1.0, 0.3},
                {0.2, 0.4, -0.6}
        };
        double[] rewards = {0.0, 0.0, 1.0};

        RewardsAwareBranchSimulation.Options options = new RewardsAwareBranchSimulation.Options();
        options.component = RewardsAwareBranchSimulation.Component.AUTO;

        RewardsAwareBranchSimulation sampler =
                new RewardsAwareBranchSimulation(q, rewards, 1, 0.75);
        RewardsAwareBranchSimulation.PathSample sample =
                sampler.sample(0, 0.0, options, new Random(19));

        int[] states = sample.getStates();
        for (int state : states) {
            assertTrue(state == 0 || state == 1);
        }
        assertEquals(0, states[0]);
        assertEquals(1, states[states.length - 1]);
        assertEquals(0.0, sample.getRewardProportion(rewards), TOL);
    }
}
