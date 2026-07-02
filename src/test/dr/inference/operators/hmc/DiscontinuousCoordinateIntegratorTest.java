/*
 * DiscontinuousCoordinateIntegratorTest.java
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

package test.dr.inference.operators.hmc;

import dr.inference.distribution.EmbeddedOrdinalLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.DiscontinuousCoordinateIntegrator;
import dr.inference.operators.hmc.DiscontinuousHamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.DiscontinuousMomentum;
import dr.math.MathUtils;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class DiscontinuousCoordinateIntegratorTest extends TestCase {

    private static final double TOL = 1.0e-12;

    public void testReflectionContinuesForRemainingTimeAfterBoundary() {
        final Parameter latent = new Parameter.Default(new double[]{0.9});
        final EmbeddedOrdinalLikelihood provider = makeProvider(latent,
                new double[]{0.0, -10.0},
                new double[]{0.0, 1.0, 2.0});
        final DiscontinuousCoordinateIntegrator integrator = makeIntegrator(provider);
        final double[] momentum = new double[]{2.0};

        final DiscontinuousCoordinateIntegrator.DetailedStepResult result =
                integrator.traceStep(momentum, 0, 0.4);

        assertTrue(result.isReflected());
        assertFalse(result.isCrossed());
        assertEquals(0.7, latent.getParameterValue(0), TOL);
        assertEquals(-2.0, momentum[0], TOL);
        assertEquals(1, result.getEvents().size());
        assertEquals("reflection", result.getEvents().get(0).getType());
        assertEquals(1.0, result.getEvents().get(0).getBoundary(), TOL);
        assertEquals(0.1, result.getEvents().get(0).getElapsedTime(), TOL);
    }

    public void testAcceptedCrossingContinuesForRemainingTimeWithReducedMomentum() {
        final Parameter latent = new Parameter.Default(new double[]{0.9});
        final EmbeddedOrdinalLikelihood provider = makeProvider(latent,
                new double[]{0.0, -1.0},
                new double[]{0.0, 1.0, 2.0});
        final DiscontinuousCoordinateIntegrator integrator = makeIntegrator(provider);
        final double[] momentum = new double[]{2.0};

        final DiscontinuousCoordinateIntegrator.DetailedStepResult result =
                integrator.traceStep(momentum, 0, 0.4);

        assertTrue(result.isCrossed());
        assertFalse(result.isReflected());
        assertEquals(1.3, latent.getParameterValue(0), TOL);
        assertEquals(1.0, momentum[0], TOL);
        assertEquals(1.0, result.getDeltaU(), TOL);
        assertEquals(1, result.getEvents().size());
        assertEquals("crossing", result.getEvents().get(0).getType());
    }

    public void testMultipleBoundariesAreHandledWithinOneStep() {
        final Parameter latent = new Parameter.Default(new double[]{0.9});
        final EmbeddedOrdinalLikelihood provider = makeProvider(latent,
                new double[]{0.0, -0.5, -0.25},
                new double[]{0.0, 1.0, 2.0, 3.0});
        final DiscontinuousCoordinateIntegrator integrator = makeIntegrator(provider);
        final double[] momentum = new double[]{2.0};

        final DiscontinuousCoordinateIntegrator.DetailedStepResult result =
                integrator.traceStep(momentum, 0, 1.4);

        assertTrue(result.isCrossed());
        assertFalse(result.isReflected());
        assertEquals(2.3, latent.getParameterValue(0), TOL);
        assertEquals(1.75, momentum[0], TOL);
        assertEquals(0.25, result.getDeltaU(), TOL);
        assertEquals(2, result.getEvents().size());
        assertEquals(1.0, result.getEvents().get(0).getBoundary(), TOL);
        assertEquals(2.0, result.getEvents().get(1).getBoundary(), TOL);
    }

    public void testAllDiscontinuousOperatorRandomizesStepSizePerProposal() {
        MathUtils.setSeed(20260701L);

        final Parameter latent = new Parameter.Default(new double[]{0.0});
        final EmbeddedOrdinalLikelihood provider = makeProvider(latent,
                new double[]{0.0},
                new double[]{-100.0, 100.0});
        final DiscontinuousHamiltonianMonteCarloOperator operator =
                new DiscontinuousHamiltonianMonteCarloOperator(
                        provider,
                        new double[]{1.0},
                        1.0,
                        1.0,
                        1,
                        1.0);

        operator.doOperation();

        final double displacement = Math.abs(latent.getParameterValue(0));
        assertTrue(displacement >= 0.5 - TOL);
        assertTrue(displacement <= 1.5 + TOL);
        assertTrue(Math.abs(displacement - 1.0) > 1E-10);
    }

    private static EmbeddedOrdinalLikelihood makeProvider(final Parameter latent,
                                                          final double[] logWeights,
                                                          final double[] cuts) {
        return new EmbeddedOrdinalLikelihood(
                latent,
                new Parameter.Default(logWeights),
                new Parameter.Default(cuts));
    }

    private static DiscontinuousCoordinateIntegrator makeIntegrator(final EmbeddedOrdinalLikelihood provider) {
        return new DiscontinuousCoordinateIntegrator(
                provider,
                new DiscontinuousMomentum(new double[]{1.0}));
    }

    public static Test suite() {
        return new TestSuite(DiscontinuousCoordinateIntegratorTest.class);
    }
}
