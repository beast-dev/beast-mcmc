/*
 * MixedDiscontinuousHamiltonianMonteCarloOperatorTest.java
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

import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.hmc.EmbeddedOrdinalParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.MixedDiscontinuousHamiltonianMonteCarloOperator;
import dr.math.MathUtils;
import dr.util.Transform;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.util.Arrays;

public class MixedDiscontinuousHamiltonianMonteCarloOperatorTest extends TestCase {

    public void testIndependentMixedTargetSampling() {
        MathUtils.setSeed(23L);

        Parameter parameter = new Parameter.Default(new double[]{0.0, 0.25});
        parameter.setId("mixedParameter");

        MixedToyTarget target = new MixedToyTarget(parameter, Math.log(0.25), Math.log(0.75));

        MixedDiscontinuousHamiltonianMonteCarloOperator operator =
                new MixedDiscontinuousHamiltonianMonteCarloOperator(
                        target,
                        target,
                        new boolean[]{false, true},
                        new double[]{1.0, 0.0},
                        new double[]{0.0, 1.0},
                        0.35,
                        8,
                        1.0
                );

        final int burnIn = 4000;
        final int sampleCount = 16000;
        int accepted = 0;
        double sumX = 0.0;
        double sumX2 = 0.0;
        int stateOneCount = 0;

        for (int i = 0; i < burnIn + sampleCount; i++) {
            parameter.storeParameterValues();
            final double oldLogDensity = target.getLogDensity();
            double logHastings = operator.operate(target.getLikelihood());
            final double newLogDensity = target.getLogDensity();
            double acceptProb = Math.min(1.0, Math.exp(newLogDensity - oldLogDensity + logHastings));
            boolean accept = MathUtils.nextDouble() < acceptProb;

            if (accept) {
                parameter.acceptParameterValues();
                operator.accept(logHastings);
            } else {
                parameter.restoreParameterValues();
                operator.reject();
            }

            if (i >= burnIn) {
                if (accept) {
                    accepted++;
                }
                double x = parameter.getParameterValue(0);
                sumX += x;
                sumX2 += x * x;
                if (target.getDiscreteState() == 1) {
                    stateOneCount++;
                }
            }
        }

        double meanX = sumX / sampleCount;
        double varianceX = sumX2 / sampleCount - meanX * meanX;
        double stateOneFrequency = stateOneCount / (double) sampleCount;
        double acceptanceRate = accepted / (double) sampleCount;

        assertEquals(0.0, meanX, 0.08);
        assertEquals(1.0, varianceX, 0.10);
        assertEquals(0.75, stateOneFrequency, 0.035);
        assertTrue("Expected non-trivial acceptance rate", acceptanceRate > 0.6);
    }

    public void testMixedProposalHasFiniteEnergyError() {
        MathUtils.setSeed(7L);

        Parameter parameter = new Parameter.Default(new double[]{0.2, 0.25});
        MixedToyTarget target = new MixedToyTarget(parameter, Math.log(0.25), Math.log(0.75));

        MixedDiscontinuousHamiltonianMonteCarloOperator operator =
                new MixedDiscontinuousHamiltonianMonteCarloOperator(
                        target,
                        target,
                        new boolean[]{false, true},
                        new double[]{1.0, 0.0},
                        new double[]{0.0, 1.0},
                        0.25,
                        4,
                        1.0
                );

        double logHastings = operator.doOperation(target.getLikelihood());
        assertFalse(Double.isNaN(logHastings));
        assertFalse(Double.isInfinite(logHastings));
    }

    public void testTransformedContinuousCoordinateStaysInsideBounds() {
        MathUtils.setSeed(17L);

        Parameter parameter = new Parameter.Default(new double[]{0.5, 0.25});
        parameter.setId("mixedParameter");
        Transform transform = new Transform.Array(Arrays.asList(Transform.LOGIT, Transform.NONE), parameter);

        MixedToyTarget target = new MixedToyTarget(parameter, Math.log(0.25), Math.log(0.75));

        MixedDiscontinuousHamiltonianMonteCarloOperator operator =
                new MixedDiscontinuousHamiltonianMonteCarloOperator(
                        target,
                        target,
                        new boolean[]{false, true},
                        new double[]{1.0, 0.0},
                        new double[]{0.0, 1.0},
                        0.5,
                        4,
                        1.0,
                        2,
                        1E-3,
                        parameter,
                        transform
                );

        for (int i = 0; i < 50; i++) {
            parameter.storeParameterValues();
            double logHastings = operator.operate(target.getLikelihood());
            assertFalse(Double.isNaN(logHastings));
            assertFalse(Double.isInfinite(logHastings));
            assertTrue(parameter.getParameterValue(0) > 0.0);
            assertTrue(parameter.getParameterValue(0) < 1.0);
            parameter.restoreParameterValues();
            operator.reject();
        }
    }

    public void testGradientCheckPassesForCorrectGradient() {
        MathUtils.setSeed(11L);

        Parameter parameter = new Parameter.Default(new double[]{0.3, 1.25});
        MixedToyTarget target = new MixedToyTarget(parameter, Math.log(0.25), Math.log(0.75));

        MixedDiscontinuousHamiltonianMonteCarloOperator operator =
                new MixedDiscontinuousHamiltonianMonteCarloOperator(
                        target,
                        target,
                        new boolean[]{false, true},
                        new double[]{1.0, 0.0},
                        new double[]{0.0, 1.0},
                        0.2,
                        4,
                        1.0,
                        5,
                        1E-3
                );

        for (int i = 0; i < 5; i++) {
            operator.doOperation(target.getLikelihood());
        }
    }

    public void testGradientCheckDetectsIncorrectGradient() {
        MathUtils.setSeed(13L);

        Parameter parameter = new Parameter.Default(new double[]{0.3, 1.25});
        MixedToyTarget target = new MixedToyTarget(parameter, Math.log(0.25), Math.log(0.75));
        WrongGradientWrapper wrongGradient = new WrongGradientWrapper(target);

        MixedDiscontinuousHamiltonianMonteCarloOperator operator =
                new MixedDiscontinuousHamiltonianMonteCarloOperator(
                        wrongGradient,
                        target,
                        new boolean[]{false, true},
                        new double[]{1.0, 0.0},
                        new double[]{0.0, 1.0},
                        0.2,
                        4,
                        1.0,
                        1,
                        1E-3
                );

        try {
            operator.doOperation(target.getLikelihood());
            fail("Expected the gradient check to detect the deliberately incorrect analytic gradient");
        } catch (RuntimeException e) {
            assertTrue(e.getMessage().contains("Gradients do not match"));
        }
    }

    public void testStepSizeIsAdaptableOnLogScale() {
        Parameter parameter = new Parameter.Default(new double[]{0.2, 0.25});
        MixedToyTarget target = new MixedToyTarget(parameter, Math.log(0.25), Math.log(0.75));

        MixedDiscontinuousHamiltonianMonteCarloOperator operator =
                new MixedDiscontinuousHamiltonianMonteCarloOperator(
                        AdaptationMode.ADAPTATION_ON,
                        target,
                        target,
                        new boolean[]{false, true},
                        new double[]{1.0, 0.0},
                        new double[]{0.0, 1.0},
                        0.25,
                        0.75,
                        4,
                        1.0,
                        0,
                        1E-3,
                        parameter,
                        null,
                        0.8
                );

        assertEquals(AdaptationMode.ADAPTATION_ON, operator.getMode());
        assertEquals(0.8, operator.getTargetAcceptanceProbability(), 1E-12);
        assertEquals(0.25, operator.getRawParameter(), 1E-12);
        assertEquals(0.75, operator.getRandomStepSizeFraction(), 1E-12);
        assertEquals(Math.log(0.25), operator.getAdaptableParameter(), 1E-12);

        operator.setAdaptableParameter(Math.log(0.1));

        assertEquals(0.1, operator.getRawParameter(), 1E-12);
    }

    public static Test suite() {
        return new TestSuite(MixedDiscontinuousHamiltonianMonteCarloOperatorTest.class);
    }

    private static class MixedToyTarget implements GradientWrtParameterProvider, DiscontinuousPotentialProvider {

        private final Parameter parameter;
        private final EmbeddedOrdinalParameter embedding;
        private final double[] logStateProbabilities;
        private final Likelihood likelihood;

        MixedToyTarget(Parameter parameter, double leftLogProbability, double rightLogProbability) {
            this.parameter = parameter;
            this.embedding = new EmbeddedOrdinalParameter(new double[]{0.0, 1.0, 2.0});
            this.logStateProbabilities = new double[]{leftLogProbability, rightLogProbability};
            this.likelihood = new Likelihood.Abstract(null) {
                @Override
                protected double calculateLogLikelihood() {
                    return getLogDensity();
                }

                @Override
                protected boolean getLikelihoodKnown() {
                    return false;
                }
            };
        }

        @Override
        public Likelihood getLikelihood() {
            return likelihood;
        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        @Override
        public int getDimension() {
            return parameter.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {
            return new double[]{-parameter.getParameterValue(0), 0.0};
        }

        @Override
        public double getLogDensity() {
            return getLogDensity(parameter.getParameterValues());
        }

        @Override
        public double getLogDensityAfterSingleCoordinateMove(int index, double proposedValue) {
            double[] values = parameter.getParameterValues();
            values[index] = proposedValue;
            return getLogDensity(values);
        }

        @Override
        public double getNextDiscontinuity(int index, double currentValue, double direction) {
            if (index == 1) {
                return embedding.getNextBoundary(currentValue, direction);
            }
            return Double.NaN;
        }

        int getDiscreteState() {
            return embedding.getStateIndex(parameter.getParameterValue(1));
        }

        private double getLogDensity(double[] values) {
            final double x = values[0];
            final double y = values[1];

            if (y < 0.0 || y > 2.0) {
                return Double.NEGATIVE_INFINITY;
            }

            final int state = embedding.getStateIndex(y);
            return -0.5 * x * x + logStateProbabilities[state];
        }
    }

    private static class WrongGradientWrapper implements GradientWrtParameterProvider {

        private final MixedToyTarget target;

        WrongGradientWrapper(MixedToyTarget target) {
            this.target = target;
        }

        @Override
        public Likelihood getLikelihood() {
            return target.getLikelihood();
        }

        @Override
        public Parameter getParameter() {
            return target.getParameter();
        }

        @Override
        public int getDimension() {
            return target.getDimension();
        }

        @Override
        public double[] getGradientLogDensity() {
            double[] gradient = target.getGradientLogDensity();
            gradient[0] += 5.0;
            return gradient;
        }
    }
}
