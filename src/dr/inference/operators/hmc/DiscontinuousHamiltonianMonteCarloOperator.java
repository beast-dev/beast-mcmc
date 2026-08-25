/*
 * DiscontinuousHamiltonianMonteCarloOperator.java
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

package dr.inference.operators.hmc;

import dr.inference.hmc.DiscontinuousPotentialProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.RewardMixturePerformanceStats;
import dr.inference.operators.SimpleMCMCOperator;

/**
 * Initial all-discontinuous HMC operator using Laplace momentum and exact
 * coordinatewise updates.
 *
 * This prototype intentionally focuses on the all-discontinuous case described
 * in Section 4.3 of Nishimura et al. and returns a zero Hastings ratio because
 * the coordinatewise integrator preserves the Hamiltonian exactly.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class DiscontinuousHamiltonianMonteCarloOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final DiscontinuousPotentialProvider provider;
    private final Parameter parameter;
    private final DiscontinuousMomentum momentumHelper;
    private final DiscontinuousCoordinateIntegrator integrator;
    private final double stepSize;
    private final double randomStepSizeFraction;
    private final int nSteps;
    private final int[] updateOrder;

    private long crossingCount = 0;
    private long reflectionCount = 0;

    public DiscontinuousHamiltonianMonteCarloOperator(DiscontinuousPotentialProvider provider,
                                                      double[] momentumScales,
                                                      double stepSize,
                                                      int nSteps,
                                                      double weight) {
        this(provider, momentumScales, stepSize, 0.0, nSteps, weight);
    }

    public DiscontinuousHamiltonianMonteCarloOperator(DiscontinuousPotentialProvider provider,
                                                      double[] momentumScales,
                                                      double stepSize,
                                                      double randomStepSizeFraction,
                                                      int nSteps,
                                                      double weight) {
        this.provider = provider;
        this.parameter = provider.getParameter();
        this.momentumHelper = new DiscontinuousMomentum(momentumScales);
        this.integrator = new DiscontinuousCoordinateIntegrator(provider, momentumHelper);
        this.stepSize = stepSize;
        this.randomStepSizeFraction = randomStepSizeFraction;
        this.nSteps = nSteps;
        this.updateOrder = new int[provider.getDimension()];

        DiscontinuousHmcUtils.validatePositiveStepSize(stepSize);
        DiscontinuousHmcUtils.validateRandomStepSizeFraction(randomStepSizeFraction);
        DiscontinuousHmcUtils.validatePositiveStepCount(nSteps);
        if (provider.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Provider and parameter dimensions must match");
        }

        for (int i = 0; i < updateOrder.length; i++) {
            updateOrder[i] = i;
        }

        setWeight(weight);
    }

    @Override
    public double doOperation() {
        provider.refresh();

        final double[] momentum = momentumHelper.drawMomentum();
        final double stepSizeThisOperation =
                DiscontinuousHmcUtils.drawStepSize(stepSize, randomStepSizeFraction);

        for (int step = 0; step < nSteps; step++) {
            DiscontinuousHmcUtils.shuffleInPlace(updateOrder);

            for (int i = 0; i < updateOrder.length; i++) {
                final int index = updateOrder[i];
                if (!provider.isDiscontinuous(index)) {
                    continue;
                }
                final DiscontinuousCoordinateIntegrator.StepResult result =
                        integrator.step(momentum, index, stepSizeThisOperation);
                RewardMixturePerformanceStats.recordDiscontinuousStep(result.isCrossed(), result.isReflected());

                if (result.isCrossed()) {
                    crossingCount++;
                    RewardMixturePerformanceStats.recordCacheClearAfterAcceptedCrossing();
                    provider.clearOperationCache(
                            RewardMixturePerformanceStats.OperationCacheClearReason.ACCEPTED_CATEGORY_CROSSING);
                } else {
                    RewardMixturePerformanceStats.recordSkippedCacheClearAfterNoCrossing();
                }
                if (result.isReflected()) {
                    reflectionCount++;
                }
            }
        }

        return 0.0;
    }

    @Override
    public String getOperatorName() {
        return "DiscontinuousHMC(" + parameter.getParameterName() + ")";
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    public double getStepSize() {
        return stepSize;
    }

    public int getStepCount() {
        return nSteps;
    }

    public double getRandomStepSizeFraction() {
        return randomStepSizeFraction;
    }

    public long getCrossingCount() {
        return crossingCount;
    }

    public long getReflectionCount() {
        return reflectionCount;
    }
}
