/*
 * MixedDiscontinuousHamiltonianMonteCarloOperator.java
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
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.markovjumps.RewardDensityDomainException;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractAdaptableOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.GeneralOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.util.Transform;

/**
 * Mixed smooth/discontinuous HMC prototype implementing the Strang splitting
 * described in Nishimura et al.:
 *
 * 1. half smooth momentum update
 * 2. half smooth position update
 * 3. exact coordinatewise discontinuous updates
 * 4. half smooth position update
 * 5. half smooth momentum update
 *
 * This version uses diagonal unit-style masses specified explicitly by arrays.
 *
 * This is a {@link GeneralOperator}: BEAST's {@code MarkovChain} dispatches it
 * through {@link #doOperation(Likelihood)}, passing the actual joint posterior
 * ({@code <joint>}/{@code <posterior>}). That joint is used directly for the
 * outer Metropolis evaluation. Like {@code HamiltonianMonteCarloOperator}, this
 * operator returns the momentum/Jacobian Hastings contribution; BEAST's acceptor
 * supplies the {@code newScore - oldScore} term from the joint likelihood. The
 * discontinuous provider still needs correct per-coordinate potential
 * differences for the coordinatewise integrator. The no-arg
 * {@link #doOperation()} is never called by {@code MarkovChain} for a
 * {@link GeneralOperator} and is intentionally unsupported.
 *
 * If {@code gradientCheckCount > 0}, the first {@code gradientCheckCount}
 * operations numerically verify {@code smoothGradientProvider}'s analytic
 * gradient against {@code joint}, restricted to the non-discontinuous
 * coordinates (see {@link GradientCheckUtils}).
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class MixedDiscontinuousHamiltonianMonteCarloOperator extends AbstractAdaptableOperator implements GeneralOperator {

    public static final double DEFAULT_TARGET_ACCEPTANCE_PROBABILITY = 0.8;

    private final GradientWrtParameterProvider smoothGradientProvider;
    private final DiscontinuousPotentialProvider discontinuousProvider;
    private final Parameter parameter;
    private final Transform transform;
    private final boolean[] discontinuousMask;
    private final double[] continuousMasses;
    private final double[] discontinuousScales;
    private double stepSize;
    private final double randomStepSizeFraction;
    private final int nSteps;
    private final int[] discontinuousOrder;
    private final DiscontinuousCoordinateIntegrator discontinuousIntegrator;
    private final boolean[] continuousMask;
    private final int gradientCheckCount;
    private final double gradientCheckTolerance;

    public MixedDiscontinuousHamiltonianMonteCarloOperator(GradientWrtParameterProvider smoothGradientProvider,
                                                           DiscontinuousPotentialProvider discontinuousProvider,
                                                           boolean[] discontinuousMask,
                                                           double[] continuousMasses,
                                                           double[] discontinuousScales,
                                                           double stepSize,
                                                           int nSteps,
                                                           double weight) {
        this(smoothGradientProvider, discontinuousProvider, discontinuousMask, continuousMasses,
                discontinuousScales, stepSize, 0.0, nSteps, weight, 0, 1E-3, null, null);
    }

    public MixedDiscontinuousHamiltonianMonteCarloOperator(GradientWrtParameterProvider smoothGradientProvider,
                                                           DiscontinuousPotentialProvider discontinuousProvider,
                                                           boolean[] discontinuousMask,
                                                           double[] continuousMasses,
                                                           double[] discontinuousScales,
                                                           double stepSize,
                                                           int nSteps,
                                                           double weight,
                                                           int gradientCheckCount,
                                                           double gradientCheckTolerance) {
        this(smoothGradientProvider, discontinuousProvider, discontinuousMask, continuousMasses,
                discontinuousScales, stepSize, 0.0, nSteps, weight, gradientCheckCount,
                gradientCheckTolerance, null, null);
    }

    public MixedDiscontinuousHamiltonianMonteCarloOperator(GradientWrtParameterProvider smoothGradientProvider,
                                                           DiscontinuousPotentialProvider discontinuousProvider,
                                                           boolean[] discontinuousMask,
                                                           double[] continuousMasses,
                                                           double[] discontinuousScales,
                                                           double stepSize,
                                                           int nSteps,
                                                           double weight,
                                                           int gradientCheckCount,
                                                           double gradientCheckTolerance,
                                                           Parameter parameter,
                                                           Transform transform) {
        this(smoothGradientProvider, discontinuousProvider, discontinuousMask, continuousMasses,
                discontinuousScales, stepSize, 0.0, nSteps, weight, gradientCheckCount,
                gradientCheckTolerance, parameter, transform);
    }

    public MixedDiscontinuousHamiltonianMonteCarloOperator(GradientWrtParameterProvider smoothGradientProvider,
                                                           DiscontinuousPotentialProvider discontinuousProvider,
                                                           boolean[] discontinuousMask,
                                                           double[] continuousMasses,
                                                           double[] discontinuousScales,
                                                           double stepSize,
                                                           double randomStepSizeFraction,
                                                           int nSteps,
                                                           double weight,
                                                           int gradientCheckCount,
                                                           double gradientCheckTolerance,
                                                           Parameter parameter,
                                                           Transform transform) {
        this(AdaptationMode.DEFAULT, smoothGradientProvider, discontinuousProvider, discontinuousMask,
                continuousMasses, discontinuousScales, stepSize, randomStepSizeFraction, nSteps, weight,
                gradientCheckCount, gradientCheckTolerance, parameter, transform,
                DEFAULT_TARGET_ACCEPTANCE_PROBABILITY);
    }

    public MixedDiscontinuousHamiltonianMonteCarloOperator(AdaptationMode adaptationMode,
                                                           GradientWrtParameterProvider smoothGradientProvider,
                                                           DiscontinuousPotentialProvider discontinuousProvider,
                                                           boolean[] discontinuousMask,
                                                           double[] continuousMasses,
                                                           double[] discontinuousScales,
                                                           double stepSize,
                                                           int nSteps,
                                                           double weight,
                                                           int gradientCheckCount,
                                                           double gradientCheckTolerance,
                                                           Parameter parameter,
                                                           Transform transform,
                                                           double targetAcceptanceProbability) {
        this(adaptationMode, smoothGradientProvider, discontinuousProvider, discontinuousMask,
                continuousMasses, discontinuousScales, stepSize, 0.0, nSteps, weight,
                gradientCheckCount, gradientCheckTolerance, parameter, transform,
                targetAcceptanceProbability);
    }

    public MixedDiscontinuousHamiltonianMonteCarloOperator(AdaptationMode adaptationMode,
                                                           GradientWrtParameterProvider smoothGradientProvider,
                                                           DiscontinuousPotentialProvider discontinuousProvider,
                                                           boolean[] discontinuousMask,
                                                           double[] continuousMasses,
                                                           double[] discontinuousScales,
                                                           double stepSize,
                                                           double randomStepSizeFraction,
                                                           int nSteps,
                                                           double weight,
                                                           int gradientCheckCount,
                                                           double gradientCheckTolerance,
                                                           Parameter parameter,
                                                           Transform transform,
                                                           double targetAcceptanceProbability) {
        super(adaptationMode, targetAcceptanceProbability);
        this.smoothGradientProvider = smoothGradientProvider;
        this.discontinuousProvider = discontinuousProvider;
        this.parameter = parameter == null ? discontinuousProvider.getParameter() : parameter;
        this.transform = transform;
        this.discontinuousMask = discontinuousMask.clone();
        this.continuousMasses = continuousMasses.clone();
        this.discontinuousScales = discontinuousScales.clone();
        this.stepSize = stepSize;
        this.randomStepSizeFraction = randomStepSizeFraction;
        this.nSteps = nSteps;
        this.gradientCheckCount = gradientCheckCount;
        this.gradientCheckTolerance = gradientCheckTolerance;

        final int dimension = this.parameter.getDimension();
        if (smoothGradientProvider.getDimension() != dimension ||
                discontinuousProvider.getDimension() != dimension ||
                discontinuousMask.length != dimension ||
                continuousMasses.length != dimension ||
                discontinuousScales.length != dimension) {
            throw new IllegalArgumentException("All mixed-HMC inputs must match the parameter dimension");
        }
        if (transform instanceof Transform.MultivariableTransform &&
                ((Transform.MultivariableTransform) transform).getDimension() != dimension) {
            throw new IllegalArgumentException("Transform dimension must match the mixed-HMC parameter dimension");
        }
        DiscontinuousHmcUtils.validatePositiveStepSize(stepSize);
        DiscontinuousHmcUtils.validateRandomStepSizeFraction(randomStepSizeFraction);
        DiscontinuousHmcUtils.validatePositiveStepCount(nSteps);

        int discontinuousCount = 0;
        for (int i = 0; i < dimension; i++) {
            if (this.discontinuousMask[i]) {
                if (!(this.discontinuousScales[i] > 0.0)) {
                    throw new IllegalArgumentException("Discontinuous scales must be positive on discontinuous coordinates");
                }
                discontinuousCount++;
            } else if (!(this.continuousMasses[i] > 0.0)) {
                throw new IllegalArgumentException("Continuous masses must be positive on continuous coordinates");
            }
        }

        this.discontinuousOrder = new int[discontinuousCount];
        this.continuousMask = new boolean[dimension];
        int k = 0;
        for (int i = 0; i < dimension; i++) {
            if (this.discontinuousMask[i]) {
                discontinuousOrder[k++] = i;
            } else {
                continuousMask[i] = true;
            }
        }

        double[] integratorScales = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            integratorScales[i] = this.discontinuousMask[i] ? this.discontinuousScales[i] : 1.0;
        }

        this.discontinuousIntegrator = new DiscontinuousCoordinateIntegrator(
                discontinuousProvider, new DiscontinuousMomentum(integratorScales));

        setWeight(weight);
    }

    @Override
    public double doOperation() {
        throw new RuntimeException("Should not be executed");
    }

    @Override
    public double doOperation(Likelihood joint) {
        discontinuousProvider.refresh();

        // A reward-mixture branch model can reach a mathematically invalid
        // state partway through a trajectory (e.g. an unconstrained smooth
        // coordinate landing arbitrarily close to an atomic reward value --
        // see REWARD_CATEGORY_DYNAMIC_ORDERING_PLAN.md). Rather than crash
        // the whole chain, treat this exactly like BEAST's other HMC
        // operators treat numerical instability during a trajectory (see
        // HamiltonianMonteCarloOperator.NumericInstabilityException): reject
        // this one proposal via -Infinity and let MarkovChain restore state,
        // leaving the chain to continue from its last valid position.
        try {
            if (getCount() < gradientCheckCount) {
                GradientCheckUtils.checkGradient(joint, parameter, smoothGradientProvider.getGradientLogDensity(),
                        continuousMask, gradientCheckTolerance, transform);
            }

            final double[] continuousMomentum = drawContinuousMomentum();
            final double[] discontinuousMomentum = drawDiscontinuousMomentum();

            final double initialAuxiliaryEnergy =
                    getContinuousKineticEnergy(continuousMomentum) +
                    getDiscontinuousKineticEnergy(discontinuousMomentum) +
                    getParameterLogJacobian();
            final double stepSizeThisOperation =
                    DiscontinuousHmcUtils.drawStepSize(stepSize, randomStepSizeFraction);

            halfStepSmoothMomentum(continuousMomentum, stepSizeThisOperation);

            for (int step = 0; step < nSteps; step++) {
                halfStepSmoothPosition(continuousMomentum, stepSizeThisOperation);
                DiscontinuousHmcUtils.shuffleInPlace(discontinuousOrder);
                for (int i = 0; i < discontinuousOrder.length; i++) {
                    discontinuousIntegrator.step(discontinuousMomentum, discontinuousOrder[i], stepSizeThisOperation);
                }
                halfStepSmoothPosition(continuousMomentum, stepSizeThisOperation);

                if (step < nSteps - 1) {
                    fullStepSmoothMomentum(continuousMomentum, stepSizeThisOperation);
                }
            }

            halfStepSmoothMomentum(continuousMomentum, stepSizeThisOperation);

            // The discontinuous landscape is intentionally frozen during the
            // trajectory, but the proposal must be scored under the embedding
            // implied by the final smooth coordinates.
            discontinuousProvider.refreshAfterPositionUpdate();
            joint.makeDirty();

            final double finalAuxiliaryEnergy =
                    getContinuousKineticEnergy(continuousMomentum) +
                    getDiscontinuousKineticEnergy(discontinuousMomentum) +
                    getParameterLogJacobian();

            return initialAuxiliaryEnergy - finalAuxiliaryEnergy;
        } catch (RewardDensityDomainException e) {
            // The exception can be thrown from deep inside a tree-likelihood
            // traversal (e.g. RewardsAwareBranchModel.computeCtsTransitionMatrices
            // via a boundary-crossing trial evaluation). Force everything dirty
            // so the next evaluation fully recomputes rather than trusting
            // whatever partial state the aborted traversal left behind, then
            // reject this proposal exactly like BEAST's other HMC operators
            // treat numerical instability (see
            // HamiltonianMonteCarloOperator.NumericInstabilityException).
            //
            // Note: MarkovChain's own post-reject sanity check ("state was not
            // correctly restored"), active only during its early full-evaluation
            // warm-up window (<mcmc fullEvaluation="..."/>, default 1000 states),
            // can still flag a discrepancy here and terminate the chain even
            // though this reject itself is handled correctly -- see
            // REWARD_CATEGORY_DYNAMIC_ORDERING_PLAN.md. Until that is root-caused,
            // XML using this operator with a reward-mixture model should set
            // fullEvaluation="0" to skip that warm-up check.
            joint.makeDirty();
            return Double.NEGATIVE_INFINITY;
        }
    }

    @Override
    public String getOperatorName() {
        return "MixedDiscontinuousHMC(" + parameter.getParameterName() + ")";
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(stepSize);
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
        stepSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return stepSize;
    }

    public double getRandomStepSizeFraction() {
        return randomStepSizeFraction;
    }

    @Override
    public String getAdaptableParameterName() {
        return "stepSize";
    }

    private double[] drawContinuousMomentum() {
        final double[] momentum = new double[parameter.getDimension()];
        for (int i = 0; i < momentum.length; i++) {
            if (!discontinuousMask[i]) {
                momentum[i] = Math.sqrt(continuousMasses[i]) * MathUtils.nextGaussian();
            }
        }
        return momentum;
    }

    private double[] drawDiscontinuousMomentum() {
        final double[] momentum = new double[parameter.getDimension()];
        for (int i = 0; i < momentum.length; i++) {
            if (discontinuousMask[i]) {
                final double sign = MathUtils.nextBoolean() ? 1.0 : -1.0;
                momentum[i] = sign * discontinuousScales[i] * MathUtils.nextExponential(1.0);
            }
        }
        return momentum;
    }

    private double getContinuousKineticEnergy(double[] momentum) {
        double total = 0.0;
        for (int i = 0; i < momentum.length; i++) {
            if (!discontinuousMask[i]) {
                total += momentum[i] * momentum[i] / continuousMasses[i];
            }
        }
        return total / 2.0;
    }

    private double getDiscontinuousKineticEnergy(double[] momentum) {
        double total = 0.0;
        for (int i = 0; i < momentum.length; i++) {
            if (discontinuousMask[i]) {
                total += Math.abs(momentum[i]) / discontinuousScales[i];
            }
        }
        return total;
    }

    private double getParameterLogJacobian() {
        return transform == null ? 0.0 :
                transform.logJacobian(parameter.getParameterValues(), 0, parameter.getDimension());
    }

    private void halfStepSmoothMomentum(double[] momentum, double operationStepSize) {
        final double[] gradient = getSmoothGradientInOperatorCoordinates();
        for (int i = 0; i < momentum.length; i++) {
            if (!discontinuousMask[i]) {
                momentum[i] += 0.5 * operationStepSize * gradient[i];
            }
        }
    }

    private void fullStepSmoothMomentum(double[] momentum, double operationStepSize) {
        final double[] gradient = getSmoothGradientInOperatorCoordinates();
        for (int i = 0; i < momentum.length; i++) {
            if (!discontinuousMask[i]) {
                momentum[i] += operationStepSize * gradient[i];
            }
        }
    }

    private void halfStepSmoothPosition(double[] momentum, double operationStepSize) {
        final double[] position = getOperatorPosition();
        for (int i = 0; i < momentum.length; i++) {
            if (!discontinuousMask[i]) {
                final double velocity = momentum[i] / continuousMasses[i];
                position[i] += 0.5 * operationStepSize * velocity;
            }
        }
        // Check before writing to the model: a NaN/Infinite position written
        // into the parameter here would only surface later as an exception
        // thrown from deep inside a tree-likelihood traversal (see
        // REWARD_CATEGORY_DYNAMIC_ORDERING_PLAN.md), by which point BEAST's
        // reject-and-restore cannot cleanly recover (the traversal's own
        // per-node dirty-bit tracking is left mid-update). Catching it here,
        // before any model state is touched, keeps the reject clean.
        checkFinite(position, "smooth position update");
        setOperatorPosition(position);
    }

    private double[] getSmoothGradientInOperatorCoordinates() {
        final double[] gradient = smoothGradientProvider.getGradientLogDensity();
        final double[] result;
        if (transform == null) {
            result = gradient;
        } else {
            result = transform.updateGradientLogDensity(gradient, parameter.getParameterValues(),
                    0, parameter.getDimension());
        }
        checkFinite(result, "smooth gradient");
        return result;
    }

    private void checkFinite(double[] values, String context) {
        for (double value : values) {
            if (Double.isNaN(value) || Double.isInfinite(value)) {
                throw new RewardDensityDomainException(
                        "Non-finite value encountered in " + context + ": " + value);
            }
        }
    }

    private double[] getOperatorPosition() {
        final double[] position = parameter.getParameterValues();
        if (transform == null) {
            return position;
        }
        return transform.transform(position, 0, position.length);
    }

    private void setOperatorPosition(double[] position) {
        final double[] rawPosition = transform == null ? position :
                transform.inverse(position, 0, position.length);
        ReadableVector.Utils.setParameter(rawPosition, parameter);
    }
}
