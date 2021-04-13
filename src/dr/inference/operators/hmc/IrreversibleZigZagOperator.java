/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class IrreversibleZigZagOperator extends AbstractZigZagOperator implements Reportable {

    public IrreversibleZigZagOperator(GradientWrtParameterProvider gradientProvider,
                                      PrecisionMatrixVectorProductProvider multiplicationProvider,
                                      PrecisionColumnProvider columnProvider,
                                      double weight, Options runtimeOptions, NativeCodeOptions nativeOptions,
                                      boolean refreshVelocity, Parameter mask,
                                      int threadCount, MassPreconditioner massPreconditioner,
                                      MassPreconditionScheduler.Type preconditionSchedulerType) {

        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, nativeOptions,
                refreshVelocity, mask, threadCount, massPreconditioner, preconditionSchedulerType);
    }


    @Override
    WrappedVector drawInitialVelocity(WrappedVector momentum) {
        if (!refreshVelocity && storedVelocity != null) {
            return storedVelocity;
        } else {
            ReadableVector mass = preconditioning.mass;

            double[] velocity = new double[mass.getDim()];

            for (int i = 0, len = mass.getDim(); i < len; ++i) {
                velocity[i] = (MathUtils.nextDouble() > 0.5) ? 1.0 : -1.0;
            }

            if (mask != null) {
                applyMask(velocity);
            }
            return new WrappedVector.Raw(velocity);
        }
    }

    @Override
    void updatePositionAndMomentum(WrappedVector position,
                                   WrappedVector velocity,
                                   WrappedVector action,
                                   WrappedVector gradient,
                                   WrappedVector momentum,
                                   double time) {
        updatePosition(position.getBuffer(), velocity.getBuffer(), time);
    }

    private MinimumTravelInformation getNextBounceNative(
            WrappedVector position,
            WrappedVector velocity,
            WrappedVector action,
            WrappedVector gradient) {

        if (TIMING) {
            timer.startTimer("getNextC++");
        }

        final MinimumTravelInformation mti = nativeZigZag.getNextIrreversibleEvent(position.getBuffer(), velocity.getBuffer(),
                action.getBuffer(), gradient.getBuffer());

        if (TIMING) {
            timer.stopTimer("getNextC++");
        }

        return mti;
    }

    private void testNative(MinimumTravelInformation firstBounce,
                            WrappedVector position,
                            WrappedVector velocity,
                            WrappedVector action,
                            WrappedVector gradient) {
        if (TIMING) {
            timer.startTimer("getNextC++");
        }

        final MinimumTravelInformation mti;

        mti = nativeZigZag.getNextIrreversibleEvent(position.getBuffer(), velocity.getBuffer(),
                action.getBuffer(), gradient.getBuffer());

        if (TIMING) {
            timer.stopTimer("getNextC++");
        }

        if (!firstBounce.equals(mti)) {
            System.err.println(mti + " ?= " + firstBounce + "\n");
            System.exit(-1);
        }
    }

    @SuppressWarnings("Duplicates")
    @Override
    MinimumTravelInformation getNextBounce(WrappedVector position,
                                           WrappedVector velocity,
                                           WrappedVector action,
                                           WrappedVector gradient,
                                           WrappedVector momentum) {

        final MinimumTravelInformation firstBounce;

        if (TIMING) {
            timer.startTimer("getNext");
        }

        if (nativeCodeOptions.useNativeFindNextBounce) {
            firstBounce = getNextBounceNative(position, velocity, action, gradient);
        } else {
            firstBounce = getNextBounceImpl(position.getBuffer(), velocity.getBuffer(), action.getBuffer(),
                    gradient.getBuffer());
        }

        if (TIMING) {
            timer.stopTimer("getNext");
        }

        if (nativeCodeOptions.testNativeFindNextBounce) {
            testNative(firstBounce, position, velocity, action, gradient);
        }

        return firstBounce;
    }

    @SuppressWarnings("Duplicates")
    private MinimumTravelInformation getNextBounceImpl(double[] p,
                                                       double[] v,
                                                       double[] a,
                                                       double[] g) {

        double minimumTime = Double.POSITIVE_INFINITY;
        int index = -1;
        Type type = Type.NONE;
        double gradientTime;

        if (NEW_WAY) {

            double[] roots = getRoots(a, g);
            double[] rootsSorted = roots.clone();
            Arrays.sort(rootsSorted);

            double minBoundaryTime = Double.POSITIVE_INFINITY;
            gradientTime = getSwitchTimeByMergedProcesses(a, g, v, roots, rootsSorted);
            int gradientIndex = getEventDimension(v, g, a, gradientTime);

            for (int i = 0, end = p.length; i < end; i++) {

                double boundaryTime = findBinaryBoundaryTime(i, p[i], v[i]);

                if (boundaryTime < minBoundaryTime) {
                    minBoundaryTime = boundaryTime;
                    index = i;

                }
            }

            if (gradientTime < minBoundaryTime) {
                index = gradientIndex;
                minimumTime = gradientTime;
                type = Type.GRADIENT;
            } else {
                minimumTime = minBoundaryTime;
                type = Type.BINARY_BOUNDARY;
            }

        } else {

            for (int i = 0, end = p.length; i < end; ++i) {

                double boundaryTime = findBinaryBoundaryTime(i, p[i], v[i]);

                if (boundaryTime < minimumTime) {
                    minimumTime = boundaryTime;
                    index = i;
                    type = Type.BINARY_BOUNDARY;
                }

                double T = MathUtils.nextExponential(1);
                gradientTime = getSwitchTime(-v[i] * g[i], v[i] * a[i], T);

                if (gradientTime < minimumTime) {
                    minimumTime = gradientTime;
                    index = i;
                    type = Type.GRADIENT;
                }
            }
        }
        return new MinimumTravelInformation(minimumTime, index, type);
    }


    private double[] getRoots(double[] action, double[] gradient) {
        // for the linear piece wise line: y = bx + a, where a = -velocity * gradient, b = velocity * action
        // root = - a / b = gradient / action
        double[] roots = new double[action.length];
        double root;
        for (int i = 0; i < action.length; i++) {
            root = gradient[i] / action[i];
            roots[i] = root >= 0 ? root : 0;
        }
        return roots;
    }

    private double getSwitchTimeByMergedProcesses(double[] action, double[] gradient, double[] velocity,
                                                  double[] roots, double[] rootsSorted) {

        // for the linear piece wise line: y = bx + a, where a = -velocity * gradient (as gradient is negated), b =
        // velocity * action
        // root = - a / b = gradient / action
        double T = MathUtils.nextExponential(1);
        double minFirstEvent = -1;
        double cumulativeS = 0;
        PiecewiseLinearEndpoints endPoints;

        if (rootsSorted[rootsSorted.length - 1] == 0) {

            endPoints = getEndpointInfo(0, 0, velocity, gradient, action, roots);
            minFirstEvent = integrateLinearFunctionToArea(endPoints, T, false);
        } else {

            int j = 1;
            double c0;
            double c1;

            double trapezoidArea;
            double residualArea;

            c0 = rootsSorted[0];

            while (j < rootsSorted.length) {
                if (rootsSorted[j] > 0) {
                    c1 = rootsSorted[j];
                    endPoints = getEndpointInfo(c0, c1, velocity, gradient, action, roots);
                    trapezoidArea = getTrapezoidArea(endPoints);
                    cumulativeS += trapezoidArea;
                    if (cumulativeS > T) {
                        residualArea = T - (cumulativeS - trapezoidArea);
                        minFirstEvent = integrateLinearFunctionToArea(endPoints, residualArea, true);
                        break;
                    } else if (j == rootsSorted.length - 1) {
                        residualArea = T - cumulativeS;
                        minFirstEvent = integrateLinearFunctionToArea(endPoints, residualArea, false);
                        break;
                    } else {
                        c0 = c1;
                        j++;
                    }
                } else {
                    j++;
                }
            }

        }

        return minFirstEvent;
    }

    private double integrateLinearFunctionToArea(PiecewiseLinearEndpoints f, double residual, boolean fromC0) {
        double slope;
        double intercept;

        if (fromC0) {
            slope = f.slope0;
            intercept = f.f0 - slope * f.c0;
            return onlyPositiveRoot(slope * 0.5, intercept, -(slope * 0.5 * f.c0 * f.c0 + intercept * f.c0 + residual));
        } else {
            slope = f.slope1;
            intercept = f.f1 - slope * f.c1;
            return onlyPositiveRoot(slope * 0.5, intercept, -(slope * 0.5 * f.c1 * f.c1 + intercept * f.c1 + residual));
        }
    }

    private double onlyPositiveRoot(double a, double b, double c) {
        return (-b + Math.sqrt(b * b - 4 * a * c)) / (2 * a);
    }

    private PiecewiseLinearEndpoints getEndpointInfo(double c0, double c1, double[] velocity, double[] gradient,
                                                     double[] action,
                                                     double[] roots) {

        double[] c0Coefficients = new double[2];
        double[] c1Coefficients = new double[2];

        for (int i = 0; i < roots.length; i++) {

            accumulateCoefficients(c0, velocity, gradient, action, roots, c0Coefficients, i);
            accumulateCoefficients(c1, velocity, gradient, action, roots, c1Coefficients, i);
        }

        return new PiecewiseLinearEndpoints(c0, c1,
                c0Coefficients[0] * c0 + c0Coefficients[1],
                c1Coefficients[0] * c1 + c1Coefficients[1],
                c0Coefficients[0], c1Coefficients[0]);
    }

    private void accumulateCoefficients(double c, double[] velocity, double[] gradient, double[] action, double[] roots,
                                        double[] coefficients, int i) {
        if ((roots[i] >= c && velocity[i] * action[i] <= 0) || (roots[i] <= c && velocity[i] * action[i] >= 0)) {
            coefficients[0] += velocity[i] * action[i];
            coefficients[1] += -velocity[i] * gradient[i];
        }
    }

    private double getTrapezoidArea(PiecewiseLinearEndpoints f) {
        return (f.f0 + f.f1) * (f.c1 - f.c0) / 2.0;
    }

    private int getEventDimension(double[] velocity, double[] gradient, double[] action, double eventTime) {

        double[] proportions = new double[velocity.length];
        double rateAtEvent;
        double rateSum = 0.0;

        for (int i = 0; i < velocity.length; i++) {

            rateAtEvent = eventTime * velocity[i] * action[i] - velocity[i] * gradient[i];
            proportions[i] = rateAtEvent > 0.0 ? rateAtEvent : 0.0;
            rateSum += proportions[i];
        }

        double u = MathUtils.nextDouble();
        double t = 0.0;
        int index = -1;

        for (int i = 0; i < proportions.length; i++) {
            t += proportions[i] / rateSum;
            if (u <= t) {
                index = i;
                break;
            }
        }

        return index;
    }

    private class PiecewiseLinearEndpoints {

        final double c0;
        final double c1;
        final double f0;
        final double f1;
        final double slope0;
        final double slope1;

        private PiecewiseLinearEndpoints(double c0, double c1, double f0, double f1, double slope0, double slope1) {
            this.c0 = c0;
            this.c1 = c1;
            this.f0 = f0;
            this.f1 = f1;
            this.slope0 = slope0;
            this.slope1 = slope1;
        }
    }

    private double getSwitchTime(double a, double b, double u) {
        // simulate T such that P(T>= t) = exp(-at-bt^2/2), using uniform random input u
        if (b > 0) {
            if (a < 0)
                return -a / b + Math.sqrt(2 * u / b);
            else       // a >= 0
                return -a / b + Math.sqrt(a * a / (b * b) + 2 * u / b);
        } else if (b == 0) {
            if (a > 0)
                return u / a;
            else
                return Double.POSITIVE_INFINITY;
        } else { // b  < 0
            if (a <= 0)
                return Double.POSITIVE_INFINITY;
            else {
                // a > 0
                double t1 = -a / b;
                if (u <= a * t1 + b * t1 * t1 / 2)
                    return t1 - Math.sqrt(t1 * t1 + 2 * u / b);
                else
                    return Double.POSITIVE_INFINITY;
            }
        }
    }

    @SuppressWarnings("Duplicates")
    @Override
    void updateDynamics(WrappedVector position,
                        WrappedVector velocity,
                        WrappedVector action,
                        WrappedVector gradient,
                        WrappedVector momentum,
                        WrappedVector column,
                        double time,
                        int[] index,
                        Type eventType) {

        final double[] p = position.getBuffer();
        final double[] v = velocity.getBuffer();
        final double[] a = action.getBuffer();
        final double[] g = gradient.getBuffer();
        final double[] c = column.getBuffer();

        final double twoV = 2 * v[index[0]];

        for (int i = 0, len = p.length; i < len; ++i) {
            final double ai = a[i];

            p[i] = p[i] + time * v[i];
            g[i] = g[i] - time * ai;
            a[i] = ai - twoV * c[i];
        }

        if (NOT_YET_IMPLEMENTED) {
            nativeZigZag.updateIrreversibleDynamics(p, v, a, g, c, time, index[0], eventType.ordinal());
        }
    }

    @Override
    public String getOperatorName() {
        return "Irreversible zig-zag operator";
    }

    static final boolean CPP_NEXT_BOUNCE = false;
    private static final boolean NEW_WAY = false;
    private static final boolean NOT_YET_IMPLEMENTED = false;
}
