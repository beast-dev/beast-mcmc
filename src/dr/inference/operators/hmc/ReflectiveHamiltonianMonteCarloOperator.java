/*
 * ReflectiveHamiltonianMonteCarloOperator.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inference.operators.AdaptationMode;
import dr.inferencexml.operators.hmc.ReflectiveHamiltonianMonteCarloOperatorParser;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.logging.Logger;


/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class ReflectiveHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator implements Reportable {

    private final GeneralBoundsProvider parameterBound;
    private boolean isAtBoundary = false;
    private long boundsRejectionCount = 0;
    private long zeroTimeCapCount = 0;
    // true only inside leapFrogGivenMomentum, so split HMC (reversiblePositionMomentumUpdate) is never checked
    private boolean checkBoundsInTrajectory = false;
    private boolean tiedConstraintsChecked = false;
    private static final boolean DEBUG = false;

    // The cap on one run of consecutive events that do not advance the drift, which stops an endless sequence of
    // them. Away from exact ties a run is at most a few events; from a tied start it is as long as the start
    // requires (about k^2 / 2 for a fully tied chain of k nodes resolved one wall per event, about pi * sqrt(r) for
    // a two-wall wedge with inverse-mass ratio r), so a start whose run exceeds the cap cannot be left.
    private static final long ZERO_TIME_EVENT_CAP = 1000000L;


    public ReflectiveHamiltonianMonteCarloOperator(AdaptationMode mode,
                                                   double weight,
                                                   GradientWrtParameterProvider gradientProvider,
                                                   Parameter parameter,
                                                   Transform transform,
                                                   Parameter maskParameter,
                                                   Options runtimeOptions,
                                                   MassPreconditioner preconditioner,
                                                   GeneralBoundsProvider bounds) {

        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner);

        if (maskParameter != null) {
            throw new IllegalArgumentException("Reflective HMC does not support masks");
        }

        if (bounds instanceof GraphicalParameterBound) {
            final GraphicalParameterBound graphicalBounds = (GraphicalParameterBound) bounds;
            graphicalBounds.checkParameterMapping();
            graphicalBounds.checkParameterMapping(parameter);
        }

        this.parameterBound = bounds;
        this.leapFrogEngine = constructLeapFrogEngine(transform);
    }

    @Override
    protected LeapFrogEngine constructLeapFrogEngine(Transform transform) {
        if (parameterBound == null) {
            return null; //will get called again
        }


        if (transform != null) {
            throw new RuntimeException("not yet implemented");
        }

        if (parameterBound instanceof GraphicalParameterBound) { //TODO: don't use 'instanceof' to deal with this.
            return new WithGraphBounds(parameter, getDefaultInstabilityHandler(), preconditioning, mask,
                    (GraphicalParameterBound) parameterBound);
        }
        return new WithMultivariateCurvedBounds(parameter, getDefaultInstabilityHandler(), preconditioning, mask,
                (BoundedSpace) parameterBound);
    }

    @Override
    public String getReport() {
        return "operator type: " + ReflectiveHamiltonianMonteCarloOperatorParser.OPERATOR_NAME + "\n" +
                "bounds-check rejections: " + boundsRejectionCount + "\n" +
                "zero-time event cap hits: " + zeroTimeCapCount + "\n\n";
    }

    public long getBoundsRejectionCount() {
        return boundsRejectionCount;
    }

    public long getZeroTimeCapCount() {
        return zeroTimeCapCount;
    }

    @Override
    public String getOperatorName() {
        return "ReflectiveHMC(" + parameter.getParameterName() + ")";
    }

    @Override
    public void reversiblePositionMomentumUpdate(WrappedVector position, WrappedVector momentum,
                                                 WrappedVector gradient, int direction, double time) {
        if (direction != 1) {
            throw new UnsupportedOperationException(
                    "Reflective HMC cannot integrate backwards in time, so it does not support NUTS");
        }
        super.reversiblePositionMomentumUpdate(position, momentum, gradient, direction, time);
    }

    @Override
    protected double leapFrogGivenMomentum(WrappedVector momentum) throws NumericInstabilityException {

        final double[] currentValues = parameter.getParameterValues();
        if (!parameterBound.isWithinBounds(currentValues)) {
            throw new IllegalStateException("The current value of " + parameter.getParameterName() +
                    " is outside the bounds, so reflective HMC cannot start from it");
        }

        if (!tiedConstraintsChecked) {
            tiedConstraintsChecked = true;
            final int tied = parameterBound.countTiedConstraints(currentValues);
            if (tied > 0) {
                Logger.getLogger("dr.inference").info("WARNING: The value of " +
                        parameter.getParameterName() + " at reflective HMC's first proposal has " + tied +
                        (tied == 1 ? " tied constraint" : " tied constraints") +
                        " (equal parent and child heights, or values exactly on a bound). Reflective HMC resolves" +
                        " such ties by reflecting at time zero; long chains of tied nodes or very unequal masses can" +
                        " take many such events, and an unusually long run of them rejects the proposal.");
            }
        }

        checkBoundsInTrajectory = true;
        try {
            return super.leapFrogGivenMomentum(momentum);
        } finally {
            checkBoundsInTrajectory = false;
        }
    }

    private static boolean isFinite(double[] values) {
        for (double value : values) {
            if (!Double.isFinite(value)) {
                return false;
            }
        }
        return true;
    }


    abstract class WithBounds extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

        WithBounds(Parameter parameter,
                   InstabilityHandler instabilityHandler,
                   MassPreconditioner preconditioning,
                   double[] mask) {
            super(parameter, instabilityHandler, preconditioning, mask);
        }

        protected abstract ReflectionEvent nextEvent(double[] position, WrappedVector momentum, double intervalLength,
                                                     boolean isAtBoundary) throws NumericInstabilityException;

        @Override
        public void updatePosition(double[] position, WrappedVector momentum,
                                   double functionalStepSize) throws NumericInstabilityException {

            long zeroTimeEvents = 0;
            double collapsedTime = 0.0;
            while (collapsedTime < functionalStepSize) {
                ReflectionEvent event = nextEvent(position, momentum, functionalStepSize - collapsedTime, isAtBoundary);

                isAtBoundary = event.doReflection(position, momentum);
                final double advancedTime = collapsedTime + event.getEventTime();

                // An event that does not advance the drift (time zero, or too short to change collapsedTime) counts
                // toward the cap; an event that advances it resets the count
                if (advancedTime > collapsedTime) {
                    zeroTimeEvents = 0;
                } else if (++zeroTimeEvents > ZERO_TIME_EVENT_CAP) {
                    ++zeroTimeCapCount;
                    throw new NumericInstabilityException();
                }
                collapsedTime = advancedTime;
            }

            if (checkBoundsInTrajectory && (!isFinite(position) || !parameterBound.isWithinBounds(position))) {
                ++boundsRejectionCount;
                throw new NumericInstabilityException();
            }
            setParameter(position);
        }


    }

    class WithMultivariateCurvedBounds extends WithBounds {

        private final BoundedSpace space;
        public final int[] defaultIndices;

        WithMultivariateCurvedBounds(Parameter parameter,
                                     InstabilityHandler instabilityHandler,
                                     MassPreconditioner preconditioning,
                                     double[] mask,
                                     BoundedSpace space) {
            super(parameter, instabilityHandler, preconditioning, mask);
            this.space = space;
            ArrayList<Integer> inds = new ArrayList<>();
            for (int i = 0; i < parameter.getDimension(); i++) {
                if (mask == null || mask[i] == 1) {
                    inds.add(i);
                }
            }
            this.defaultIndices = new int[inds.size()];
            for (int i = 0; i < inds.size(); i++) {
                defaultIndices[i] = inds.get(i);
            }
        }


        @Override
        protected ReflectionEvent nextEvent(double[] position, WrappedVector momentum, double intervalLength,
                                            boolean isAtBoundary) throws NumericInstabilityException {
            double[] velocity = preconditioning.getVelocity(momentum);
            double timeToReflection = space.forwardDistanceToBoundary(position, velocity, isAtBoundary);

            if (DEBUG) {
                System.out.println("Time to reflection: " + timeToReflection);
                System.out.println("Interval length: " + intervalLength);
            }


            if (timeToReflection > intervalLength) {
                return new ReflectionEvent(ReflectionType.None, intervalLength, Double.NaN, new int[0]);
            } else {
                double[] boundaryPosition = new double[position.length];
                for (int i = 0; i < position.length; i++) {
                    boundaryPosition[i] = position[i] + timeToReflection * velocity[i];
                }
                double[] normalVector = space.getNormalVectorAtBoundary(boundaryPosition);
                return new ReflectionEvent(ReflectionType.MultivariateReflection, timeToReflection,
                        intervalLength - timeToReflection,
                        boundaryPosition, normalVector, defaultIndices);
            }
        }
    }


    class WithGraphBounds extends WithBounds {

        final private GraphicalParameterBound graphicalParameterBound;

        protected WithGraphBounds(Parameter parameter,
                                  InstabilityHandler instabilityHandler,
                                  MassPreconditioner preconditioning,
                                  double[] mask,
                                  GraphicalParameterBound graphicalParameterBound) {

            super(parameter, instabilityHandler, preconditioning, mask);

            this.graphicalParameterBound = graphicalParameterBound;

        }


        @Override
        protected ReflectionEvent nextEvent(double[] position, WrappedVector momentum, double intervalLength,
                                            boolean isAtBoundary) {
            ReflectionEvent reflectionEventAtFixedBound = firstReflectionAtFixedBounds(position, momentum, intervalLength);
            ReflectionEvent collisionEvent = firstCollision(position, momentum, intervalLength);
            return (reflectionEventAtFixedBound.getEventTime() < collisionEvent.getEventTime()) ? reflectionEventAtFixedBound : collisionEvent;
        }

        private ReflectionEvent firstCollision(double[] position, ReadableVector momentum, double intervalLength) {

            final int dim = position.length;
            double firstCollisionTime = intervalLength;
            double collisionLocation = -1.0;
            ReflectionType type = ReflectionType.None;
            int index1 = -1;
            int index2 = -1;
            for (int i = 0; i < dim; i++) {
                final double velocity1 = preconditioning.getVelocity(i, momentum);
                int[] connectedParameterIndices = graphicalParameterBound.getConnectedParameterIndices(i);
                if (connectedParameterIndices != null) {
                    for (int j : graphicalParameterBound.getConnectedParameterIndices(i)) {
                        if (j > i) {
                            final double velocity2 = preconditioning.getVelocity(j, momentum);
                            // The wall x_i == x_j is hit only when the motion points out of the constrained space
                            // (closing > 0); on the wall (gap == 0), or past it by rounding, that is an event at time 0.
                            // Motion inward or along the wall is no event, and so is a NaN gap (a NaN coordinate, or
                            // two infinite ones of the same sign).
                            final boolean iBelow = graphicalParameterBound.isBelow(i, j);
                            final double gap = iBelow ? position[j] - position[i] : position[i] - position[j];
                            final double closing = iBelow ? velocity1 - velocity2 : velocity2 - velocity1;
                            if (closing > 0) {
                                final double collisionTime;
                                if (gap > 0) {
                                    collisionTime = (position[j] - position[i]) / (velocity1 - velocity2);
                                } else if (gap <= 0) {
                                    collisionTime = 0.0;
                                } else {
                                    continue;
                                }
                                if (collisionTime < firstCollisionTime) {
                                    firstCollisionTime = collisionTime;
                                    collisionLocation = collisionTime * velocity1 + position[i];
                                    index1 = i;
                                    index2 = j;
                                    type = ReflectionType.Collision;
                                }
                            }
                        }
                    }
                }
            }
            return new ReflectionEvent(type, firstCollisionTime, collisionLocation, new int[]{index1, index2});
        }

        private ReflectionEvent firstReflectionAtFixedBounds(double[] position, ReadableVector momentum, double intervalLength) {
            final int dim = position.length;
            double firstReflection = intervalLength;
            double reflectionLocation = -1.0;
            ReflectionType type = ReflectionType.None;
            int reflectionIndex = -1;
            for (int i = 0; i < dim; i++) {
                final double velocity = preconditioning.getVelocity(i, momentum);
                // The velocity picks the bound that can be hit; on that bound, or past it by rounding, the
                // reflection is an event at time 0. An infinite bound gives an infinite time, so no event, and a
                // NaN gap (a NaN coordinate, or an infinite one at an infinite bound) is no event either.
                if (velocity > 0) {
                    final double upperBound = graphicalParameterBound.getFixedUpperBound(i);
                    final double gap = upperBound - position[i];
                    final double reflectionTime;
                    if (gap > 0) {
                        reflectionTime = (upperBound - position[i]) / velocity;
                    } else if (gap <= 0) {
                        reflectionTime = 0.0;
                    } else {
                        continue;
                    }
                    if (reflectionTime < firstReflection) {
                        firstReflection = reflectionTime;
                        type = ReflectionType.Reflection;
                        reflectionIndex = i;
                        reflectionLocation = upperBound;
                    }
                } else if (velocity < 0) {
                    final double lowerBound = graphicalParameterBound.getFixedLowerBound(i);
                    final double gap = position[i] - lowerBound;
                    final double reflectionTime;
                    if (gap > 0) {
                        reflectionTime = (lowerBound - position[i]) / velocity;
                    } else if (gap <= 0) {
                        reflectionTime = 0.0;
                    } else {
                        continue;
                    }
                    if (reflectionTime < firstReflection) {
                        firstReflection = reflectionTime;
                        type = ReflectionType.Reflection;
                        reflectionIndex = i;
                        reflectionLocation = lowerBound;
                    }
                }
            }
            return new ReflectionEvent(type, firstReflection, reflectionLocation, new int[]{reflectionIndex});
        }


    }

    class ReflectionEvent {
        private final ReflectionType type;
        private final double eventTime;
        private final double[] eventLocation;
        private final int[] indices;
        private final double[] normalVector;
        private final double remainingTime;

        ReflectionEvent(ReflectionType type,
                        double eventTime,
                        double remainingTime,
                        double[] eventLocation,
                        double[] normalVector,
                        int[] indices) {
            this.type = type;
            this.eventTime = eventTime;
            this.indices = indices;
            this.eventLocation = eventLocation;
            this.normalVector = normalVector;
            this.remainingTime = remainingTime;
        }

        ReflectionEvent(ReflectionType type,
                        double eventTime,
                        double eventLocation,
                        int[] indices) {
            this(type, eventTime, Double.NaN, new double[]{eventLocation}, null, indices);
        }


        public double getEventTime() {
            return eventTime;
        }

        public ReflectionType getType() {
            return type;
        }

        public boolean doReflection(double[] position, WrappedVector momentum) {

            if (DEBUG) {
                System.out.println("time: " + eventTime);
                System.out.print("start: ");
                System.out.println(new WrappedVector.Raw(position));
                System.out.println(momentum);
            }

            type.doReflection(position, preconditioning, momentum, eventLocation, indices, normalVector, eventTime, remainingTime);

            if (DEBUG) {
                System.out.print("end: ");
                System.out.println(new WrappedVector.Raw(position));
                System.out.println(momentum);
            }

            return type.isAtBoundary();
        }

    }

    enum ReflectionType {

        MultivariateReflection {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation[], int[] indices, double[] normalVector, double time,
                              double remainingTime) {


                updatePosition(position, preconditioning, momentum, time);
                double vn = 0;
                double nn = 0;

                for (int i : indices) {
                    vn += momentum.get(i) * normalVector[i];
                    nn += normalVector[i] * normalVector[i];
                }

                double c = 2 * vn / nn;

                for (int i : indices) {
                    momentum.set(i, momentum.get(i) - c * normalVector[i]);
                    position[i] = eventLocation[i];
                }


            }
        },
        Reflection {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation[], int[] indices, double[] normalVector, double time,
                              double remainingTime) {
                updatePosition(position, preconditioning, momentum, time);
                momentum.set(indices[0], -momentum.get(indices[0]));
                position[indices[0]] = eventLocation[0];  // TODO: XJ: they should be the same, change it to assert?
            }
        },
        Collision {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation[], int[] indices, double[] normalVector, double time,
                              double remainingTime) {
                updatePosition(position, preconditioning, momentum, time);
                ReadableVector updatedMomentum = preconditioning.doCollision(indices, momentum);

                for (int index : indices) {
                    momentum.set(index, updatedMomentum.get(index));
                    position[index] = eventLocation[0];   // TODO: XJ: they should be the same, change it to assert?
                }

            }
        },
        None {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation[], int[] indices, double[] normalVector, double time,
                              double remainginTime) {

                updatePosition(position, preconditioning, momentum, time);

            }

            @Override
            public boolean isAtBoundary() {
                return false;
            }
        };

        void updatePosition(double[] position, MassPreconditioner preconditioning, WrappedVector momentum, double time) {
            final int dim = position.length;
            for (int i = 0; i < dim; i++) {
                position[i] += preconditioning.getVelocity(i, momentum) * time;
            }
        }

        abstract void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                                   double eventLocation[], int[] indices, double[] normalVector, double time, double remainingTime);


        public boolean isAtBoundary() {
            return true;
        }

    }

}
