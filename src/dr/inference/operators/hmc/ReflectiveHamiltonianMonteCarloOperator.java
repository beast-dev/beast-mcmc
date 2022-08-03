/*
 * DiscontinuousHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.bss.Utils;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inference.operators.AdaptationMode;
import dr.inferencexml.operators.hmc.ReflectiveHamiltonianMonteCarloOperatorParser;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;

import java.util.ArrayList;


/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class ReflectiveHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator implements Reportable {

    private final GeneralBoundsProvider parameterBound;
    private boolean isAtBoundary = false;
    private static final boolean DEBUG = false;


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
        return "operator type: " + ReflectiveHamiltonianMonteCarloOperatorParser.OPERATOR_NAME + "\n\n";
    }

    @Override
    public String getOperatorName() {
        return "ReflectiveHMC(" + parameter.getParameterName() + ")";
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

            double collapsedTime = 0.0;
            while (collapsedTime < functionalStepSize) {
                ReflectionEvent event = nextEvent(position, momentum, functionalStepSize - collapsedTime, isAtBoundary);

                isAtBoundary = event.doReflection(position, momentum);
                collapsedTime += event.getEventTime();
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

        private boolean isReflected(double position, double intendedNewPosition, double bound) {
            if (position > bound) {
                return intendedNewPosition <= bound;
            } else if (position < bound) {
                return intendedNewPosition >= bound;
            } else {
                return false;
            }
        }

        private boolean isCollision(double position1, double intendedPosition1,
                                    double position2, double intendedPosition2) {
            if (position1 > position2) {
                return intendedPosition1 <= intendedPosition2;
            } else if (position1 < position2) {
                return intendedPosition1 >= intendedPosition2;
            } else {
                return false;
            }
        }

        private ReflectionEvent firstCollision(double[] position, ReadableVector momentum, double intervalLength) {

            final int dim = position.length;
            double[] intendedNewPosition = getIntendedPosition(position, momentum, intervalLength);
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
                            if (isCollision(position[i], intendedNewPosition[i], position[j], intendedNewPosition[j])) {
                                final double collisionTime = (position[j] - position[i]) / (velocity1 - velocity2);
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

        private double[] getIntendedPosition(double[] position, ReadableVector momentum, double intervalLength) {
            final int dim = position.length;
            double[] intendedNewPosition = new double[dim];
            for (int i = 0; i < dim; i++) {
                final double velocity = preconditioning.getVelocity(i, momentum);
                intendedNewPosition[i] = position[i] + intervalLength * velocity;
            }
            return intendedNewPosition;
        }

        private ReflectionEvent firstReflectionAtFixedBounds(double[] position, ReadableVector momentum, double intervalLength) {
            final int dim = position.length;
            double[] intendedNewPosition = getIntendedPosition(position, momentum, intervalLength);
            double firstReflection = intervalLength;
            double reflectionLocation = -1.0;
            ReflectionType type = ReflectionType.None;
            int reflectionIndex = -1;
            for (int i = 0; i < dim; i++) {
                final double velocity = preconditioning.getVelocity(i, momentum);
                final double upperBound = graphicalParameterBound.getFixedUpperBound(i);
                final double lowerBound = graphicalParameterBound.getFixedLowerBound(i);
                if (isReflected(position[i], intendedNewPosition[i], upperBound)) {
                    final double reflectionTime = (upperBound - position[i]) / velocity;
                    if (reflectionTime < 0.0) {
                        throw new RuntimeException("Check isReflected() function plz.");
                    } else if (reflectionTime < firstReflection) {
                        firstReflection = reflectionTime;
                        type = ReflectionType.Reflection;
                        reflectionIndex = i;
                        reflectionLocation = upperBound;
                    }
                } else if (isReflected(position[i], intendedNewPosition[i], lowerBound)) {
                    final double reflectionTime = (lowerBound - position[i]) / velocity;
                    if (reflectionTime < 0.0) {
                        throw new RuntimeException("Check isReflected() function plz.");
                    } else if (reflectionTime < firstReflection) {
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
                Utils.printArray(position);
                System.out.println(momentum);
            }

            type.doReflection(position, preconditioning, momentum, eventLocation, indices, normalVector, eventTime, remainingTime);

            if (DEBUG) {
                System.out.print("end: ");
                Utils.printArray(position);
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
                position[indices[0]] = eventLocation[0];
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
                    position[index] = eventLocation[0];
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
