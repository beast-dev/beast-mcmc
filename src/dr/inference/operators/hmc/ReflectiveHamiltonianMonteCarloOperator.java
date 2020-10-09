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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.GraphicalParameterBound;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inferencexml.operators.hmc.ReflectiveHamiltonianMonteCarloOperatorParser;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;
import dr.xml.Reportable;


/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class ReflectiveHamiltonianMonteCarloOperator extends HamiltonianMonteCarloOperator implements Reportable {

    private final GraphicalParameterBound treeParameterBound;


    public ReflectiveHamiltonianMonteCarloOperator(AdaptationMode mode,
                                                   double weight,
                                                   GradientWrtParameterProvider gradientProvider,
                                                   Parameter parameter,
                                                   Transform transform,
                                                   Parameter maskParameter,
                                                   Options runtimeOptions,
                                                   MassPreconditioner preconditioner,
                                                   GraphicalParameterBound graphicalParameterBound) {

        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner);

        this.treeParameterBound = graphicalParameterBound;
        this.leapFrogEngine = constructLeapFrogEngine(transform);
    }

    @Override
    protected LeapFrogEngine constructLeapFrogEngine(Transform transform) {
        return new WithGraphBounds(parameter, getDefaultInstabilityHandler(), preconditioning, mask, treeParameterBound);
    }

    @Override
    public String getReport() {
        return "operator type: " + ReflectiveHamiltonianMonteCarloOperatorParser.OPERATOR_NAME + "\n\n";
    }

    @Override
    public String getOperatorName() {
        return "ReflectiveHMC(" + parameter.getParameterName() + ")";
    }


    class WithGraphBounds extends HamiltonianMonteCarloOperator.LeapFrogEngine.Default {

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
        public void updatePosition(double[] position, WrappedVector momentum,
                                   double functionalStepSize) {

            double collapsedTime = 0.0;
            while (collapsedTime < functionalStepSize) {
                ReflectionEvent event = nextEvent(position, momentum, functionalStepSize - collapsedTime);
                event.doReflection(position, momentum);
                collapsedTime += event.getEventTime();
            }
            setParameter(position);
        }

        private ReflectionEvent nextEvent(double[] position, WrappedVector momentum, double intervalLength) {
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
            return new ReflectionEvent(type, firstCollisionTime, collisionLocation, intervalLength, new int[]{index1, index2});
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
            return new ReflectionEvent(type, firstReflection, reflectionLocation, intervalLength, new int[]{reflectionIndex});
        }


    }

    class ReflectionEvent {
        private final ReflectionType type;
        private final double eventTime;
        private final double eventLocation;
        private final double intervalLength;
        private final int[] indices;

        ReflectionEvent(ReflectionType type,
                        double eventTime,
                        double eventLocation,
                        double intervalLength,
                        int[] indices) {
            this.type = type;
            this.eventTime = eventTime;
            this.intervalLength = intervalLength;
            this.indices = indices;
            this.eventLocation = eventLocation;
        }

        public double getEventTime() {
            return eventTime;
        }

        public ReflectionType getType() {
            return type;
        }

        public void doReflection(double[] position, WrappedVector momentum) {
            type.doReflection(position, preconditioning, momentum, eventLocation, indices, eventTime);
        }

    }

    enum ReflectionType {
        Reflection {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation, int[] indices, double time) {
                updatePosition(position, preconditioning, momentum, time);
                momentum.set(indices[0], -momentum.get(indices[0]));
                position[indices[0]] = eventLocation;
            }
        },
        Collision {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation, int[] indices, double time) {
                updatePosition(position, preconditioning, momentum, time);
                ReadableVector updatedMomentum = preconditioning.doCollision(indices, momentum);

                for (int index : indices) {
                    momentum.set(index, updatedMomentum.get(index));
                    position[index] = eventLocation;
                }

            }
        },
        None {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation, int[] indices, double time) {
                updatePosition(position, preconditioning, momentum, time);
            }
        };

        void updatePosition(double[] position, MassPreconditioner preconditioning, WrappedVector momentum, double time) {
            final int dim = position.length;
            for (int i = 0; i < dim; i++) {
                position[i] += preconditioning.getVelocity(i, momentum) * time;
            }
        }

        abstract void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                                   double eventLocation, int[] indices, double time);
    }

}
