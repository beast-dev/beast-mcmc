/*
 * DiscontinuousCoordinateIntegrator.java
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Coordinatewise discontinuous integrator corresponding to Algorithm 1 in the
 * Nishimura et al. paper.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class DiscontinuousCoordinateIntegrator {

    private static final double TIME_TOLERANCE = 1.0e-14;
    private static final int MAX_BOUNDARY_EVENTS = 100000;

    protected final DiscontinuousPotentialProvider provider;
    protected final Parameter parameter;
    protected final DiscontinuousMomentum momentumHelper;

    public DiscontinuousCoordinateIntegrator(DiscontinuousPotentialProvider provider,
                                             DiscontinuousMomentum momentumHelper) {
        this.provider = provider;
        this.parameter = provider.getParameter();
        this.momentumHelper = momentumHelper;

        if (provider.getDimension() != momentumHelper.getDimension()) {
            throw new IllegalArgumentException("Provider and momentum dimensions must match");
        }
    }

    public StepResult step(double[] momentum, int index, double stepSize) {
        return traceStep(momentum, index, stepSize).toStepResult();
    }

    public DetailedStepResult traceStep(double[] momentum, int index, double stepSize) {
        final double initialValue = parameter.getParameterValue(index);
        final double initialMomentum = momentum[index];
        final List<SubStepEvent> events = new ArrayList<>();

        if (DiscontinuousMomentum.sign(initialMomentum) == 0.0 || stepSize == 0.0) {
            return new DetailedStepResult(index, initialValue, initialValue, initialMomentum, initialMomentum,
                    0.0, false, false, events);
        }

        double currentValue = initialValue;
        double currentMomentum = initialMomentum;
        double remainingTime = stepSize;
        double acceptedDeltaU = 0.0;
        boolean crossed = false;
        boolean reflected = false;
        int eventCount = 0;

        while (remainingTime > TIME_TOLERANCE && DiscontinuousMomentum.sign(currentMomentum) != 0.0) {
            if (++eventCount > MAX_BOUNDARY_EVENTS) {
                throw new IllegalStateException("Exceeded maximum discontinuous-HMC boundary events");
            }

            final double velocity = momentumHelper.getVelocity(currentMomentum, index);
            final double direction = DiscontinuousMomentum.sign(velocity);
            final double targetValue = currentValue + remainingTime * velocity;
            final double boundary = provider.getNextDiscontinuity(index, currentValue, direction);

            if (!isReachableBoundary(currentValue, targetValue, boundary, direction)) {
                currentValue = targetValue;
                remainingTime = 0.0;
                break;
            }

            final double timeToBoundary = (boundary - currentValue) / velocity;
            if (!(timeToBoundary >= 0.0) || timeToBoundary > remainingTime + TIME_TOLERANCE) {
                currentValue = targetValue;
                remainingTime = 0.0;
                break;
            }

            remainingTime = Math.max(0.0, remainingTime - timeToBoundary);

            final double deltaU = provider.getPotentialDifferenceAcrossBoundary(index, boundary, direction);
            final double availableKinetic = momentumHelper.getKineticEnergy(currentMomentum, index);

            if (Double.isFinite(deltaU) && availableKinetic >= deltaU) {
                final double newMomentum =
                        momentumHelper.updateMomentumAfterCrossing(currentMomentum, index, deltaU);
                events.add(new SubStepEvent("crossing", currentValue, boundary, boundary,
                        currentMomentum, newMomentum, deltaU, timeToBoundary));
                acceptedDeltaU += deltaU;
                currentMomentum = newMomentum;
                currentValue = representativeValuePastBoundary(boundary, direction);
                crossed = true;
            } else {
                final double newMomentum = momentumHelper.reflectMomentum(currentMomentum);
                events.add(new SubStepEvent("reflection", currentValue, boundary, boundary,
                        currentMomentum, newMomentum, deltaU, timeToBoundary));
                currentMomentum = newMomentum;
                currentValue = representativeValuePastBoundary(boundary, -direction);
                reflected = true;
            }
        }

        parameter.setParameterValue(index, currentValue);
        momentum[index] = currentMomentum;

        return new DetailedStepResult(index, initialValue, currentValue, initialMomentum, currentMomentum,
                acceptedDeltaU, crossed, reflected, events);
    }

    private static boolean isReachableBoundary(double currentValue,
                                               double targetValue,
                                               double boundary,
                                               double direction) {
        if (Double.isNaN(boundary)) {
            return false;
        }
        if (direction > 0.0) {
            return boundary > currentValue && boundary <= targetValue;
        } else if (direction < 0.0) {
            return boundary < currentValue && boundary >= targetValue;
        }
        return false;
    }

    private static double representativeValuePastBoundary(double boundary, double direction) {
        return Math.nextAfter(boundary, boundary + direction);
    }

    public static class StepResult {
        private final int index;
        private final double oldPosition;
        private final double newPosition;
        private final double oldMomentum;
        private final double newMomentum;
        private final double deltaU;
        private final boolean crossed;
        private final boolean reflected;

        StepResult(int index, double oldPosition, double newPosition,
                   double oldMomentum, double newMomentum,
                   double deltaU, boolean crossed, boolean reflected) {
            this.index = index;
            this.oldPosition = oldPosition;
            this.newPosition = newPosition;
            this.oldMomentum = oldMomentum;
            this.newMomentum = newMomentum;
            this.deltaU = deltaU;
            this.crossed = crossed;
            this.reflected = reflected;
        }

        public int getIndex() {
            return index;
        }

        public double getOldPosition() {
            return oldPosition;
        }

        public double getNewPosition() {
            return newPosition;
        }

        public double getOldMomentum() {
            return oldMomentum;
        }

        public double getNewMomentum() {
            return newMomentum;
        }

        public double getDeltaU() {
            return deltaU;
        }

        public boolean isCrossed() {
            return crossed;
        }

        public boolean isReflected() {
            return reflected;
        }
    }

    public static class DetailedStepResult extends StepResult {
        private final List<SubStepEvent> events;

        DetailedStepResult(int index, double oldPosition, double newPosition,
                           double oldMomentum, double newMomentum,
                           double deltaU, boolean crossed, boolean reflected,
                           List<SubStepEvent> events) {
            super(index, oldPosition, newPosition, oldMomentum, newMomentum, deltaU, crossed, reflected);
            this.events = Collections.unmodifiableList(new ArrayList<>(events));
        }

        public List<SubStepEvent> getEvents() {
            return events;
        }

        StepResult toStepResult() {
            return new StepResult(getIndex(), getOldPosition(), getNewPosition(),
                    getOldMomentum(), getNewMomentum(), getDeltaU(), isCrossed(), isReflected());
        }
    }

    public static class SubStepEvent {
        private final String type;
        private final double startPosition;
        private final double endPosition;
        private final double boundary;
        private final double startMomentum;
        private final double endMomentum;
        private final double deltaU;
        private final double elapsedTime;

        SubStepEvent(String type, double startPosition, double endPosition, double boundary,
                     double startMomentum, double endMomentum, double deltaU, double elapsedTime) {
            this.type = type;
            this.startPosition = startPosition;
            this.endPosition = endPosition;
            this.boundary = boundary;
            this.startMomentum = startMomentum;
            this.endMomentum = endMomentum;
            this.deltaU = deltaU;
            this.elapsedTime = elapsedTime;
        }

        public String getType() {
            return type;
        }

        public double getStartPosition() {
            return startPosition;
        }

        public double getEndPosition() {
            return endPosition;
        }

        public double getBoundary() {
            return boundary;
        }

        public double getStartMomentum() {
            return startMomentum;
        }

        public double getEndMomentum() {
            return endMomentum;
        }

        public double getDeltaU() {
            return deltaU;
        }

        public double getElapsedTime() {
            return elapsedTime;
        }
    }
}
