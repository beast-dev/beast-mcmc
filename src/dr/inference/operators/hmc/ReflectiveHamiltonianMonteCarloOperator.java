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
import java.util.Arrays;
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

    /**
     * How the engine resolves walls hit at the same instant. GENERALIZED applies Gamma, the reflection across the span
     * of the normals of every wall hit at that instant, in one event; SEQUENTIAL keeps the earlier rule of one wall
     * per event (a collision wins a tie with a bound, the lowest index wins among bounds, the first pair enumerated
     * among collisions). The two differ only when two or more walls are hit at exactly the same time.
     */
    public enum ReflectionRule {
        GENERALIZED("generalized"),
        SEQUENTIAL("sequential");

        private final String name;

        ReflectionRule(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public static ReflectionRule parse(String name) {
            for (ReflectionRule rule : values()) {
                if (rule.name.equalsIgnoreCase(name)) {
                    return rule;
                }
            }
            throw new IllegalArgumentException("Unknown reflection rule '" + name + "'; use 'generalized' or 'sequential'");
        }
    }

    public static final ReflectionRule DEFAULT_REFLECTION_RULE = ReflectionRule.GENERALIZED;
    // all events in one drift, time-0 events included; a T3 caterpillar 1e-9 off its apex with 1000 nodes and strongly
    // unequal masses (log-sd 2) took 2.0 to 2.8 million positive-time events in its first drift
    public static final long DEFAULT_MAX_EVENTS_PER_STEP = 10000000L;

    private final ReflectionRule reflectionRule;
    private final long maxEventsPerStep;

    // counters, cumulative over the run
    private long eventCount = 0;            // wall events (time-0 events included)
    private long nonAdvancingEventCount = 0;
    private long wallCount = 0;             // walls reflected, summed over events (|J| per event)
    private long multiWallEventCount = 0;   // events with two or more walls (|J| >= 2)
    private int largestWallSet = 0;         // largest |J| in one event
    private long largestEventsInDrift = 0;
    private long maxEventsRejectionCount = 0;


    public ReflectiveHamiltonianMonteCarloOperator(AdaptationMode mode,
                                                   double weight,
                                                   GradientWrtParameterProvider gradientProvider,
                                                   Parameter parameter,
                                                   Transform transform,
                                                   Parameter maskParameter,
                                                   Options runtimeOptions,
                                                   MassPreconditioner preconditioner,
                                                   GeneralBoundsProvider bounds) {
        this(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner,
                bounds, DEFAULT_REFLECTION_RULE, DEFAULT_MAX_EVENTS_PER_STEP);
    }

    public ReflectiveHamiltonianMonteCarloOperator(AdaptationMode mode,
                                                   double weight,
                                                   GradientWrtParameterProvider gradientProvider,
                                                   Parameter parameter,
                                                   Transform transform,
                                                   Parameter maskParameter,
                                                   Options runtimeOptions,
                                                   MassPreconditioner preconditioner,
                                                   GeneralBoundsProvider bounds,
                                                   ReflectionRule reflectionRule,
                                                   long maxEventsPerStep) {

        super(mode, weight, gradientProvider, parameter, transform, maskParameter, runtimeOptions, preconditioner);

        if (maxEventsPerStep < 1) {
            throw new IllegalArgumentException("maxEventsPerStep must be at least 1");
        }
        this.reflectionRule = reflectionRule;
        this.maxEventsPerStep = maxEventsPerStep;

        if (maskParameter != null) {
            throw new IllegalArgumentException("Reflective HMC does not support masks");
        }

        if (bounds instanceof GraphicalParameterBound) {
            final GraphicalParameterBound graphicalBounds = (GraphicalParameterBound) bounds;
            graphicalBounds.checkParameterMapping();
            graphicalBounds.checkParameterMapping(parameter);
            if (!(preconditioner instanceof MassPreconditioner.NoPreconditioning
                    || preconditioner instanceof MassPreconditioner.DiagonalPreconditioning)) {
                for (int i = 0; i < parameter.getDimension(); i++) {
                    if (graphicalBounds.getConnectedParameterIndices(i) != null) {
                        throw new IllegalArgumentException("Reflective HMC with " + parameter.getParameterName() +
                                "'s graph bounds needs no preconditioning or a diagonal one; a dense mass matrix" +
                                " cannot reflect off the walls between connected coordinates yet");
                    }
                }
            }
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
                "zero-time event cap hits: " + zeroTimeCapCount + "\n" +
                "maxEventsPerStep rejections: " + maxEventsRejectionCount + "\n" +
                "reflection rule: " + reflectionRule.getName() + "\n" +
                "wall events: " + eventCount + " (" + nonAdvancingEventCount + " not advancing the drift, " +
                multiWallEventCount + " with two or more walls; largest wall set " + largestWallSet + ")\n" +
                "most events in one drift: " + largestEventsInDrift + "\n\n";
    }

    public ReflectionRule getReflectionRule() {
        return reflectionRule;
    }

    public long getMaxEventsPerStep() {
        return maxEventsPerStep;
    }

    public long getEventCount() {
        return eventCount;
    }

    public long getNonAdvancingEventCount() {
        return nonAdvancingEventCount;
    }

    public long getWallCount() {
        return wallCount;
    }

    public long getMultiWallEventCount() {
        return multiWallEventCount;
    }

    public int getLargestWallSet() {
        return largestWallSet;
    }

    public long getLargestEventsInDrift() {
        return largestEventsInDrift;
    }

    public long getMaxEventsRejectionCount() {
        return maxEventsRejectionCount;
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
            long eventsInDrift = 0;
            double collapsedTime = 0.0;
            while (collapsedTime < functionalStepSize) {
                ReflectionEvent event = nextEvent(position, momentum, functionalStepSize - collapsedTime, isAtBoundary);

                isAtBoundary = event.doReflection(position, momentum);
                final double advancedTime = collapsedTime + event.getEventTime();

                final int walls = event.getWallCount();
                if (walls > 0) {
                    ++eventCount;
                    wallCount += walls;
                    if (walls > 1) {
                        ++multiWallEventCount;
                    }
                    if (walls > largestWallSet) {
                        largestWallSet = walls;
                    }
                    if (++eventsInDrift > largestEventsInDrift) {
                        largestEventsInDrift = eventsInDrift;
                    }
                    if (eventsInDrift > maxEventsPerStep) {
                        ++maxEventsRejectionCount;
                        throw new NumericInstabilityException();
                    }
                }

                // An event that does not advance the drift (time zero, or too short to change collapsedTime) counts
                // toward the cap; an event that advances it resets the count
                if (advancedTime > collapsedTime) {
                    zeroTimeEvents = 0;
                } else {
                    ++nonAdvancingEventCount;
                    if (++zeroTimeEvents > ZERO_TIME_EVENT_CAP) {
                        ++zeroTimeCapCount;
                        throw new NumericInstabilityException();
                    }
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
            if (reflectionRule == ReflectionRule.GENERALIZED) {
                return nextEventGeneralized(position, momentum, intervalLength);
            }
            ReflectionEvent reflectionEventAtFixedBound = firstReflectionAtFixedBounds(position, momentum, intervalLength);
            ReflectionEvent collisionEvent = firstCollision(position, momentum, intervalLength);
            return (reflectionEventAtFixedBound.getEventTime() < collisionEvent.getEventTime()) ? reflectionEventAtFixedBound : collisionEvent;
        }


        // work space of the generalized rule
        private double[] velocity = new double[0];
        private int[] pairFirst = new int[4];
        private int[] pairSecond = new int[4];
        private int pairCount;
        private int[] boundIndex = new int[4];
        private double[] boundValue = new double[4];
        private boolean[] boundIsUpper = new boolean[4];
        private int boundCount;

        private void addPair(int i, int j) {
            if (pairCount == pairFirst.length) {
                pairFirst = Arrays.copyOf(pairFirst, 2 * pairCount);
                pairSecond = Arrays.copyOf(pairSecond, 2 * pairCount);
            }
            pairFirst[pairCount] = i;
            pairSecond[pairCount] = j;
            ++pairCount;
        }

        private void addBound(int i, double value, boolean isUpper) {
            if (boundCount == boundIndex.length) {
                boundIndex = Arrays.copyOf(boundIndex, 2 * boundCount);
                boundValue = Arrays.copyOf(boundValue, 2 * boundCount);
                boundIsUpper = Arrays.copyOf(boundIsUpper, 2 * boundCount);
            }
            boundIndex[boundCount] = i;
            boundValue[boundCount] = value;
            boundIsUpper[boundCount] = isUpper;
            ++boundCount;
        }

        /**
         * The generalized rule (06 sec. 3.2): every wall the state moves into (exit rate > 0) gets a time, 0 when the
         * state is on it or past it by rounding; J is the set of walls at the smallest time t*, by exact equality.
         * An event exactly at the end of the drift is still an event. With one wall in J the event is the same
         * Collision or Reflection, with the same arithmetic, as under the sequential rule.
         */
        private ReflectionEvent nextEventGeneralized(double[] position, WrappedVector momentum, double intervalLength) {
            final int dim = position.length;
            if (velocity.length != dim) {
                velocity = new double[dim];
            }
            for (int i = 0; i < dim; i++) {
                velocity[i] = preconditioning.getVelocity(i, momentum);
            }
            double tStar = Double.POSITIVE_INFINITY;
            pairCount = 0;
            boundCount = 0;

            for (int i = 0; i < dim; i++) {
                final int[] connected = graphicalParameterBound.getConnectedParameterIndices(i);
                if (connected == null) {
                    continue;
                }
                final double velocity1 = velocity[i];
                for (int j : connected) {
                    if (j > i) {
                        final double velocity2 = velocity[j];
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
                                continue; // NaN gap
                            }
                            if (collisionTime < tStar) {
                                tStar = collisionTime;
                                pairCount = 0;
                                boundCount = 0;
                                addPair(i, j);
                            } else if (collisionTime == tStar) {
                                addPair(i, j);
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < dim; i++) {
                final double v = velocity[i];
                if (v > 0) {
                    final double upperBound = graphicalParameterBound.getFixedUpperBound(i);
                    final double gap = upperBound - position[i];
                    final double reflectionTime;
                    if (gap > 0) {
                        reflectionTime = (upperBound - position[i]) / v;
                    } else if (gap <= 0) {
                        reflectionTime = 0.0;
                    } else {
                        continue;
                    }
                    if (reflectionTime < tStar) {
                        tStar = reflectionTime;
                        pairCount = 0;
                        boundCount = 0;
                        addBound(i, upperBound, true);
                    } else if (reflectionTime == tStar) {
                        addBound(i, upperBound, true);
                    }
                } else if (v < 0) {
                    final double lowerBound = graphicalParameterBound.getFixedLowerBound(i);
                    final double gap = position[i] - lowerBound;
                    final double reflectionTime;
                    if (gap > 0) {
                        reflectionTime = (lowerBound - position[i]) / v;
                    } else if (gap <= 0) {
                        reflectionTime = 0.0;
                    } else {
                        continue;
                    }
                    if (reflectionTime < tStar) {
                        tStar = reflectionTime;
                        pairCount = 0;
                        boundCount = 0;
                        addBound(i, lowerBound, false);
                    } else if (reflectionTime == tStar) {
                        addBound(i, lowerBound, false);
                    }
                }
            }

            if (pairCount + boundCount == 0 || tStar > intervalLength) {
                return new ReflectionEvent(ReflectionType.None, intervalLength, -1.0, new int[]{-1, -1});
            }
            if (pairCount == 1 && boundCount == 0) {
                final int i = pairFirst[0];
                return new ReflectionEvent(ReflectionType.Collision, tStar, tStar * velocity[i] + position[i],
                        new int[]{i, pairSecond[0]});
            }
            if (pairCount == 0 && boundCount == 1) {
                return new ReflectionEvent(ReflectionType.Reflection, tStar, boundValue[0], new int[]{boundIndex[0]});
            }
            return new GammaEvent(tStar, Arrays.copyOf(pairFirst, pairCount), Arrays.copyOf(pairSecond, pairCount),
                    Arrays.copyOf(boundIndex, boundCount), Arrays.copyOf(boundValue, boundCount),
                    Arrays.copyOf(boundIsUpper, boundCount));
        }

        /**
         * Gamma on two or more walls hit at the same instant (06 sec. 3.3, cluster form). The edge walls in J join
         * nodes into clusters; a cluster with a fixed-bound wall in J is pinned (every momentum in it flips), and a
         * free cluster keeps its centre-of-mass velocity and reverses every relative velocity:
         * p_i <- 2 m_i vbar - p_i, vbar = sum p / sum m. A free cluster of two nodes uses doCollision, so its arithmetic
         * is that of a single collision. Each cluster is snapped to one height (06 sec. 3.6): the bound for a pinned
         * cluster, the advanced height of the lower index for a single pair, else the mean, clamped against the
         * cluster's fixed bounds and its neighbours outside it. Needs a diagonal mass for clusters of three or more.
         */
        class GammaEvent extends ReflectionEvent {
            private final double time;
            private final int[] first;
            private final int[] second;
            private final int[] bounds;
            private final double[] values;
            private final boolean[] upper;

            GammaEvent(double time, int[] first, int[] second, int[] bounds, double[] values, boolean[] upper) {
                super(ReflectionType.Gamma, time, Double.NaN, new int[0]);
                this.time = time;
                this.first = first;
                this.second = second;
                this.bounds = bounds;
                this.values = values;
                this.upper = upper;
            }

            @Override
            int getWallCount() {
                return first.length + bounds.length;
            }

            @Override
            public boolean doReflection(double[] position, WrappedVector momentum) {
                final int dim = position.length;
                for (int k = 0; k < dim; k++) {
                    position[k] += preconditioning.getVelocity(k, momentum) * time;
                }
                ensureWorkSpace(dim);
                ++stamp;

                // clusters: union-find over the nodes of the walls in J (arrays indexed by node, valid when seen)
                for (int a = 0; a < first.length; a++) {
                    touch(first[a]);
                    touch(second[a]);
                    final int ri = root(first[a]);
                    final int rj = root(second[a]);
                    if (ri != rj) {
                        ufParent[Math.max(ri, rj)] = Math.min(ri, rj);
                    }
                }
                for (int b : bounds) {
                    touch(b);
                }

                // one pass over the touched nodes numbers the clusters and counts their members
                int clusters = 0;
                for (int s = 0; s < touchedCount; s++) {
                    final int r = root(touched[s]);
                    if (clusterStamp[r] != stamp) {
                        clusterStamp[r] = stamp;
                        clusterOf[r] = clusters;
                        clusterSize[clusters] = 0;
                        clusterPairs[clusters] = 0;
                        clusterPinned[clusters] = false;
                        clusterHasLower[clusters] = false;
                        clusterLower[clusters] = Double.NEGATIVE_INFINITY;
                        clusterUpper[clusters] = Double.POSITIVE_INFINITY;
                        ++clusters;
                    }
                    ++clusterSize[clusterOf[r]];
                }
                for (int a = 0; a < first.length; a++) {
                    final int c = clusterOf[root(first[a])];
                    ++clusterPairs[c];
                    clusterPair[c] = a;
                }
                for (int b = 0; b < bounds.length; b++) {
                    final int c = clusterOf[root(bounds[b])];
                    clusterPinned[c] = true;
                    if (upper[b]) {
                        clusterUpper[c] = Math.min(clusterUpper[c], values[b]);
                    } else {
                        clusterHasLower[c] = true;
                        clusterLower[c] = Math.max(clusterLower[c], values[b]);
                    }
                }
                // members of cluster c are memberList[clusterStart[c] .. clusterStart[c] + clusterSize[c]), sorted
                int offset = 0;
                for (int c = 0; c < clusters; c++) {
                    clusterStart[c] = offset;
                    clusterFill[c] = offset;
                    offset += clusterSize[c];
                }
                Arrays.sort(touched, 0, touchedCount);
                for (int s = 0; s < touchedCount; s++) {
                    final int node = touched[s];
                    memberList[clusterFill[clusterOf[root(node)]]++] = node;
                }

                for (int c = 0; c < clusters; c++) {
                    final int from = clusterStart[c];
                    final int to = from + clusterSize[c];
                    final boolean pinned = clusterPinned[c];
                    final boolean singlePair = !pinned && clusterSize[c] == 2 && clusterPairs[c] == 1;

                    // snap
                    double height;
                    if (pinned) {
                        height = clusterHasLower[c] ? clusterLower[c] : clusterUpper[c];
                    } else if (singlePair) {
                        height = position[first[clusterPair[c]]];
                    } else {
                        // the common height if the members already share one (always so at a time-0 event), else the mean
                        height = position[memberList[from]];
                        boolean equal = true;
                        double sum = 0.0;
                        for (int s = from; s < to; s++) {
                            sum += position[memberList[s]];
                            equal &= position[memberList[s]] == height;
                        }
                        if (!equal) {
                            height = sum / clusterSize[c];
                        }
                    }
                    double lo = Double.NEGATIVE_INFINITY;
                    double hi = Double.POSITIVE_INFINITY;
                    for (int s = from; s < to; s++) {
                        final int m = memberList[s];
                        lo = Math.max(lo, graphicalParameterBound.getFixedLowerBound(m));
                        hi = Math.min(hi, graphicalParameterBound.getFixedUpperBound(m));
                        final int[] connected = graphicalParameterBound.getConnectedParameterIndices(m);
                        if (connected != null) {
                            for (int j : connected) {
                                if (nodeStamp[j] != stamp || clusterOf[root(j)] != c) {
                                    if (graphicalParameterBound.isBelow(m, j)) {
                                        hi = Math.min(hi, position[j]);
                                    } else {
                                        lo = Math.max(lo, position[j]);
                                    }
                                }
                            }
                        }
                    }
                    if (lo <= hi) {
                        height = Math.min(Math.max(height, lo), hi);
                    }
                    for (int s = from; s < to; s++) {
                        position[memberList[s]] = height;
                    }

                    // Gamma on the cluster
                    if (pinned) {
                        for (int s = from; s < to; s++) {
                            final int m = memberList[s];
                            momentum.set(m, -momentum.get(m));
                        }
                    } else if (singlePair) {
                        final int i = first[clusterPair[c]];
                        final int j = second[clusterPair[c]];
                        final ReadableVector updated = preconditioning.doCollision(new int[]{i, j}, momentum);
                        momentum.set(i, updated.get(i));
                        momentum.set(j, updated.get(j));
                    } else {
                        double sumP = 0.0;
                        double sumM = 0.0;
                        for (int s = from; s < to; s++) {
                            final int m = memberList[s];
                            inverseMassWork[s] = inverseMassOf(m, dim);
                            sumP += momentum.get(m);
                            sumM += 1.0 / inverseMassWork[s];
                        }
                        final double vBar = sumP / sumM;
                        for (int s = from; s < to; s++) {
                            final int m = memberList[s];
                            momentum.set(m, 2.0 * vBar / inverseMassWork[s] - momentum.get(m));
                        }
                    }
                }
                return true;
            }
        }

        // work space of GammaEvent, indexed by node; an entry is valid only when its stamp equals the current one
        private int stamp = 0;
        private int[] nodeStamp = new int[0];
        private int[] ufParent = new int[0];
        private int[] touched = new int[0];
        private int touchedCount;
        private int[] clusterStamp = new int[0];
        private int[] clusterOf = new int[0];
        private int[] clusterSize = new int[0];
        private int[] clusterPairs = new int[0];
        private int[] clusterPair = new int[0];
        private boolean[] clusterPinned = new boolean[0];
        private boolean[] clusterHasLower = new boolean[0];
        private double[] clusterLower = new double[0];
        private double[] clusterUpper = new double[0];
        private int[] clusterStart = new int[0];
        private int[] clusterFill = new int[0];
        private int[] memberList = new int[0];
        private double[] inverseMassWork = new double[0];

        private void ensureWorkSpace(int dim) {
            if (nodeStamp.length != dim) {
                nodeStamp = new int[dim];
                ufParent = new int[dim];
                touched = new int[dim];
                clusterStamp = new int[dim];
                clusterOf = new int[dim];
                clusterSize = new int[dim];
                clusterPairs = new int[dim];
                clusterPair = new int[dim];
                clusterPinned = new boolean[dim];
                clusterHasLower = new boolean[dim];
                clusterLower = new double[dim];
                clusterUpper = new double[dim];
                clusterStart = new int[dim];
                clusterFill = new int[dim];
                memberList = new int[dim];
                inverseMassWork = new double[dim];
                stamp = 0;
            }
            touchedCount = 0;
        }

        private void touch(int node) {
            if (nodeStamp[node] != stamp) {
                nodeStamp[node] = stamp;
                ufParent[node] = node;
                touched[touchedCount++] = node;
            }
        }

        private int root(int node) {
            int r = node;
            while (ufParent[r] != r) {
                r = ufParent[r];
            }
            while (ufParent[node] != r) {
                final int next = ufParent[node];
                ufParent[node] = r;
                node = next;
            }
            return r;
        }

        private WrappedVector probe = null;

        private double inverseMassOf(int index, int dim) {
            if (!(preconditioning instanceof MassPreconditioner.NoPreconditioning
                    || preconditioning instanceof MassPreconditioner.DiagonalPreconditioning)) {
                throw new UnsupportedOperationException("Reflective HMC's generalized rule needs a diagonal mass" +
                        " matrix when three or more nodes are tied at once");
            }
            if (probe == null || probe.getDim() != dim) {
                probe = new WrappedVector.Raw(new double[dim]);
            }
            probe.set(index, 1.0);
            final double inverseMass = preconditioning.getVelocity(index, probe);
            probe.set(index, 0.0);
            return inverseMass;
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

        /** the number of walls this event reflects off (|J|); 0 for the end of the drift */
        int getWallCount() {
            return type == ReflectionType.None ? 0 : 1;
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
        // Gamma on several walls at once; WithGraphBounds.GammaEvent does the reflection itself
        Gamma {
            @Override
            void doReflection(double[] position, MassPreconditioner preconditioning, WrappedVector momentum,
                              double eventLocation[], int[] indices, double[] normalVector, double time,
                              double remainingTime) {
                throw new IllegalStateException("Gamma events reflect through GammaEvent.doReflection");
            }
        },

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
