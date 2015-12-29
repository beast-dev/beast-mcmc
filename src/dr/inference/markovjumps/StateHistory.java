/*
 * StateHistory.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.markovjumps;

import dr.evolution.datatype.DataType;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A class to represent the complete state history of a continuous-time Markov chain in the
 * interval [0,T].
 * <p/>
 * This work is supported by NSF grant 0856099
 *
 * @author Marc A. Suchard
 */

public class StateHistory {

//    private StateHistory(int startingState, int stateCount) {
//        this(0.0, startingState, stateCount);
//    }

    protected StateHistory(double startingTime, int startingState, int stateCount) {
        stateList = new ArrayList<StateChange>();
        stateList.add(new StateChange(startingTime, startingState));
        this.stateCount = stateCount;
        finalized = false;
    }

    public void addChange(StateChange stateChange) {
        checkFinalized(false);
        stateList.add(stateChange);
    }

    public void addEndingState(StateChange stateChange) {
        checkFinalized(false);
        stateList.add(stateChange);
        finalized = true;
    }

    public int[] getJumpCounts() {
        int[] counts = new int[stateCount * stateCount];
        accumulateSufficientStatistics(counts, null);
        return counts;
    }

    public double[] getWaitingTimes() {
        double[] times = new double[stateCount];
        accumulateSufficientStatistics(null, times);
        return times;
    }

    public double getTotalRegisteredCounts(double[] register) {
        int[] counts = getJumpCounts();
//        double total = 0;
//        for (int i = 0; i < counts.length; i++) {
//            total += counts[i] * register[i];
//        }
//        return total;
        return dotProduct(counts, register);
    }

    private double dotProduct(int[] a, double[] b) {
        double total = 0;
        final int length = a.length;
        for (int i = 0; i < length; i++) {
            total += a[i] * b[i];
        }
        return total;
    }

    public double getTotalReward(double[] register) {
        double[] times = getWaitingTimes();
        double total = 0;
        for (int i = 0; i < times.length; i++) {
            total += times[i] * register[i]; // stateCount length vector
        }
        return total;
    }

    public void accumulateSufficientStatistics(int[] counts, double[] times) {
        checkFinalized(true);
        int nJumps = getNumberOfJumps();

        StateChange initialState = stateList.get(0);
        int currentState = initialState.getState();
        double currentTime = initialState.getTime();

        for (int i = 1; i <= nJumps; i++) {

            StateChange nextStateChange = stateList.get(i);
            int nextState = nextStateChange.getState();
            double nextTime = nextStateChange.getTime();

            if (counts != null) {
                counts[currentState * stateCount + nextState]++;
            }
            if (times != null) {
                times[currentState] += (nextTime - currentTime);
            }
            currentState = nextState;
            currentTime = nextTime;
        }
        if (times != null) { // Add last waiting time
            StateChange finalState = stateList.get(nJumps + 1);
            times[currentState] += (finalState.getTime() - currentTime);
        }
    }

    public int getNumberOfJumps() {
        checkFinalized(true);
        return stateList.size() - 2; // Discount starting and ending states
    }

    private void checkFinalized(boolean isTrue) {
        if (isTrue != finalized) {
            throw new IllegalAccessError("StateHistory " + (finalized ? "is" : "is not" + " finalized"));
        }
    }

    public int getStartingState() {
        return stateList.get(0).getState();
    }

    public int getEndingState() {
        checkFinalized(true);
        return stateList.get(stateList.size() - 1).getState();
    }

    public double getStartingTime() {
        return stateList.get(0).getTime();
    }

    public double getEndingTime() {
        checkFinalized(true);
        return stateList.get(stateList.size() - 1).getTime();
    }

    public void rescaleTimesOfEvents(double inStartTime, double inEndTime) {

        final double scale = (inEndTime - inStartTime) / (getEndingTime() - getStartingTime());

        StateChange currentStateChange = stateList.get(0);
        double oldCurrentTime = currentStateChange.getTime();
        currentStateChange.setTime(inStartTime);
        double newCurrentTime = inStartTime;

        for (int i = 1; i < stateList.size(); ++i) {
            StateChange nextStateChange = stateList.get(i);
            double oldNextTime = nextStateChange.getTime();
            double oldTimeDiff = oldNextTime - oldCurrentTime;

            double newNextTime = oldTimeDiff * scale + newCurrentTime;
            nextStateChange.setTime(newNextTime);

            oldCurrentTime = oldNextTime;
            newCurrentTime = newNextTime;

        }
    }

    public String toString() {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < stateList.size(); i++) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(stateList.get(i));
        }
        sb.append("]");
        return sb.toString();
    }

    public static void main(String[] args) {
        System.err.println("Testing time rescaling:");
        StateHistory stateHistory = new StateHistory(1, 1, 4);
        StateChange stateChange;
        stateChange = new StateChange(2, 2);
        stateHistory.addChange(stateChange);
        stateChange = new StateChange(5, 2);
        stateHistory.addEndingState(stateChange);

        System.err.println("Initial history: " + stateHistory);

        stateHistory.rescaleTimesOfEvents(8.0, 0.0);
        System.err.println("Rescale history: " + stateHistory);

        stateHistory.rescaleTimesOfEvents(0.0, 4.0);
        System.err.println("Rescale history: " + stateHistory);
    }

    public StateHistory filterChanges(double[] register) {

        if (getNumberOfJumps() == 0) {
            return this;
        }

        StateChange currentState = stateList.get(0);
        StateHistory newHistory = new StateHistory(currentState.getTime(), currentState.getState(), stateCount);

        for (int i = 1; i < stateList.size() - 1; ++i) {
            StateChange nextState = stateList.get(i);
            if (register[currentState.getState() * stateCount + nextState.getState()] == 1) {
                nextState = nextState.clone();
                nextState.setPreviousState(currentState.getState());
                newHistory.addChange(nextState);
            }
            currentState = nextState;
        }
        newHistory.addEndingState(stateList.get(stateList.size() - 1));
        // This function can produce inconsistent histories when not all changes are reported.
        isFiltered = true;
        return newHistory;
    }

    public double getLogLikelihood(final double[] infinitesimalRates, final int stateCount) {
        checkFinalized(true);

        // TODO This function needs testing
        double logLikelihood = 0.0;
        final int totalChanges = getNumberOfJumps();

        int currentState = stateList.get(0).getState();
        double currentTime = stateList.get(0).getTime();
        for (int i = 1; i < totalChanges; ++i) {
            int nextState = stateList.get(i).getState();
            double nextTime = stateList.get(i).getTime();

            // Exponential pdf and destination choice
            logLikelihood += Math.log(infinitesimalRates[currentState * stateCount + nextState])
                    + infinitesimalRates[currentState * stateCount + currentState] * (nextTime - currentTime);
            // terms involving Math.log(\lambda_{ii}) cancel

            currentState = nextState;
            currentTime = nextTime;
        }

        final int lastState = stateList.get(stateList.size() - 1).getState();
        final double lastTime = stateList.get(stateList.size() - 1).getTime();

        assert (lastState == currentState);
        assert (lastTime >= currentTime);

        // No event in last interval
        logLikelihood += infinitesimalRates[currentState * stateCount + currentState] * (lastTime - currentTime);

        return logLikelihood;
    }

    public String toStringChanges(int site, DataType dataType) {
        return toStringChanges(site, dataType, true);
    }

    public String toStringChanges(int site, DataType dataType, boolean wrap) {
        StringBuilder sb = wrap ? new StringBuilder("{") : new StringBuilder();
        // site number gets put into each and every event string
//        sb.append(site).append(",");
        int currentState = stateList.get(0).getState();
        boolean firstChange = true;
        for (int i = 1; i < stateList.size() - 1; i++) {  // TODO Code review: should this really be size() - 1?
            int nextState = stateList.get(i).getState();
            if (isFiltered) {
                currentState = stateList.get(i).getPreviousState();
            }
            if (nextState != currentState) {
                if (!firstChange) {
                    sb.append(",");
                }
                double time = stateList.get(i).getTime(); // + startTime;
                addEventToStringBuilder(sb, dataType.getCode(currentState), dataType.getCode(nextState), time, site);
                firstChange = false;
                currentState = nextState;
            }
        }
        if (wrap) {
            sb.append("}"); // Always returns an array of arrays
        }
        return sb.toString();
    }

    public static void addEventToStringBuilder(StringBuilder sb, String source, String dest, double time, int site) {
        // AR changed this to match an attribute array:
        sb.append("{");
        if (site > 0) {
            sb.append(site).append(",");
        }
        sb.append(time).append(",").append(source).append(",").append(dest).append("}");
    }

    public static StateHistory simulateConditionalOnEndingState(double startingTime,
                                                                int startingState,
                                                                double endingTime,
                                                                int endingState,
                                                                double[] lambda,
                                                                int stateCount) {
        throw new RuntimeException("Impossible to simulate a conditioned CTMC in StateHistory");
    }

    public static StateHistory simulateUnconditionalOnEndingState(double startingTime,
                                                                  int startingState,
                                                                  double endingTime,
                                                                  double[] lambda,
                                                                  int stateCount) {

        StateHistory history = new StateHistory(startingTime, startingState, stateCount);
        double[] multinomial = new double[stateCount];

        double currentTime = startingTime;
        int currentState = startingState;

        while (currentTime < endingTime) {

            double currentRate = -lambda[currentState * stateCount + currentState];
            double waitingTime = MathUtils.nextExponential(currentRate);

            currentTime += waitingTime;
            if (currentTime < endingTime) { // Simulate a jump
                System.arraycopy(lambda, currentState * stateCount, multinomial, 0, stateCount);
                multinomial[currentState] = 0;
                currentState = MathUtils.randomChoicePDF(multinomial); // Does not need to be normalized

                history.addChange(new StateChange(currentTime, currentState));
            }
        }

        history.addEndingState(new StateChange(endingTime, currentState));

        return history;
    }


    private int stateCount;
    private List<StateChange> stateList;
    private boolean finalized;
    private boolean isFiltered = false;

}
