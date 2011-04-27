package dr.inference.markovjumps;

import dr.evolution.datatype.DataType;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A class to represent the complete state history of a continuous-time Markov chain in the
 * interval [0,T].
 *
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
            times[currentState] += ( finalState.getTime() - currentTime);
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
        return stateList.get(stateList.size()-1).getState();
    }

    public double getStartingTime() {
        return stateList.get(0).getTime();
    }

    public double getEndingTime() {
        checkFinalized(true);
        return stateList.get(stateList.size()-1).getState();
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

    public String toStringChanges(DataType dataType, double startTime) {
        StringBuilder sb = new StringBuilder("{");
        int currentState = stateList.get(0).getState();
        boolean firstChange = true;
        for (int i = 1; i < stateList.size(); i++) {
            int nextState = stateList.get(i).getState();
            if (nextState != currentState) {
                if (!firstChange) {
                    sb.append(",");
                }
                double time = stateList.get(i).getTime() + startTime;
                addEventToStringBuilder(sb, dataType.getCode(currentState), dataType.getCode(nextState), time);
                firstChange = false;
                currentState = nextState;
            }
        }
        sb.append("}");
        return sb.toString();
    }

    public static void addEventToStringBuilder(StringBuilder sb, String source, String dest, double time) {
        sb.append("[").append(source).append(":").append(dest).append(":").append(time).append("]");
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
}
