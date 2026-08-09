package dr.evomodel.epidemiology;

import dr.math.MathUtils;

public abstract class StochasticSimulator {

    protected CompartmentalModel compartmentalModel;
    protected double cutOff;
    protected int numGridPoints;
    protected int numSpecies;
    protected int numReactionChannels;
    protected int[][] vMatrix;
    protected double intervalWidth;

    public StochasticSimulator(CompartmentalModel compartmentalModel) {
        this.compartmentalModel = compartmentalModel;
        this.cutOff = compartmentalModel.cutOff;
        this.numGridPoints = compartmentalModel.numGridPoints;
        this.numSpecies = compartmentalModel.numSpecies;
        this.numReactionChannels = compartmentalModel.numReactionChannels;
        this.vMatrix = compartmentalModel.vMatrix;
        this.intervalWidth = cutOff/numGridPoints;
    }

    public abstract void simulateTrajectory();

    protected SimulationState initializeSimulation(double T) {

        int nextRecordIndex = numGridPoints-1;

        // set default compartment counts for time intervals that completely precede origin
        while (nextRecordIndex * intervalWidth >= T) {
            compartmentalModel.setDefaultCompartmentCounts(nextRecordIndex);
            nextRecordIndex--;
        }

        // set initial compartment counts for time interval that contains origin
        compartmentalModel.setOriginTimeCompartmentCounts(nextRecordIndex);

        SimulationState state = new SimulationState();

        // Initialize time for forward time stochastic simulation. Start at 0.0 simulate for total time of T.
        // simulationTime = 0.0 corresponds to time of origin
        // model time is "backward time" that increases into past, but simulation time is "forward time"
        state.simulationTime = 0.0;

        // start time (in forward time) of next interval that needs to have compartment counts set
        // index of this interval will correspond to nextRecordIndex
        // set compartment counts for this interval to whatever simulated values are at nextIntervalStartTime
        state.nextIntervalStartTime = T-nextRecordIndex*intervalWidth;
        // from now on, increase nextIntervalStartTime by simply adding intervalWidth

        // keep track of current compartment counts (needed for simulation)
        state.currentCounts = new double[numSpecies];

        for (int s = 0; s < numSpecies; s++) {
            state.currentCounts[s] = compartmentalModel.compartmentCounts.get(s).getParameterValue(nextRecordIndex);
        }
        nextRecordIndex--;
        state.nextRecordIndex = nextRecordIndex;

        return state;
    }

    protected int sampleReactionChannel(double[] reactionInt, double reactionIntSum) {
        double r = MathUtils.nextDouble();
        double threshold = r*reactionIntSum;
        double cumulative = 0.0;

        for (int i = 0; i < reactionInt.length; i++) {
            cumulative = cumulative + reactionInt[i];
            if (threshold < cumulative) {
                return i;
            }
        }
        return reactionInt.length - 1;
    }

    protected double sumIntensities(double[] reactionInt) {
        double r0 = 0;
        for (int c = 0; c < numReactionChannels; c++) {
            r0 = r0 + reactionInt[c];
        }
        return r0;
    }

    protected void recordCompartmentCountsUpTo(SimulationState state, double candidateTime) {
        while (candidateTime > state.nextIntervalStartTime && state.nextRecordIndex >= 0) {
            for (int s = 0; s < numSpecies; s++) {
                compartmentalModel.compartmentCounts.get(s).setParameterValue(state.nextRecordIndex, state.currentCounts[s]);
            }
            state.nextRecordIndex--;
            state.nextIntervalStartTime = state.nextIntervalStartTime + intervalWidth;
        }
    }

    protected static class SimulationState {

        // Time for forward time stochastic simulation. Will start at 0.0 simulate for total time of T.
        // simulationTime = 0.0 corresponds to time of origin
        // model time is "backward time" that increases into past, but simulation time is "forward time"
        double simulationTime;

        // next index of compartmentalModel compartmentCounts parameter that needs to be set
        // start with last index, furthest into past and proceed until we reach index 0
        // which corresponds to most recent sampling time
        int nextRecordIndex;

        // start time (in forward time) of next interval that needs to have compartment counts set
        // index of this interval will correspond to nextRecordIndex
        // set compartment counts for this interval to whatever simulated values are at nextIntervalStartTime
        double nextIntervalStartTime;

        // keep track of current compartment counts (needed for simulation)
        double[] currentCounts;

        // reaction intensities
        double[] reactionInt;

        // Maximum number of times that a reaction with a positive intensity can fire before
        // exhausting one of its reactants
        double[] maxFiringTimes;
    }
}
