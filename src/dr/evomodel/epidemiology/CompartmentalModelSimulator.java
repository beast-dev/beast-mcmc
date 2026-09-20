package dr.evomodel.epidemiology;

import dr.math.MathUtils;

public abstract class CompartmentalModelSimulator {

    public CompartmentalModel compartmentalModel;
    protected double cutOff;
    protected int numGridPoints;
    protected int numSpecies;
    protected int numReactionChannels;
    protected int[][] vMatrix;
    protected double intervalWidth;

    // boolean to keep track of lineage count constraint
    protected boolean lineageConstraintViolated = false;
    // index (i,j) of array corresponds to (pathogen, shared trajectory index)
    // shared trajectory index runs from 0 to numGridPoints-1
    // Will be null if there is no constraint to check
    protected int[][] lineageCounts = null;

    public CompartmentalModelSimulator(CompartmentalModel compartmentalModel) {
        this.compartmentalModel = compartmentalModel;
        this.cutOff = compartmentalModel.cutOff;
        this.numGridPoints = compartmentalModel.numGridPoints;
        this.numSpecies = compartmentalModel.numSpecies;
        this.numReactionChannels = compartmentalModel.numReactionChannels;
        this.vMatrix = compartmentalModel.vMatrix;
        this.intervalWidth = cutOff/numGridPoints;
    }

    public abstract void simulateTrajectory();

    //protected abstract double[] getPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau, double simTime);

    protected SimulationState initializeSimulation() {

        // make sure this is reset at start of each simulation
        lineageConstraintViolated = false;
        int nextRecordIndex = numGridPoints-1;
        double oldestOrigin = compartmentalModel.getOldestOrigin();
        //System.out.println("oldestOrigin from initializeSimulation(): " + oldestOrigin);

        // set default compartment counts for time intervals that completely precede origin
        while (nextRecordIndex * intervalWidth > oldestOrigin) {
            compartmentalModel.setDefaultCompartmentCounts(nextRecordIndex);
            nextRecordIndex--;
        }
        // set initial compartment counts for time interval that contains origin
        compartmentalModel.setOriginTimeCompartmentCounts(nextRecordIndex);

        SimulationState state = new SimulationState();
        // Initialize time for forward time stochastic simulation.
        // simulation starts at time of oldest origin and continues until more recent of the most recent sampling dates.
        // model time is "backward time" that increases into past, but simulation time is "forward time"
        // Time = 0.0 corresponds to time of cutOff
        state.simulationTime = cutOff - oldestOrigin;
        // keep track of current compartment counts (needed for simulation)
        state.currentCounts = new double[numSpecies];

        for (int s = 0; s < numSpecies; s++) {
            state.currentCounts[s] = compartmentalModel.compartmentCounts.get(s).getParameterValue(nextRecordIndex);
        }
        nextRecordIndex--;
        state.nextRecordIndex = nextRecordIndex;
        // start time (in forward time) of next interval that needs to have compartment counts set
        // index of this interval will correspond to nextRecordIndex
        // set compartment counts for this interval to whatever simulated values are at nextIntervalStartTime
        state.nextIntervalStartTime = cutOff-nextRecordIndex*intervalWidth;
        // from now on, increase nextIntervalStartTime by simply adding intervalWidth

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
        //System.out.println("recordCompartmentCountsUpTo: candidateTime=" + candidateTime
        //        + " nextIntervalStartTime=" + state.nextIntervalStartTime
        //        + " nextRecordIndex=" + state.nextRecordIndex
        //        + " SS=" + state.currentCounts[0]);

        while (candidateTime > state.nextIntervalStartTime && state.nextRecordIndex >= 0) {
            //System.out.println("Recording grid point " + state.nextRecordIndex +
            //        " at simulationTime=" + state.simulationTime +
            //        " candidateTime=" + candidateTime +
            //        " nextIntervalStartTime=" + state.nextIntervalStartTime +
            //        " IS=" + state.currentCounts[4] +
            //        " SI=" + state.currentCounts[1]);


            for (int s = 0; s < numSpecies; s++) {
                compartmentalModel.compartmentCounts.get(s).setParameterValue(state.nextRecordIndex, state.currentCounts[s]);
            }
            //System.out.println("in recordCompartmentCountsUpTo");
            if(!checkLineageCountConstraint(state.nextRecordIndex, state.currentCounts)) {
                lineageConstraintViolated = true;
                //System.out.println("Lineage constraint violated. simTime=" + state.simulationTime + " IS=" + state.currentCounts[4]);
                // set nextRecordIndex to -1 so that the while loop in simulateTrajectory() will be exited
                state.nextRecordIndex = -1;
                return;
            }

            state.nextRecordIndex--;
            state.nextIntervalStartTime = state.nextIntervalStartTime + intervalWidth;
        }
    }

    protected static class SimulationState {

        // Time for forward time stochastic simulation. Will start at time of oldest origin and simulate for total time
        // equivalent to time between oldest origin and more recent of the most recent sampling dates.
        // simulationTime = 0.0 corresponds to time of cutOff.
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

    public void setLineageCounts(int[][] lineageCounts) {
        this.lineageCounts = lineageCounts;
    }

    public boolean isLineageConstraintViolated() {
        return lineageConstraintViolated;
    }

    // Returns true if the lineage count constraint is satisfied (if the number of infected individuals
    // for a given pathogen trajectory is >= the lineage count of the corresponding pathogen tree)
    private boolean checkLineageCountConstraint(int index, double[] currentCounts){
        if(lineageCounts == null){
            return true;
        }
        int[] numInfected = compartmentalModel.getLineageCountConstraintCounts(currentCounts);
        for(int i = 0; i < numInfected.length; i++){
            if(numInfected[i] < lineageCounts[i][index]){
               // System.out.println("Constraint violated at index " + index +
               //         " pathogen " + i +
               //         " infectedCount=" + numInfected[i] +
               //         " lineageCount=" + lineageCounts[i][index]);
                // print recorded infected counts for pathogen 1 at all grid points so far
               // System.out.println("Recorded IS counts at grid points:");
               // for (int k = index; k <= 5; k++) {
               //     System.out.println("  index " + k + ": IS=" +
                //            compartmentalModel.compartmentCounts.get(4).getParameterValue(k));
                //}
                return false;
            }
        }
        //System.out.println("index: " + index);
        //System.out.println("currentCounts: " + Arrays.toString(currentCounts));
        return true;
    }

    public void resetLineageConstraintViolated() {
        lineageConstraintViolated = false;
    }
}
