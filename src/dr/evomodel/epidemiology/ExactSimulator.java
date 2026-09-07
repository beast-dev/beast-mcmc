package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class ExactSimulator extends StochasticSimulator {

    protected Parameter elapsedTime;

    public ExactSimulator(CompartmentalModel compartmentalModel) {
        super(compartmentalModel);
    }

    public ExactSimulator(CompartmentalModel compartmentalModel, Parameter elapsedTime) {
        super(compartmentalModel);
        this.elapsedTime = elapsedTime;
    }

    public void simulateTrajectory() {

        long startTime = System.nanoTime();

        //double T = compartmentalModel.getOldestOrigin();

        SimulationState state = initializeSimulation();

        state.reactionInt = compartmentalModel.getReactionIntensities(state.currentCounts, state.simulationTime);

        while (state.nextRecordIndex >= 0) {
            runOneSSAStep(state);
        }

        long endTime = System.nanoTime();
        if(elapsedTime != null) {
            double elapsedTimeInSeconds = (endTime - startTime) / 1e9;
            System.out.println("Elapsed time: " + elapsedTimeInSeconds + " seconds");
            elapsedTime.setParameterValue(0, elapsedTimeInSeconds);
        }
    }

    private void runOneSSAStep(SimulationState state) {

        double r0 = sumIntensities(state.reactionInt);

        // find time to next reaction
        double timeToReaction = -Math.log(MathUtils.nextDouble()) / r0;

        // if next reaction occurs after nextIntervalStartTime, record current compartment counts for next interval
        recordCompartmentCountsUpTo(state, state.simulationTime + timeToReaction);

        int sampledReactionChannel = sampleReactionChannel(state.reactionInt, r0);

        // update simulationTime, current compartment counts and reaction intensities
        state.simulationTime = state.simulationTime + timeToReaction;
        state.currentCounts = compartmentalModel.introduceSecondPathogen(state.simulationTime, state.currentCounts);

        for (int s = 0; s < numSpecies; s++) {
            state.currentCounts[s] = state.currentCounts[s] + vMatrix[s][sampledReactionChannel];
        }
        state.reactionInt = compartmentalModel.getReactionIntensities(state.currentCounts, state.simulationTime);
    }

}