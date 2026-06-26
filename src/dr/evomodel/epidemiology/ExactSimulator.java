package dr.evomodel.epidemiology;

import dr.math.MathUtils;

public class ExactSimulator extends StochasticSimulator {

    public ExactSimulator(CompartmentalModel compartmentalModel) {
        super(compartmentalModel);
    }

    public void simulateTrajectory() {

        // set up time interval vector
        // duration for which we need to simulate trajectory
        //double T = compartmentalModel.origin.getParameterValue(0);
        double T = compartmentalModel.getOldestOrigin();

        // next index of compartmentalModel compartmentCounts parameter that needs to be set
        // start with last index, furthest into past and proceed until we reach index 0
        // which corresponds to most recent sampling time
        int nextRecordIndex = numGridPoints-1;
        // set default compartment counts for time intervals that completely precede origin
        // each interval has width cutOff/numGridPoints
        // last interval goes from ((numGridPoints-1)/numGridPoints)*cutOff to cutOff
        while (nextRecordIndex * intervalWidth >= T) {
            compartmentalModel.setDefaultCompartmentCounts(nextRecordIndex);
            nextRecordIndex--;
        }

        // set initial compartment counts for time interval that contains origin
        compartmentalModel.setOriginTimeCompartmentCounts(nextRecordIndex);

        // Initialize time for forward time stochastic simulation. Start at 0.0 simulate for total time of T.
        // simulationTime = 0.0 corresponds to time of origin
        // model time is "backward time" that increases into past, but simulation time is "forward time"
        double simulationTime = 0.0;
        // start time (in forward time) of next interval that needs to have compartment counts set
        // index of this interval will correspond to nextRecordIndex
        // set compartment counts for this interval to whatever simulated values are at nextIntervalStartTime
        double nextIntervalStartTime = T - nextRecordIndex*intervalWidth;
        // from now on, increase nextIntervalStartTime by simply adding intervalWidth

        // keep track of current compartment counts (needed for simulation)
        double[] currentCounts = new double[numSpecies];

        for (int s = 0; s < numSpecies; s++) {
            currentCounts[s] = compartmentalModel.compartmentCounts.get(s).getParameterValue(nextRecordIndex);
        }
        nextRecordIndex--;
        double[] reactionInt = compartmentalModel.getReactionIntensities(currentCounts);

        double timeToReaction;
        int sampledReactionChannel;
        double r0;

        while (nextRecordIndex >= 0) {
            if (nextRecordIndex < 0) {
                break;
            }

            // perform one step of Gillespie stochastic simulation algorithm
            r0 = 0;
            for (int c = 0; c < numReactionChannels; c++) {
                r0 = r0 + reactionInt[c];
            }

            //find time to next reaction
            timeToReaction = -Math.log(MathUtils.nextDouble()) / r0;

            // if next reaction occurs after nextIntervalStartTime, record current compartment counts for next interval
            while ((simulationTime + timeToReaction > nextIntervalStartTime) && nextRecordIndex >= 0) {
                for (int s = 0; s < numSpecies; s++) {
                    compartmentalModel.compartmentCounts.get(s).setParameterValue(nextRecordIndex, currentCounts[s]);
                }
                nextRecordIndex--;
                nextIntervalStartTime = nextIntervalStartTime + intervalWidth;
            }

            sampledReactionChannel = sampleReactionChannel(reactionInt, r0);

            //update simulationTime, current compartment counts and reaction intensities
            simulationTime = simulationTime + timeToReaction;

            currentCounts = compartmentalModel.introduceSecondPathogen(simulationTime, currentCounts);

            for (int s = 0; s < numSpecies; s++) {
                currentCounts[s] = currentCounts[s] + vMatrix[s][sampledReactionChannel];
            }
            reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
        }
    }
}
