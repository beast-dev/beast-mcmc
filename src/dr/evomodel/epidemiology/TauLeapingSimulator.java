package dr.evomodel.epidemiology;

import cern.jet.random.Poisson;
import dr.math.MathUtils;
import java.util.ArrayList;
import java.util.List;

public class TauLeapingSimulator extends StochasticSimulator {

    protected double epsilon;
    protected int criticalNumber;

    public TauLeapingSimulator(CompartmentalModel compartmentalModel, double epsilon, int criticalNumber) {
        super(compartmentalModel);
        this.epsilon = epsilon;
        this.criticalNumber = criticalNumber;
    }

    // Implements hybrid tau-leaping/SSA algorithm with step size selection as outlined by Cao et al. (2006)
    public void simulateTrajectory() {

        // set up time interval vector
        // duration for which we need to simulate trajectory
        //double T = compartmentalModel.originOne.getParameterValue(0);
        double T = compartmentalModel.getOldestOrigin();

        // next index of compartmentalModel compartmentCounts parameter that needs to be set
        // start with last index, furthest into past and proceed until we reach index 0
        // which corresponds to most recent sampling time
        int nextRecordIndex = numGridPoints - 1;
        // set default compartment counts for time intervals that completely precede origin
        // each interval has width of cutOff/numGridPoints
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
        double nextIntervalStartTime = T - nextRecordIndex * intervalWidth;
        // from now on, increase nextIntervalStartTime by simply adding intervalWidth

        // keep track of current compartment counts (needed for simulation)
        double[] currentCounts = new double[numSpecies];

        for (int s = 0; s < numSpecies; s++) {
            currentCounts[s] = compartmentalModel.compartmentCounts.get(s).getParameterValue(nextRecordIndex);
        }
        nextRecordIndex--;

        List<Integer> critical = new ArrayList<>(numReactionChannels);
        List<Integer> noncritical = new ArrayList<>(numReactionChannels);

        // Simulate until we have set compartment counts for all intervals
        // instead of while(simulationTime < T), keep track of nextRecordIndex
        while (nextRecordIndex >= 0) {

            int[] g = compartmentalModel.getHighestOrdersOfReactions();
            double[] e = new double[numSpecies];
            for (int s = 0; s < numSpecies; s++) {
                e[s] = epsilon/g[s];
            }

            // Step 1
            double[] reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
            // Maximum number of times that a reaction with a positive intensity can fire before
            // exhausting one of its reactants
            double[] maxFiringTimes = getMaxFiringTimes(currentCounts, reactionInt);

            for (int c = 0; c < numReactionChannels; c++) {
                if (maxFiringTimes[c] < criticalNumber) {
                    critical.add(c);
                } else {
                    noncritical.add(c);
                }
            }

            // Step 2
            double tauPrime = computeTauPrime(noncritical, reactionInt, e, currentCounts);

            boolean tauPrimeTooSmall = true;

            while (tauPrimeTooSmall && nextRecordIndex >= 0) {

                // step 3
                double r0 = 0;
                for (int c = 0; c < numReactionChannels; c++) {
                    r0 = r0 + reactionInt[c];
                }

                if (tauPrime < 10.0 / r0) {
                    // abandon tau-leaping temporarily and execute
                    // 100 steps of the Gillespie single reaction stochastic simulation algorithm
                    for (int k = 0; k < 100; k++) {

                        if (nextRecordIndex < 0) {
                            break;
                        }

                        // perform one step of Gillespie stochastic simulation algorithm

                        //find time to next reaction
                        double timeToReaction = -Math.log(MathUtils.nextDouble()) / r0;

                        // if next reaction occurs after nextIntervalStartTime, record current compartment counts for next interval
                        while ((simulationTime + timeToReaction > nextIntervalStartTime) && nextRecordIndex >= 0) {
                            for (int s = 0; s < numSpecies; s++) {
                                compartmentalModel.compartmentCounts.get(s).setParameterValue(nextRecordIndex, currentCounts[s]);
                            }
                            nextRecordIndex--;
                            nextIntervalStartTime = nextIntervalStartTime + intervalWidth;
                        }

                        int sampledReactionChannel = sampleReactionChannel(reactionInt, r0);

                        //update simulationTime and current compartment counts
                        //double previousSimulationTime = simulationTime;
                        simulationTime = simulationTime + timeToReaction;

                        currentCounts = compartmentalModel.introduceSecondPathogen(
                                //previousSimulationTime,
                                simulationTime,
                                currentCounts);

                        for (int s = 0; s < numSpecies; s++) {
                            currentCounts[s] = currentCounts[s] + vMatrix[s][sampledReactionChannel];
                        }
                        reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
                        r0 = 0;
                        for (int c = 0; c < numReactionChannels; c++) {
                            r0 = r0 + reactionInt[c];
                        }
                    }
                    // end of 100 SSA steps

                    // Step 1 (again)
                    // Step 1 has been removed from beginning of loop so that we do not have to repeat

                    maxFiringTimes = getMaxFiringTimes(currentCounts, reactionInt);
                    critical.clear();
                    noncritical.clear();

                    for (int c = 0; c < numReactionChannels; c++) {
                        if (maxFiringTimes[c] < criticalNumber) {
                            critical.add(c);
                        } else {
                            noncritical.add(c);
                        }
                    }

                    // Step 2 (again)
                    tauPrime = computeTauPrime(noncritical, reactionInt, e, currentCounts);

                } else {
                    // tauPrime is not too small, and we can proceed to Step 4
                    tauPrimeTooSmall = false;
                }
            }


            if (nextRecordIndex >= 0) {

                // checks for negative counts or counts below a minimal threshold for certain species
                boolean hasMinimalCounts = false;

                // store candidate updated counts and candidate tau
                // will only be used if all counts meet minimal threshold
                double[] updatedCounts = new double[numSpecies];
                double tau = 0;

                while (!hasMinimalCounts && nextRecordIndex >= 0) {

                    // step 4

                    // calculate a second option for tau
                    double r0Critical = 0.0;
                    for (int index = 0; index < critical.size(); index++) {
                        r0Critical += reactionInt[critical.get(index)];
                    }

                    double tauDoublePrime;

                    if (r0Critical > 0.0) {
                        double u = MathUtils.nextDouble();
                        tauDoublePrime = -Math.log(u) / r0Critical;
                    } else {
                        tauDoublePrime = Double.POSITIVE_INFINITY;
                    }

                    // step 5
                    // determine tau
                    tau = Math.min(tauPrime, tauDoublePrime);

                    reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
                    double[] poissonIntensities = getPoissonIntensities(currentCounts, reactionInt, tau);

                    double[] numFirings = new double[numReactionChannels];

                    if (tauPrime < tauDoublePrime) {

                        // for all critical reactions, no reaction occurs (Step 5(a))
                        // for all noncritical reactions, generate new number of "firings" as a Poisson sample.
                        for (int j = 0; j < numReactionChannels; j++) {
                            if (maxFiringTimes[j] < criticalNumber) {
                                numFirings[j] = 0;
                            } else {
                                numFirings[j] = Poisson.staticNextInt(poissonIntensities[j]);
                            }
                        }

                    } else {
                        // Step 5(b), one critical reaction occurs, other critical reactions are 0,
                        // noncritical reactions are poisson sample
                        // the critical reaction that does occur is sampled same way as SSA

                        // first sample the noncritical reactions
                        for (int j = 0; j < numReactionChannels; j++) {
                            if (maxFiringTimes[j] < criticalNumber) {
                                numFirings[j] = 0;
                            } else {
                                numFirings[j] = Poisson.staticNextInt(poissonIntensities[j]);
                            }
                        }

                        // now choose jc among critical reactions only
                        if (r0Critical > 0.0) {
                            // draw uniform sample from (0, r0Critical)
                            double u = MathUtils.nextDouble() * r0Critical;
                            double cumulative = 0.0;
                            // jc is the index of the reaction that occurs, start at a value that it can't be
                            int jc = -1;
                            for (int idx = 0; idx < critical.size(); idx++) {
                                int j = critical.get(idx);
                                cumulative += reactionInt[j];
                                if (u < cumulative) {
                                    jc = j;
                                    break;
                                }
                            }
                            // exactly one critical reaction fires
                            if (jc >= 0) {
                                numFirings[jc] = 1;
                            }
                        }
                    }

                    updatedCounts = compartmentalModel.getUpdatedCompartmentCounts(currentCounts, numFirings);

                    // Do we need to check for anything beyond negative counts?
                    hasMinimalCounts = compartmentalModel.hasMinimalCounts(updatedCounts);

                    // Step 6
                    // If there are negative counts, reduce tauPrime by half, run 100 steps of SSA
                    // and then go through Steps 4 - 6 again (this will be done because the loop condition !hasMinimalCounts is true)
                    if (!hasMinimalCounts) {
                        tauPrime = tauPrime * 0.5;

                        // step 3 (again)
                        for (int k = 0; k < 100; k++) {

                            if (nextRecordIndex < 0) {
                                break;
                            }

                            // perform one step of Gillespie stochastic simulation algorithm
                            // update all arguments accordingly
                            double r0 = 0;
                            for (int c = 0; c < numReactionChannels; c++) {
                                r0 = r0 + reactionInt[c];
                            }

                            //find time to next reaction
                            double timeToReaction = -Math.log(MathUtils.nextDouble())/r0;

                            // if next reaction occurs after nextIntervalStartTime, record current compartment counts for next interval
                            while ((simulationTime + timeToReaction > nextIntervalStartTime) && nextRecordIndex >= 0) {
                                for (int s = 0; s < numSpecies; s++) {
                                    compartmentalModel.compartmentCounts.get(s).setParameterValue(nextRecordIndex, currentCounts[s]);
                                }
                                nextRecordIndex--;
                                nextIntervalStartTime = nextIntervalStartTime + intervalWidth;
                            }

                            int sampledReactionChannel = sampleReactionChannel(reactionInt, r0);

                            //update simulationTime and current compartment counts
                            //double previousSimulationTime = simulationTime;
                            simulationTime = simulationTime + timeToReaction;

                            currentCounts = compartmentalModel.introduceSecondPathogen(
                                    //previousSimulationTime,
                                    simulationTime,
                                    currentCounts);

                            for (int s = 0; s < numSpecies; s++) {
                                currentCounts[s] = currentCounts[s] + vMatrix[s][sampledReactionChannel];
                            }

                            reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
                        }
                        // end of 100 SSA steps

                    }

                    // Update before going back through loop
                    //reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
                    maxFiringTimes = getMaxFiringTimes(currentCounts, reactionInt);
                    critical.clear();
                    noncritical.clear();

                    for (int c = 0; c < numReactionChannels; c++) {
                        if (maxFiringTimes[c] < criticalNumber) {
                            critical.add(c);
                        } else {
                            noncritical.add(c);
                        }
                    }

                }

                // We can proceed and update the currentCounts and simulationTime and other parameters as necessary
                //double previousSimulationTime = simulationTime;
                simulationTime = simulationTime + tau;

                currentCounts = compartmentalModel.introduceSecondPathogen(
                        //previousSimulationTime,
                        simulationTime,
                        currentCounts);

                // if next reaction occurs after nextIntervalStartTime, record current compartment counts for next interval
                while ((simulationTime > nextIntervalStartTime) && nextRecordIndex >= 0) {
                    for (int s = 0; s < numSpecies; s++) {
                        compartmentalModel.compartmentCounts.get(s).setParameterValue(nextRecordIndex, currentCounts[s]);
                    }
                    nextRecordIndex--;
                    nextIntervalStartTime = nextIntervalStartTime + intervalWidth;
                }

                for (int s = 0; s < numSpecies; s++) {
                    currentCounts[s] = updatedCounts[s];
                }

            }
        }
    }

    protected double computeTauPrime(List<Integer> noncritical, double[] reactionInt, double[] e, double[] currentCounts) {

        double tauPrime = Double.POSITIVE_INFINITY;

        if (!noncritical.isEmpty()) {
            double[] muHat = new double[numSpecies];
            double[] sigmaSqHat = new double[numSpecies];

            for (int species = 0; species < numSpecies; species++) {
                double mu = 0.0;
                double s2 = 0.0;

                for (int idx = 0; idx < noncritical.size(); idx++) {
                    int reaction = noncritical.get(idx);
                    double vij = vMatrix[species][reaction];
                    mu = mu + vij * reactionInt[reaction];
                    s2 = s2 + (vij * vij) * reactionInt[reaction];
                }

                muHat[species] = mu;
                sigmaSqHat[species] = s2;
            }

            // numerator
            double[] bound = new double[numSpecies];

            for (int species = 0; species < numSpecies; species++) {
                bound[species] = Math.max(e[species] * currentCounts[species], 1.0);
            }

            for (int species = 0; species < numSpecies; species++) {
                // can't let denominator be 0
                double denom1 = Math.max(Math.abs(muHat[species]), 1e-16);
                double denom2 = Math.max(sigmaSqHat[species], 1e-16);
                double tau1 = bound[species] / denom1;
                double tau2 = (bound[species] * bound[species]) / denom2;
                double candidate = Math.min(tau1, tau2);

                if (candidate < tauPrime) {
                    tauPrime = candidate;
                }
            }
        }
        return tauPrime;
    }

    protected double[] getMaxFiringTimes(double[] currentCounts, double[] r){
        double[] returnVal = new double[numReactionChannels];
        for (int c = 0; c < numReactionChannels; c++) {
            returnVal[c] = Double.POSITIVE_INFINITY;
            for(int i = 0; i < numSpecies; i++) {
                if(vMatrix[i][c] < 0) {
                    double candidate = currentCounts[i] / Math.abs(vMatrix[i][c]);
                    if (r[c] > 0) {
                        returnVal[c] = Math.min(returnVal[c], candidate);
                    }
                }
            }
        }
        return returnVal;
    }

    protected double[] getTauLeapingPoissonIntensities(double[] reactionInt, double tau){
        double[] returnVal = new double[numReactionChannels];
        // for standard tau leaping
        for(int r = 0; r < numReactionChannels; r++) {
            returnVal[r] = reactionInt[r]*tau;
        }
        return returnVal;
    }

    protected double[] getPoissonIntensities(double[] currentCounts, double[] reactionInt, double tau) {
        return getTauLeapingPoissonIntensities(reactionInt, tau);
    }

}
