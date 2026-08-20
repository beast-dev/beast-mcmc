package dr.evomodel.epidemiology;

import cern.jet.random.Poisson;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import java.util.ArrayList;
import java.util.List;

public class TauLeapingSimulator extends StochasticSimulator {

    protected double epsilon;
    protected int criticalNumber;
    protected Parameter elapsedTime;

    // Count the number of times SSA is used instead of Tau or SAL
    // Number of separate times the simulator reverts to SSA
    protected long ssaCount = 0;
    // Number of individual SSA reaction steps performed
    protected long ssaStepCount = 0;
    // Number of times SSA is used due to invalid counts
    protected long invalidCountsCount = 0;
    // number of times SSA used due to tauprime too small
    protected long smallTauCount = 0;
    // number of times a tau step is accepted
    protected long successfulLeapCount = 0;

    public TauLeapingSimulator(CompartmentalModel compartmentalModel, double epsilon, int criticalNumber) {
        super(compartmentalModel);
        this.epsilon = epsilon;
        this.criticalNumber = criticalNumber;
    }

    public TauLeapingSimulator(CompartmentalModel compartmentalModel, double epsilon, int criticalNumber, Parameter elapsedTime) {
        this(compartmentalModel, epsilon, criticalNumber);
        this.elapsedTime = elapsedTime;
    }

    // Implements hybrid tau-leaping/SSA algorithm with step size selection as outlined by Cao et al. (2006)
    public void simulateTrajectory() {

        // To monitor how long it takes to simulate one trajectory
        long startTime = System.nanoTime();

        // set up time interval vector
        // duration for which we need to simulate trajectory
        double T = compartmentalModel.getOldestOrigin();

        // initialize values that we will have to keep track of
        SimulationState state = initializeSimulation(T);

        List<Integer> critical = new ArrayList<>(numReactionChannels);
        List<Integer> noncritical = new ArrayList<>(numReactionChannels);

        // Simulate until we have set compartment counts for all intervals
        // instead of while(simulationTime < T), keep track of nextRecordIndex
        while (state.nextRecordIndex >= 0) {
            runOneLeapingIteration(state, critical, noncritical);
        }

        printSimulationSummary();

        long endTime = System.nanoTime();
        double elapsedTimeInSeconds = (endTime - startTime) / 1e9;
        System.out.println("Elapsed time: " + elapsedTimeInSeconds + " seconds");
        elapsedTime.setParameterValue(0, elapsedTimeInSeconds);
    }

    /*
    private SimulationState initializeSimulation(double T) {

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
    */

    // A full run through Steps 1 - 6 of algorithm described in Cao et al. paper
    private void runOneLeapingIteration(SimulationState state, List<Integer> critical, List<Integer> noncritical) {

        while (true) {
            double[] e = computeEpsilonVector();

            // identifies critical and noncritical reactions
            stepOne(state, critical, noncritical);

            // computes tauPrime
            double tauPrime = stepTwo(state, noncritical, e);

            // proceeds if tauPrime is large enough, otherwise temporarily abandons tau-leaping,
            // performs 100 SSA steps, and then returns to Step 1
            tauPrime = stepThree(state, critical, noncritical, e, tauPrime);

            if (state.nextRecordIndex < 0) {
                return;
            }

            // Attempt Steps 4-6. May need to return to Step 3
            LeapOutcome outcome = attemptStepsFourThroughSix(state, critical, noncritical, tauPrime);

            if (outcome.status == LeapStatus.END_OF_TRAJECTORY_REACHED) {
                return;
            }

            if (outcome.status == LeapStatus.RESTART_FROM_STEP_ONE) {
                // Need to restart from Step 1, so skip the rest and go through loop again
                continue;
            }

            // We don't have to restart from Step 1, and we did not reach the end of the trajectory
            // So we have successfully exectued a tau leap
            executeLeap(state, outcome.leapResult);
            return;
        }
    }

    private double[] computeEpsilonVector() {
        int[] g = compartmentalModel.getHighestOrdersOfReactions();
        double[] e = new double[numSpecies];
        for (int s = 0; s < numSpecies; s++) {
            e[s] = epsilon / g[s];
        }
        return e;
    }

    private void stepOne(SimulationState state, List<Integer> critical, List<Integer> noncritical) {
        state.reactionInt = compartmentalModel.getReactionIntensities(state.currentCounts, state.simulationTime);
        // Maximum number of times that a reaction with a positive intensity can fire before
        // exhausting one of its reactants
        state.maxFiringTimes = getMaxFiringTimes(state.currentCounts, state.reactionInt);
        critical.clear();
        noncritical.clear();
        appendReactionClassification(state, critical, noncritical);
    }

    private void appendReactionClassification(SimulationState state, List<Integer> critical, List<Integer> noncritical) {
        for (int c = 0; c < numReactionChannels; c++) {
            if (state.maxFiringTimes[c] < criticalNumber) {
                critical.add(c);
            } else {
                noncritical.add(c);
            }
        }
    }

    private double stepTwo(SimulationState state, List<Integer> noncritical, double[] e) {
        return computeTauPrime(noncritical, state.reactionInt, e, state.currentCounts);
    }

    private double stepThree(SimulationState state, List<Integer> critical, List<Integer> noncritical,
                                                     double[] e, double tauPrime) {
        boolean tauPrimeTooSmall = true;

        while (tauPrimeTooSmall && state.nextRecordIndex >= 0) {

            double r0 = sumIntensities(state.reactionInt);

            if (!isTauPrimeTooSmall(tauPrime, r0, 10.0)) {
                tauPrimeTooSmall = false;
                // if tauPrime is not too small, skip the rest
                continue;
            }
            // if tauPrime is indeed too small, then
            // abandon tau-leaping temporarily and run 100 steps of Gillespie SSA algorithm
            ssaCount++;
            smallTauCount++;

            // DEBUGGING: Print a note when SSA is used
            //System.out.printf(
            //        getClass().getSimpleName() + " -> SSA (tauPrime too small): tauPrime=%f threshold=%f time=%f%n",
            //        tauPrime, 10.0 / r0, state.simulationTime
            //);

            runSSASteps(state, 100, true);

            if (state.nextRecordIndex < 0) {
                break;
            }

            // Do Step 1  and Step 2 again, after doing 100 SSA steps
            state.maxFiringTimes = getMaxFiringTimes(state.currentCounts, state.reactionInt);
            classifyReactions(state, critical, noncritical);
            tauPrime = computeTauPrime(noncritical, state.reactionInt, e, state.currentCounts);
        }
        return tauPrime;
    }

    /*
    private double sumIntensities(double[] reactionInt) {
        double r0 = 0;
        for (int c = 0; c < numReactionChannels; c++) {
            r0 = r0 + reactionInt[c];
        }
        return r0;
    }
    */


    // Checks if tauPrime is less than some small multiple (default 10) of 1/r0
    private boolean isTauPrimeTooSmall(double tauPrime, double r0, double multiple) {
        return tauPrime < multiple / r0;
    }

    // runs Gillespie stochastic simulation algorithm for certain number of steps
    // countSteps is set to true if we want to keep track of number of SSA steps taken
    private void runSSASteps(SimulationState state, int numSteps, boolean countSteps) {

        for (int k = 0; k < numSteps; k++) {

            if (state.nextRecordIndex < 0) {
                break;
            }

            // perform one step of Gillespie stochastic simulation algorithm

            double r0 = sumIntensities(state.reactionInt);

            //find time to next reaction
            double timeToReaction = -Math.log(MathUtils.nextDouble()) / r0;

            // if next reaction occurs after nextIntervalStartTime, record current compartment
            // counts for next interval
            recordCompartmentCountsUpTo(state, state.simulationTime + timeToReaction);

            int sampledReactionChannel = sampleReactionChannel(state.reactionInt, r0);

            if (countSteps) {
                // Add to total number of SSA steps taken
                ssaStepCount++;
            }

            // update simulationTime and current compartment counts
            state.simulationTime = state.simulationTime + timeToReaction;
            state.currentCounts = compartmentalModel.introduceSecondPathogen(state.simulationTime, state.currentCounts);

            for (int s = 0; s < numSpecies; s++) {
                state.currentCounts[s] = state.currentCounts[s] + vMatrix[s][sampledReactionChannel];
            }

            state.reactionInt = compartmentalModel.getReactionIntensities(state.currentCounts, state.simulationTime);

            // check if we need this
            //r0 = 0;
            //for (int c = 0; c < numReactionChannels; c++) {
            //    r0 = r0 + state.reactionInt[c];
            //}

        }
    }

    /*
    private void recordCompartmentCountsUpTo(SimulationState state, double candidateTime) {
        while (candidateTime > state.nextIntervalStartTime && state.nextRecordIndex >= 0) {
            for (int s = 0; s < numSpecies; s++) {
                compartmentalModel.compartmentCounts.get(s).setParameterValue(state.nextRecordIndex, state.currentCounts[s]);
            }
            state.nextRecordIndex--;
            state.nextIntervalStartTime = state.nextIntervalStartTime + intervalWidth;
        }
    }
    */

    private void classifyReactions(SimulationState state, List<Integer> critical, List<Integer> noncritical) {
        critical.clear();
        noncritical.clear();
        appendReactionClassification(state, critical, noncritical);
    }


    private LeapOutcome attemptStepsFourThroughSix(SimulationState state, List<Integer> critical, List<Integer> noncritical, double tauPrime) {

        // checks for negative counts or counts below a minimal threshold for certain species
        boolean hasMinimalCounts = false;

        // store candidate updated counts and candidate tau
        // will only be used if all counts meet minimal threshold
        double[] updatedCounts = new double[numSpecies];
        double tau = 0;

        while (!hasMinimalCounts && state.nextRecordIndex >= 0) {

            StepFourResult stepFourResult = stepFour(state, critical);
            StepFiveResult stepFiveResult = stepFive(state, critical, tauPrime,
                    stepFourResult.tauDoublePrime, stepFourResult.r0Critical);

            tau = stepFiveResult.tau;

            updatedCounts = compartmentalModel.getUpdatedCompartmentCounts(state.currentCounts, stepFiveResult.numFirings);
            hasMinimalCounts = compartmentalModel.hasMinimalCounts(updatedCounts);

            if (!hasMinimalCounts) {

                // we don't have minimal counts, so reduce tauPrime by half and return to Step 3
                // code below will check if tauPrime is large enough to continue to Step 4
                // or if we have to execute 100 SSA steps and then return to Step 1
                StepSixResult stepSixResult = stepSix(state, tauPrime);
                tauPrime = stepSixResult.newTauPrime;

                if (stepSixResult.returnToStepOne) {
                    LeapOutcome outcome = new LeapOutcome();
                    outcome.status = LeapStatus.RESTART_FROM_STEP_ONE;
                    // will return this outcome and skip the rest of the code
                    return outcome;
                }
                // If we don't have to return to Step 1, then that means
                // tauPrime is large enough,and we proceed to Step 4 again
                // we accomplish this by computing updates below and going through the loop again
            }

            // Update before going back through loop (or before exiting once hasMinimalCounts is true)
            state.maxFiringTimes = getMaxFiringTimes(state.currentCounts, state.reactionInt);
            classifyReactions(state, critical, noncritical);
        }

        LeapOutcome outcome = new LeapOutcome();
        if (state.nextRecordIndex < 0) {
            outcome.status = LeapStatus.END_OF_TRAJECTORY_REACHED;
        } else {
            outcome.status = LeapStatus.SUCCESS;
            LeapResult result = new LeapResult();
            result.tau = tau;
            result.updatedCounts = updatedCounts;
            outcome.leapResult = result;
        }
        return outcome;
    }


    private StepFourResult stepFour(SimulationState state, List<Integer> critical) {

        double r0Critical = 0.0;
        for (int index = 0; index < critical.size(); index++) {
            r0Critical += state.reactionInt[critical.get(index)];
        }

        // calculate a second option for tau
        double tauDoublePrime;

        if (r0Critical > 0.0) {
            double u = MathUtils.nextDouble();
            tauDoublePrime = -Math.log(u) / r0Critical;
        } else {
            tauDoublePrime = Double.POSITIVE_INFINITY;
        }

        StepFourResult result = new StepFourResult();
        result.r0Critical = r0Critical;
        result.tauDoublePrime = tauDoublePrime;
        return result;
    }


    private StepFiveResult stepFive(SimulationState state, List<Integer> critical,
                                            double tauPrime, double tauDoublePrime, double r0Critical) {
        // determine tau
        double tau = Math.min(tauPrime, tauDoublePrime);

        state.reactionInt = compartmentalModel.getReactionIntensities(state.currentCounts, state.simulationTime);
        double[] poissonIntensities = getPoissonIntensities(state.currentCounts, state.reactionInt, tau);

        double[] numFirings = new double[numReactionChannels];

        // Generate number of firings for each noncritical reaction as Poisson sample
        // For each critical reaction, no reaction occurs
        for (int j = 0; j < numReactionChannels; j++) {
            if (state.maxFiringTimes[j] < criticalNumber) {
                numFirings[j] = 0;
            } else {
                numFirings[j] = Poisson.staticNextInt(poissonIntensities[j]);
            }
        }

        // Step 5(b), one critical reaction occurs, other critical reactions are 0 (already done above)
        // noncritical reactions are poisson sample (already done above)
        // the critical reaction that does occur is sampled same way as SSA
        if (tauPrime >= tauDoublePrime) {
            // Step 5(b): exactly one critical reaction also fires, sampled the same way as SSA
            int jc = sampleCriticalReaction(state, critical, r0Critical);
            // exactly one critical reaction fires
            if (jc >= 0) {
                numFirings[jc] = 1;
            }
        }

        StepFiveResult result = new StepFiveResult();
        result.tau = tau;
        result.numFirings = numFirings;
        return result;
    }

    private int sampleCriticalReaction(SimulationState state, List<Integer> critical, double r0Critical) {
        if (r0Critical <= 0.0) {
            return -1;
        }
        double u = MathUtils.nextDouble() * r0Critical;
        double cumulative = 0.0;
        for (int idx = 0; idx < critical.size(); idx++) {
            int j = critical.get(idx);
            cumulative += state.reactionInt[j];
            if (u < cumulative) {
                return j;
            }
        }
        return -1;
    }

    private StepSixResult stepSix(SimulationState state, double tauPrime) {

        double newTauPrime = tauPrime * 0.5;
        ssaCount++;
        invalidCountsCount++;

        //System.out.printf(
        //        getClass().getSimpleName() + " -> SSA (invalid counts): tauPrime=%f newTauPrime=%f time=%f%n",
        //        newTauPrime * 2.0, newTauPrime, state.simulationTime
        //);

        double r0 = sumIntensities(state.reactionInt);

        StepSixResult result = new StepSixResult();
        result.newTauPrime = newTauPrime;

        // tauPrime has been halved, and now we return to Step 3
        // Check if it is too small, and if so, run 100 SSA steps and then return to Step 1
        if (isTauPrimeTooSmall(newTauPrime, r0, 10.0)) {
            runSSASteps(state, 100, true);
            result.returnToStepOne = true;
        } else {
            // if tauPrime is not too small, proceed to Step 4
            result.returnToStepOne = false;
        }
        return result;
    }

    private void executeLeap(SimulationState state, LeapResult leap) {

        state.simulationTime = state.simulationTime + leap.tau;
        successfulLeapCount++;

        state.currentCounts = compartmentalModel.introduceSecondPathogen(state.simulationTime, state.currentCounts);

        recordCompartmentCountsUpTo(state, state.simulationTime);

        for (int s = 0; s < numSpecies; s++) {
            state.currentCounts[s] = leap.updatedCounts[s];
        }
    }



    // Implements hybrid tau-leaping/SSA algorithm with step size selection as outlined by Cao et al. (2006)
    /*
    public void simulateTrajectory() {

        // To monitor how long it takes to simulate one trajectory
        long startTime = System.nanoTime();

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

                    // Add to count of times SSA is used
                    ssaCount++;
                    smallTauCount++;

                    // DEBUGGING: Print a note when SSA is used
                    System.out.printf(
                            getClass().getSimpleName() + " -> SSA (tauPrime too small): tauPrime=%f threshold=%f time=%f%n",
                            tauPrime,
                            10.0 / r0,
                            simulationTime
                    );

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

                        // Add to total number of SSA steps taken
                        ssaStepCount++;

                        //update simulationTime and current compartment counts
                        //double previousSimulationTime = simulationTime;
                        simulationTime = simulationTime + timeToReaction;

                        currentCounts = compartmentalModel.introduceSecondPathogen(
                                //previousSimulationTime,
                                simulationTime,
                                currentCounts);

                        // DEBUGGING

                        //double before = 0.0;
                        //for (int s = 0; s < numSpecies; s++) {
                        //    before += currentCounts[s];
                        //}

                        //System.out.println("SSA reaction = " + sampledReactionChannel);
                        //System.out.println("SSA total before = " + before);

                        // END DEBUGGING

                        for (int s = 0; s < numSpecies; s++) {
                            currentCounts[s] = currentCounts[s] + vMatrix[s][sampledReactionChannel];
                        }
                        reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
                        r0 = 0;
                        for (int c = 0; c < numReactionChannels; c++) {
                            r0 = r0 + reactionInt[c];
                        }

                        // DEBUGGING

                        //double after = 0.0;
                        //for (int s = 0; s < numSpecies; s++) {
                        //    after += currentCounts[s];
                        //}

                        //System.out.println("SSA total after = " + after);

                        // END DEBUGGING
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

                        // Add to total times SSA is used and total times used due to invalid counts
                        ssaCount++;
                        invalidCountsCount++;


                        // DEBUGGING: print a note when the counts are invalid
                        System.out.printf(
                                getClass().getSimpleName() + " -> SSA (invalid counts): tauPrime=%f newTauPrime=%f time=%f%n",
                                tauPrime * 2.0,
                                tauPrime,
                                simulationTime
                        );

                        double r0 = 0;
                        for (int c = 0; c < numReactionChannels; c++) {
                            r0 = r0 + reactionInt[c];
                        }

                        if (tauPrime < 10.0 / r0) {


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

                                // DEBUGGING
                                //double before = 0.0;
                                //for (int s = 0; s < numSpecies; s++) {
                                //    before += currentCounts[s];
                                //}

                                //System.out.println("SSA reaction = " + sampledReactionChannel);
                                //System.out.println("SSA total before = " + before);
                                // END DEBUGGING

                                for (int s = 0; s < numSpecies; s++) {
                                    currentCounts[s] = currentCounts[s] + vMatrix[s][sampledReactionChannel];
                                }

                                reactionInt = compartmentalModel.getReactionIntensities(currentCounts);
                            }
                            // end of 100 SSA steps
                        }else{
                            // Proceed to Step 4 (again)
                            // need to fill this in with code that does everything from Step 4 onward again
                        }

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
                successfulLeapCount++;

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

        System.out.printf(
                "%s SSA summary: Successful leaps: %d, fallbacks=%d, SSA steps=%d, " +
                        "small-tau fallbacks=%d, invalid-count fallbacks=%d%n",
                getClass().getSimpleName(),
                successfulLeapCount,
                ssaCount,
                ssaStepCount,
                smallTauCount,
                invalidCountsCount
        );
        long endTime = System.nanoTime();
        double elapsedTimeInSeconds = (endTime - startTime)/1e9;
        System.out.println("Elapsed time: " + elapsedTimeInSeconds + " seconds");
        elapsedTime.setParameterValue(0, elapsedTimeInSeconds);
    }
    */

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

    // get all SSA useage numbers
    public long getSSACount() {
        return ssaCount;
    }

    public long getSSAStepCount() {
        return ssaStepCount;
    }

    public long getSmallTauCount() {
        return smallTauCount;
    }

    public long getInvalidCountsCount() {
        return invalidCountsCount;
    }

    /*
    private static class SimulationState {

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
    */

    private static class LeapResult {
        // updated tau of successful tau leap
        double tau;
        // updated counts of successful tau leap
        double[] updatedCounts;
    }

    private enum LeapStatus { SUCCESS, RESTART_FROM_STEP_ONE, END_OF_TRAJECTORY_REACHED }

    private static class LeapOutcome {
        LeapStatus status;
        // update leap result with updated counts and tau if we have a success
        LeapResult leapResult;
    }

    // Keeps track of both things computed in Step 4
    private static class StepFourResult {
        // sum of reaction intensities for all critical reactions
        double r0Critical;
        // second candidate time leap
        double tauDoublePrime;
    }

    // keeps track of both things computed in Step 5
    private static class StepFiveResult {
        // which tau candidate to set as tau
        double tau;
        // number of firings for each reaction
        double[] numFirings;
    }

    // keeps track of possible results of Step 6
    private static class StepSixResult {
        // new tauPrime computed, then return to Step 3 and evaluate
        double newTauPrime;
        // return to Step 1 if the tau leap was successful and we need to continue simulation
        boolean returnToStepOne;
    }

    private void printSimulationSummary() {
        System.out.printf(
                "%s SSA summary: Successful leaps: %d, fallbacks=%d, SSA steps=%d, " +
                        "small-tau fallbacks=%d, invalid-count fallbacks=%d%n",
                getClass().getSimpleName(),
                successfulLeapCount,
                ssaCount,
                ssaStepCount,
                smallTauCount,
                invalidCountsCount
        );
    }


}