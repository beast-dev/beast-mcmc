package dr.app.beagle.evomodel.substmodel;

import dr.inference.markovjumps.*;
import dr.inference.model.Model;

/**
 * A class extension for implementing Markov chain-induced counting processes (markovjumps)
 * via uniformization in BEAST using BEAGLE
 * <p/>
 * This work is supported by NSF grant 0856099
 * <p/>
 * Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 * Journal of Mathematical Biology, 56, 391-412.
 * <p/>
 * Rodrigue N, Philippe H and Lartillot N (2006) Uniformization for sampling realizations of Markov processes:
 * applications to Bayesian implementations of codon substitution models. Bioinformatics, 24, 56-62.
 * <p/>
 * Hobolth A and Stone E (2009) Simulation from endpoint-conditioned, continuous-time Markov chains on a finite
 * state space, with applications to molecular evolution. Annals of Applied Statistics, 3, 1204-1231.
 *
 * @author Marc Suchard
 * @author Vladimir Minin
 */

public class UniformizedSubstitutionModel extends MarkovJumpsSubstitutionModel {

    public UniformizedSubstitutionModel(SubstitutionModel substModel) {
        this(substModel, MarkovJumpsType.COUNTS);
    }

    public UniformizedSubstitutionModel(SubstitutionModel substModel, MarkovJumpsType type) {
        this(substModel, type, 1);
    }

    public UniformizedSubstitutionModel(SubstitutionModel substModel, MarkovJumpsType type, int numSimulants) {
        super(substModel, type);
        this.numSimulants = numSimulants;
        updateSubordinator = true;
    }

    protected void setupStorage() {
        registration = new double[stateCount * stateCount];
        tmp = new double[stateCount * stateCount];
    }

    private void constructSubordinator() {
        substModel.getInfinitesimalMatrix(tmp);
        subordinator = new SubordinatedProcess(tmp, stateCount);
        updateSubordinator = false;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == substModel) {
            updateSubordinator = true;
        }
    }

    public void computeCondStatMarkovJumps(double time,
                                           double[] countMatrix) {

        throw new IllegalArgumentException("Not implemented for UniformizedSubstitutionModel");
    }

    public void computeCondStatMarkovJumps(double time,
                                           double[] transitionProbs,
                                           double[] countMatrix) {

        throw new IllegalArgumentException("Not implemented for UniformizedSubstitutionModel");
    }

    public void computeJointStatMarkovJumps(double time,
                                            double[] countMatrix) {

        throw new IllegalArgumentException("Not implemented for UniformizedSubstitutionModel");
    }

    public double computeCondStatMarkovJumps(int startingState,
                                             int endingState,
                                             double time) {

        substModel.getTransitionProbabilities(time, tmp);
        return computeCondStatMarkovJumps(startingState, endingState, time,
                tmp[startingState * stateCount + endingState]);
    }

    public double computeCondStatMarkovJumps(int startingState,
                                             int endingState,
                                             double time,
                                             double transitionProbability) {

        if (updateSubordinator) {
            constructSubordinator();
        }

        double total = 0;
        for (int i = 0; i < numSimulants; i++) {
            StateHistory history = UniformizedStateHistory.simulateConditionalOnEndingState(
                    0.0,
                    startingState,
                    time,
                    endingState,
                    transitionProbability,
                    stateCount,
                    subordinator
            );
            if (type == MarkovJumpsType.COUNTS) {
                total += history.getTotalRegisteredCounts(registration);
            } else {
                total += history.getTotalReward(registration);
            }
        }
        return total / (double) numSimulants;
    }

    private final int numSimulants;
    private boolean updateSubordinator;
    private SubordinatedProcess subordinator;

    private double[] tmp;
}
