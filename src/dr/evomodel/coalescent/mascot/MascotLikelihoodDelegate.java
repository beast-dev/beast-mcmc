/*
 * MascotLikelihoodDelegate.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

/**
 * MASCOT-native likelihood delegate boundary.
 */
public interface MascotLikelihoodDelegate {

    double calculateLikelihood(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                               boolean checkProbabilities);

    double calculateLikelihoodAndDerivatives(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                             double[] gradientOut, double[] clockGradientOut,
                                             double[] ancestralStateScoresOut, boolean checkProbabilities);

    /**
     * Runs a taped forward pass and prepares the chosen stochastic mapper on it (see {@link MascotStochasticMapper}).
     * An {@link StochasticMapSampler#H_TRANSFORM} request falls back to {@link StochasticMapSampler#GRID} when the RK4
     * steps are too stiff for it.
     * {@code nodeHeights} holds every node's tree height by node number. Returns the forward pass's log likelihood.
     */
    double prepareStochasticMap(MascotPreparedInput prepared, double[] theta, double[] branchRates,
                                double[] nodeHeights, boolean checkProbabilities, StochasticMapSampler sampler);

    /** One joint draw from the most recent {@link #prepareStochasticMap} preparation. */
    MascotStochasticMap drawStochasticMap();

    final class Result {
        public final double logLikelihood;
        public final double[] gradient;
        public final double[] clockGradient;
        public final double[] rootProbabilities;
        public final int[] activeLineages;
        public final double[] ancestralStateScores;

        Result(double logLikelihood, double[] gradient, double[] clockGradient,
               double[] rootProbabilities, int[] activeLineages, double[] ancestralStateScores) {
            this.logLikelihood = logLikelihood;
            this.gradient = gradient;
            this.clockGradient = clockGradient;
            this.rootProbabilities = rootProbabilities;
            this.activeLineages = activeLineages;
            this.ancestralStateScores = ancestralStateScores;
        }
    }

    class NumericalException extends RuntimeException {
        public NumericalException(String message) {
            super(message);
        }
    }

    /** A pass would exceed its step budget (MascotStepRule); the likelihood there is -infinity, as for any numerical failure. */
    final class StepBudgetException extends NumericalException {
        public StepBudgetException(String message) {
            super(message);
        }
    }
}
